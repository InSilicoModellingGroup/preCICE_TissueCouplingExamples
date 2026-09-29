#ifndef BDM_ADAPTER_H_
#define BDM_ADAPTER_H_

/**
 * @file bdm_adapter.h
 * @brief BioDynaMo-preCICE Adapter Library
 *
 * This header-only library provides coupling between BioDynaMo agent-based
 * simulations and external solvers (e.g., OpenFOAM) via preCICE.
 *
 * Usage:
 * @code
 *   #include "bdm_adapter.h"
 *
 *   // Define your cell type with required interface
 *   class MyCell : public bdm::Cell {
 *       void SetFluidVelocity(const bdm::Real3& v);
 *       const bdm::Real3& GetFluidVelocity() const;
 *       void SetPressure(double p);
 *       double GetPressure() const;
 *       void SetPressureGradient(const std::array<double, 3>& g);
 *       std::array<double, 3> GetPressureGradient() const;
 *       bool IsStationary() const;
 *   };
 *
 *   // Create adapter with your cell type
 *   bdm_adapter::PreciceAdapter<MyCell> adapter(params);
 *   adapter.Initialize();
 *   while (adapter.IsCouplingOngoing()) {
 *       adapter.FetchInputData(dt);
 *       // ... simulation step ...
 *       adapter.Advance(dt);
 *   }
 *   adapter.Finalize();
 * @endcode
 */

// ============================================================================
// Standard Library Includes
// ============================================================================
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cctype>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <unordered_map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include <random>

// ============================================================================
// External Dependencies
// ============================================================================
#include "biodynamo.h"
#include "precice/precice.hpp"

// ============================================================================
// Configuration Headers
// ============================================================================
#include "coupling_configuration.h"
#include "parameters.h"

// ============================================================================
// Utility Headers
// ============================================================================
#include "mesh_mapping.h"
#include "boundary_geometry.h"
#include "cell_interface.h"

namespace bdm_adapter {

namespace detail {
template <typename T, typename = void>
struct has_set_velocity : std::false_type {};

template <typename T>
struct has_set_velocity<T, std::void_t<decltype(std::declval<T*>()->SetVelocity(std::declval<bdm::Real3>()))>>
    : std::true_type {};

template <typename T, typename = void>
struct has_get_diameter : std::false_type {};

template <typename T>
struct has_get_diameter<T, std::void_t<decltype(std::declval<const T&>().GetDiameter())>>
    : std::true_type {};
}  // namespace detail

// ============================================================================
// Utility structs for mesh and bounds
// ============================================================================

/**
 * @brief Main preCICE adapter for BioDynaMo simulations
 *
 * This templated class handles all preCICE interactions:
 * - Participant initialization and lifecycle
 * - Mesh setup and data exchange
 * - Reading/writing field data (velocity, pressure, etc.)
 *
 * @tparam CellType The BioDynaMo cell type (must have required interface methods)
 * @tparam Config   The coupling configuration type
 */
template<typename CellType, typename Config = DefaultCouplingConfig>
class PreciceAdapter {
public:
    using ConfigType = Config;

    /**
     * @brief Construct adapter with given parameters
     * @param params Adapter configuration parameters
     */
    explicit PreciceAdapter(const AdapterParameters& params = GetDefaultAdapterParameters())
        : params_(params)
        , config_()
    {
        SetupProfiling();
    }

    ~PreciceAdapter() {
        if (participant_) {
            Finalize();
        }
    }

    void ReadFluidAccelerationData(double relative_time) {
        if (!fluid_acceleration_available_) return;

        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;

        // Use shared read-buffer built by FetchInputData; fallback to rebuild if called standalone.
        if (!read_buf_valid_) RebuildReadBuffer();
        const auto& coords = read_coords_buf_;
        const auto& uids   = read_uid_buf_;

        if (uids.empty()) return;

        std::vector<double> values(uids.size() * 3, 0.0);
        auto coords_span = precice::span<const double>(coords.data(), coords.size());
        auto values_span = precice::span<double>(values.data(), values.size());

        try {
            participant_->mapAndReadData(mesh_name, params_.precice.fluid_acceleration_data_name,
                                        coords_span, relative_time, values_span);

            size_t i = 0;
            for (const auto& uid : uids) {
                auto* cell = dynamic_cast<CellType*>(rm->GetAgent(uid));
                if (cell) {
                    bdm::Real3 accel = {values[i*3], values[i*3+1], values[i*3+2]};
                    if constexpr (has_fluid_acceleration_v<CellType>) {
                        cell->SetFluidAcceleration(accel);
                    }
                }
                ++i;
            }
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter", "Failed to read DUcDt: ", e.what());
        }
    }

    void ReadFluidVorticityData(double relative_time) {
        if (!fluid_vorticity_available_) return;

        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;

        // Use shared read-buffer built by FetchInputData; fallback to rebuild if called standalone.
        if (!read_buf_valid_) RebuildReadBuffer();
        const auto& coords = read_coords_buf_;
        const auto& uids   = read_uid_buf_;

        if (uids.empty()) return;

        std::vector<double> values(uids.size() * 3, 0.0);
        auto coords_span = precice::span<const double>(coords.data(), coords.size());
        auto values_span = precice::span<double>(values.data(), values.size());

        try {
            participant_->mapAndReadData(mesh_name, params_.precice.fluid_vorticity_data_name,
                                        coords_span, relative_time, values_span);

            size_t i = 0;
            for (const auto& uid : uids) {
                auto* cell = dynamic_cast<CellType*>(rm->GetAgent(uid));
                if (cell) {
                    bdm::Real3 vorticity = {values[i*3], values[i*3+1], values[i*3+2]};
                    if constexpr (has_fluid_vorticity_v<CellType>) {
                        cell->SetFluidVorticity(vorticity);
                    }
                }
                ++i;
            }
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter", "Failed to read curlUc: ", e.what());
        }
    }

    // =========================================================================
    // Lifecycle Management
    // =========================================================================

    /**
     * @brief Initialize preCICE participant and meshes
     *
     * Must be called before any data exchange. Sets up:
     * - preCICE participant
     * - Volume meshes
     * - Wall meshes
     * - Data field registration
     */
    void Initialize() {
        if (is_initialized_) return;

        // Create preCICE participant
        participant_ = std::make_unique<precice::Participant>(
            params_.precice.participant_name,
            params_.precice.config_file,
            0, 1  // rank, size (single process for now)
        );

        // Set mesh access regions BEFORE initialize() - required for receiving mesh data
        for (const auto& mesh_name : params_.mesh.default_volume_meshes) {
            try {
                int dimensions = participant_->getMeshDimensions(mesh_name);
                SetWideAccessRegion(mesh_name, dimensions);
            } catch (const std::exception& e) {
                bdm::Log::Warning("PreciceAdapter", "Could not set access region for '", mesh_name, "': ", e.what());
            }
        }
        for (const auto& mesh_name : params_.mesh.default_wall_meshes) {
            try {
                int dimensions = participant_->getMeshDimensions(mesh_name);
                SetWideAccessRegion(mesh_name, dimensions);
            } catch (const std::exception& e) {
                bdm::Log::Warning("PreciceAdapter", "Could not set access region for '", mesh_name, "': ", e.what());
            }
        }
        if (!params_.mesh.boundary_triangle_mesh.empty()) {
            try {
                const int dimensions = participant_->getMeshDimensions(
                    params_.mesh.boundary_triangle_mesh);
                SetWideAccessRegion(params_.mesh.boundary_triangle_mesh, dimensions);
            } catch (const std::exception& e) {
                bdm::Log::Warning("PreciceAdapter",
                    "Could not set access region for exact boundary mesh '",
                    params_.mesh.boundary_triangle_mesh, "': ", e.what());
            }
        }

        // Set access regions for inlet/outlet meshes if configured
        for (const auto& mesh_name : params_.mesh.default_inlet_meshes) {
            try {
                int dimensions = participant_->getMeshDimensions(mesh_name);
                SetWideAccessRegion(mesh_name, dimensions);
            } catch (const std::exception& e) {
                bdm::Log::Warning("PreciceAdapter", "Could not set access region for inlet mesh '", mesh_name, "': ", e.what());
            }
        }
        for (const auto& mesh_name : params_.mesh.default_outlet_meshes) {
            try {
                int dimensions = participant_->getMeshDimensions(mesh_name);
                SetWideAccessRegion(mesh_name, dimensions);
            } catch (const std::exception& e) {
                bdm::Log::Warning("PreciceAdapter", "Could not set access region for outlet mesh '", mesh_name, "': ", e.what());
            }
        }

        // Initialize preCICE (must be called after setMeshAccessRegion but before accessing mesh data).
        //
        // NOTE: preCICE may emit "Mapping distance not available due to empty partition" here.
        // This is EXPECTED for just-in-time (JIT) mapping (receive-mesh with api-access="true").
        // The JIT mapper has no pre-registered target/query coordinates at initialize() time;
        // target coordinates are supplied per-call via mapAndReadData().  The first real JIT
        // read happens in ReadVelocityData() / FetchInputData(), at which point agents already
        // exist at inlet-cell-centre positions and the coordinate set is non-empty.
        participant_->initialize();

        // Get timestep after initialization
        dt_ = participant_->getMaxTimeStepSize();

        // Now initialize meshes (after preCICE initialize)
        InitializeVolumeMeshes();
        InitializeWallMeshes();
        InitializeInletOutletMeshes();

        // Setup data availability flags
        SetupDataAvailability();

        is_initialized_ = true;
        bdm::Log::Info("PreciceAdapter", "Initialized successfully. dt = ", dt_);
    }

    /**
     * @brief Finalize preCICE coupling
     */
    void Finalize() {
        if (participant_) {
            participant_->finalize();
            participant_.reset();
        }
        is_initialized_ = false;
    }

    /**
     * @brief Advance the coupling by one timestep
     * @param dt Timestep size
     */
    void Advance(double dt) {
        if (!participant_ || !participant_->isCouplingOngoing()) return;

        participant_->advance(dt);
        dt_ = participant_->getMaxTimeStepSize();
        if (participant_->isTimeWindowComplete()) {
            coupling_iteration_ = 0;
        } else {
            ++coupling_iteration_;
        }
    }

    /**
     * @brief Check if coupling is still ongoing
     */
    bool IsCouplingOngoing() const {
        return participant_ && participant_->isCouplingOngoing();
    }

    /**
     * @brief Get current timestep size
     */
    double GetTimeStepSize() const { return dt_; }

    /**
     * @brief Get current coupling iteration
     */
    int GetCouplingIteration() const { return coupling_iteration_; }

    /**
     * @brief Check if checkpoint reading is required (implicit coupling)
     */
    bool RequiresReadingCheckpoint() const {
        return participant_ && participant_->requiresReadingCheckpoint();
    }

    /**
     * @brief Check if checkpoint writing is required (implicit coupling)
     */
    bool RequiresWritingCheckpoint() const {
        return participant_ && participant_->requiresWritingCheckpoint();
    }

    /**
     * @brief Check if current time window is complete
     */
    bool IsTimeWindowComplete() const {
        return participant_ && participant_->isTimeWindowComplete();
    }

    /**
     * @brief Get domain bounds from the received FVMmesh (cell centers)
     *
     * Returns the bounding box computed from the mesh vertices received
     * from OpenFOAM via preCICE. Note: FVMmesh vertices are cell centers,
     * which are slightly inside the actual domain walls.
     *
     * For actual wall positions, use GetActualDomainBounds() which uses
     * WallSurfaceMesh face centers.
     *
     * @param mesh_name Name of the mesh to get bounds from (default: primary volume mesh)
     * @return BoundingBox with min/max corners from the mesh cell centers
     */
    BoundingBox GetDomainBoundsFromMesh(const std::string& mesh_name = "") const {
        const std::string& target_mesh = mesh_name.empty() ?
            params_.mesh.primary_volume_mesh : mesh_name;

        auto it = volume_meshes_.find(target_mesh);
        if (it != volume_meshes_.end() && it->second.is_initialized) {
            bdm::Log::Info("PreciceAdapter", "Cell center bounds from '", target_mesh, "': [",
                          it->second.bounds.min[0], ",", it->second.bounds.min[1], ",", it->second.bounds.min[2],
                          "] to [",
                          it->second.bounds.max[0], ",", it->second.bounds.max[1], ",", it->second.bounds.max[2], "]");
            return it->second.bounds;
        }

        bdm::Log::Warning("PreciceAdapter", "Mesh '", target_mesh, "' not found, returning default bounds");
        return BoundingBox{};
    }

    /**
     * @brief Utility to check if a BoundingBox has valid values
     */
    bool IsValidBounds(const BoundingBox& b) const {
        return !std::isnan(b.min[0]) && !std::isnan(b.max[0]);
    }

    /**
     * @brief Read OpenFOAM polyMesh points and compute exact domain bounds
     */
    BoundingBox ReadPolyMeshBounds(const std::string& foam_case_dir) const {
        BoundingBox bounds;

        const std::string points_path = foam_case_dir + "/constant/polyMesh/points";
        std::ifstream points_file(points_path);
        if (!points_file) {
            bdm::Log::Warning("PreciceAdapter", "Could not open polyMesh points at '", points_path, "'");
            return bounds;
        }

        std::string line;
        bool started = false;
        while (std::getline(points_file, line)) {
            if (!started) {
                if (line.size() == 1 && line[0] == '(') {
                    started = true;
                }
                continue;
            }

            if (!line.empty() && line[0] == ')') {
                break;
            }

            for (char& c : line) {
                if (c == '(' || c == ')') c = ' ';
            }
            std::stringstream ss(line);
            double x, y, z;
            if (!(ss >> x >> y >> z)) {
                continue;
            }

            if (!IsValidBounds(bounds)) {
                bounds.min = {x, y, z};
                bounds.max = bounds.min;
            } else {
                bounds.min[0] = std::min(bounds.min[0], x);
                bounds.min[1] = std::min(bounds.min[1], y);
                bounds.min[2] = std::min(bounds.min[2], z);
                bounds.max[0] = std::max(bounds.max[0], x);
                bounds.max[1] = std::max(bounds.max[1], y);
                bounds.max[2] = std::max(bounds.max[2], z);
            }
        }

        if (IsValidBounds(bounds)) {
            bdm::Log::Info("PreciceAdapter", "polyMesh bounds: [", bounds.min[0], ",", bounds.min[1], ",", bounds.min[2],
                           "] -> [", bounds.max[0], ",", bounds.max[1], ",", bounds.max[2], "]");
        } else {
            bdm::Log::Warning("PreciceAdapter", "polyMesh bounds not initialized (no points parsed)");
        }
        return bounds;
    }

    /**
     * @brief Get actual domain bounds from WallSurfaceMesh
     *
     * Returns the bounding box computed from wall face centers, which represent
     * the actual domain boundaries (not cell centers like FVMmesh).
     * This is more accurate for wall collision detection and BioDynaMo domain setup.
     *
     * @return BoundingBox with min/max corners from wall face centers
     */
    BoundingBox GetActualDomainBounds() const {
        BoundingBox bounds;
        bool first = true;

        for (const auto& [mesh_name, mesh] : wall_meshes_) {
            if (!mesh.is_initialized || mesh.geometry.empty()) continue;

            // Wall geometry contains face centers (3 components per face)
            size_t num_faces = mesh.geometry.size() / 3;
            for (size_t i = 0; i < num_faces; ++i) {
                double x = mesh.geometry[i * 3];
                double y = mesh.geometry[i * 3 + 1];
                double z = mesh.geometry[i * 3 + 2];

                if (first) {
                    bounds.min = {x, y, z};
                    bounds.max = {x, y, z};
                    first = false;
                } else {
                    bounds.min[0] = std::min(bounds.min[0], x);
                    bounds.min[1] = std::min(bounds.min[1], y);
                    bounds.min[2] = std::min(bounds.min[2], z);
                    bounds.max[0] = std::max(bounds.max[0], x);
                    bounds.max[1] = std::max(bounds.max[1], y);
                    bounds.max[2] = std::max(bounds.max[2], z);
                }
            }
        }

        if (!first) {
            bdm::Log::Info("PreciceAdapter", "Actual domain bounds from WallSurfaceMesh: [",
                          bounds.min[0], ",", bounds.min[1], ",", bounds.min[2],
                          "] to [",
                          bounds.max[0], ",", bounds.max[1], ",", bounds.max[2], "]");
        } else {
            bdm::Log::Warning("PreciceAdapter", "No wall mesh data available, returning default bounds");
        }

        return bounds;
    }

    /**
     * @brief Check if wall mesh data is available for actual domain bounds
     */
    bool HasWallMeshGeometry() const {
        for (const auto& [name, mesh] : wall_meshes_) {
            if (mesh.is_initialized && !mesh.geometry.empty()) {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Check if mesh has valid domain bounds
     */
    bool HasValidDomainBounds(const std::string& mesh_name = "") const {
        const std::string& target_mesh = mesh_name.empty() ?
            params_.mesh.primary_volume_mesh : mesh_name;

        auto it = volume_meshes_.find(target_mesh);
        return it != volume_meshes_.end() && it->second.is_initialized &&
               it->second.NumVertices() > 0;
    }

    /**
     * @brief Compute global domain bounds from all available meshes
     *
     * Aggregates bounds from volume meshes (cell centers) and wall meshes
     * (face centers) to derive a unified domain bounding box. Updates
     * `domain_bounds_` and logs the result for diagnostics.
     */
    void ComputeAndSetGlobalDomainBounds() {
        BoundingBox global{};
        bool first = true;

        auto accumulate_bounds = [&](const BoundingBox& b, const std::string& source) {
            bdm::Log::Info("PreciceAdapter", "  Bounds from '", source, "': [",
                          b.min[0], ",", b.min[1], ",", b.min[2], "] -> [",
                          b.max[0], ",", b.max[1], ",", b.max[2], "]");
            if (first) {
                global = b;
                first = false;
                return;
            }
            global.min[0] = std::min(global.min[0], b.min[0]);
            global.min[1] = std::min(global.min[1], b.min[1]);
            global.min[2] = std::min(global.min[2], b.min[2]);
            global.max[0] = std::max(global.max[0], b.max[0]);
            global.max[1] = std::max(global.max[1], b.max[1]);
            global.max[2] = std::max(global.max[2], b.max[2]);
        };

        bdm::Log::Info("PreciceAdapter", "Computing global domain bounds from meshes...");

        // First, try exact bounds from polyMesh (OpenFOAM geometry)
        BoundingBox poly_bounds = ReadPolyMeshBounds(params_.mesh.foam_case_directory);
        if (IsValidBounds(poly_bounds)) {
            accumulate_bounds(poly_bounds, "polyMesh");
        }

        for (const auto& [name, mesh] : volume_meshes_) {
            if (mesh.is_initialized && mesh.NumVertices() > 0) {
                bdm::Log::Info("PreciceAdapter", "  Volume mesh '", name, "' has ", mesh.NumVertices(), " vertices");
                // Log first few vertices for debugging
                size_t num_to_log = std::min(size_t(5), mesh.NumVertices());
                for (size_t i = 0; i < num_to_log; ++i) {
                    auto v = mesh.GetVertex(i);
                    bdm::Log::Info("PreciceAdapter", "    Vertex ", i, ": (", v[0], ", ", v[1], ", ", v[2], ")");
                }
                accumulate_bounds(mesh.bounds, name);
            }
        }

        for (const auto& [name, mesh] : wall_meshes_) {
            if (mesh.is_initialized && !mesh.geometry.empty()) {
                BoundingBox b;
                b.min = b.max = {mesh.geometry[0], mesh.geometry[1], mesh.geometry[2]};
                size_t n = mesh.geometry.size() / 3;
                bdm::Log::Info("PreciceAdapter", "  Wall mesh '", name, "' has ", n, " faces");
                // Log first few face centers for debugging
                size_t num_to_log = std::min(size_t(5), n);
                for (size_t i = 0; i < num_to_log; ++i) {
                    bdm::Log::Info("PreciceAdapter", "    Face ", i, " center: (",
                                  mesh.geometry[3*i], ", ", mesh.geometry[3*i+1], ", ", mesh.geometry[3*i+2], ")");
                }
                for (size_t i = 0; i < n; ++i) {
                    double x = mesh.geometry[3 * i];
                    double y = mesh.geometry[3 * i + 1];
                    double z = mesh.geometry[3 * i + 2];
                    b.min[0] = std::min(b.min[0], x);
                    b.min[1] = std::min(b.min[1], y);
                    b.min[2] = std::min(b.min[2], z);
                    b.max[0] = std::max(b.max[0], x);
                    b.max[1] = std::max(b.max[1], y);
                    b.max[2] = std::max(b.max[2], z);
                }
                accumulate_bounds(b, name);
            }
        }

        if (!first) {
            domain_bounds_ = global;
            auto size = global.Size();
            bdm::Log::Info("PreciceAdapter", "Global domain bounds set to: [", global.min[0],
                           ",", global.min[1], ",", global.min[2], "] -> [",
                           global.max[0], ",", global.max[1], ",", global.max[2], "]");
            bdm::Log::Info("PreciceAdapter", "  Domain extents: X=", size[0], ", Y=", size[1], ", Z=", size[2]);
        } else {
            bdm::Log::Warning("PreciceAdapter", "ComputeAndSetGlobalDomainBounds: no meshes available to determine bounds");
        }
    }

    // =========================================================================
    // Data Exchange
    // =========================================================================

    /**
     * @brief Fetch all input data from preCICE
     *
     * Reads velocity, pressure, and optionally temperature data
     * from preCICE and applies it to all dynamic agents.
     *
     * @param relative_time Relative read time within timestep
     */
    void FetchInputData(double relative_time = 0.0) {
        if (!is_initialized_) {
            Initialize();
        }

        if (!participant_) {
            ApplyCachedFieldDataToAgents();
            return;
        }

        if (!participant_->isCouplingOngoing()) {
            ApplyCachedFieldDataToAgents();
            return;
        }

        // preCICE 3 `mapAndReadData` takes a relative_read_time in
        // [0, getMaxTimeStepSize()].  Semantics:
        //   • 0        → read data at the START of the current window.
        //                For explicit coupling this is the MOST RECENT
        //                fully-completed write from the producer.  On
        //                window 1 it returns the initialize="true"
        //                buffer; on window k>1 it returns window k-1's
        //                end data.  Always usable, never extrapolated.
        //   • dt       → read data at the END of the current window.
        //                Requires the producer to have already written
        //                this window's data — NOT TRUE in parallel-
        //                explicit one-way coupling, where both
        //                participants advance concurrently and the
        //                producer only commits its window k value when
        //                advance() completes.  Returns zeros or a
        //                stale/extrapolated value in that case.
        // We therefore now interpret any NON-NEGATIVE relative_time
        // literally (0 is a legitimate, preferred value for explicit
        // couplings).  Legacy "use default dt_" behaviour is kept for
        // callers that pass a negative sentinel.
        double read_time;
        if (relative_time < 0.0) {
            read_time = dt_ > 0.0 ? dt_ : participant_->getMaxTimeStepSize();
            dt_ = read_time;
        } else {
            read_time = relative_time;  // 0 == start of window (explicit default)
            // keep dt_ up to date for any downstream consumers that read it
            if (dt_ <= 0.0) {
                dt_ = participant_->getMaxTimeStepSize();
            }
        }

        ++fetch_input_call_count_;
        const bool should_read_fields = ShouldReadFieldData();
        const bool should_read_walls = ShouldReadWallData();

        if (!should_read_fields && !field_freeze_notice_emitted_) {
            bdm::Log::Info("PreciceAdapter",
                           "Skipping CFD field read and reusing previously sampled values"
                           " (coupling.field_read_interval=",
                           params_.coupling.field_read_interval,
                           ", fetch_call=", fetch_input_call_count_, ")");
            field_freeze_notice_emitted_ = true;
        }
        if (!should_read_walls && !wall_freeze_notice_emitted_) {
            bdm::Log::Info("PreciceAdapter",
                           "Skipping wall mesh read and reusing cached wall geometry"
                           " (coupling.static_mesh=",
                           params_.coupling.static_mesh ? "true" : "false",
                           ", coupling.wall_read_interval=",
                           params_.coupling.wall_read_interval,
                           ", fetch_call=", fetch_input_call_count_, ")");
            wall_freeze_notice_emitted_ = true;
        }

        // Read field data using runtime flags from params.csv.
        // This keeps one compiled cells participant usable across multiple OpenFOAM participants.
        if (should_read_fields) {
            // Build shared coord/uid buffer ONCE for all Read*Data calls below.
            // Eliminates N-1 redundant ForEachAgent passes when multiple fields
            // are active (e.g. velocity + pressure + TKE + ... in aneurysm).
            RebuildReadBuffer();

            if (params_.coupling.enable_velocity) {
                ReadVelocityData(read_time);
            }
            ReadFluidAccelerationData(read_time);
            ReadFluidVorticityData(read_time);
            if (params_.coupling.enable_pressure) {
                ReadPressureData(read_time);
                if (params_.coupling.enable_pressure_gradient) {
                    ReadPressureGradientData(read_time);
                }
            }
            if (params_.coupling.enable_temperature) {
                ReadTemperatureData(read_time);
            }
            if (params_.coupling.enable_oxygen) {
                ReadOxygenData(read_time);
            }
            if (params_.coupling.enable_h2o2) {
                ReadH2O2Data(read_time);
            }
            if (params_.coupling.enable_no2) {
                ReadNO2Data(read_time);
            }

            // Read turbulence fields if enabled
            if (params_.coupling.enable_turbulence) {
                ReadTKEData(read_time);
                ReadEpsilonData(read_time);
                ReadOmegaData(read_time);
                ReadNutData(read_time);
            }

            // Invalidate shared buffer so it cannot be reused outside FetchInputData.
            read_buf_valid_ = false;

            field_data_read_once_ = true;

            if (ShouldUseCachedFieldResampling()) {
                CacheCurrentFieldSnapshot(read_time);
            }
        } else {
            ApplyCachedFieldDataToAgents();
        }

        // Read wall mesh data
        if (should_read_walls) {
            for (auto& [mesh_name, mesh] : wall_meshes_) {
                if (mesh.is_initialized) {
                    ReadWallMeshData(mesh_name, mesh, read_time);
                }
            }
            wall_data_read_once_ = true;
        }

        if (!domain_bounds_computed_ || should_read_walls) {
            ComputeAndSetGlobalDomainBounds();
            domain_bounds_computed_ = true;
        }
    }

    /**
     * @brief Cache CFD fields on mesh vertices for later standalone reuse.
     *
     * This snapshot enables repeated BioDynaMo-only iterations with moving
     * agents while the CFD field remains frozen.
     */
    void EnsureCachedFieldSnapshot(double relative_time = 0.0) {
        if (!participant_) {
            return;
        }
        CacheCurrentFieldSnapshot(relative_time);
    }

    /**
     * @brief Re-sample cached CFD fields at current agent positions.
     *
     * Use this during standalone continuation (after preCICE finalization)
     * or when field reads are intentionally throttled.
     */
    void RefreshCachedFieldData() {
        ApplyCachedFieldDataToAgents();
    }

    /**
     * @brief Map the current carrier velocity to every existing dynamic agent.
     *
     * This is intentionally a live JIT read while coupling is active. Newly
     * injected agents did not exist during FetchInputData() at the beginning
     * of the window and therefore cannot be initialized by merely reapplying
     * an optional cache (the cache is normally absent when fields are read
     * every window). Without this read, a new massless tracer starts with
     * zero fluid velocity and remains at its inlet owner-cell center.
     */
    void RefreshCurrentVelocityData(double relative_time = 0.0) {
        if (participant_ && participant_->isCouplingOngoing() &&
            params_.coupling.enable_velocity && velocity_available_) {
            RebuildReadBuffer();
            ReadVelocityData(relative_time);
            read_buf_valid_ = false;
            return;
        }

        // After coupling finalization, the only available carrier field is
        // the deliberately frozen mesh snapshot.
        ApplyCachedFieldDataToAgents();
    }

    /**
     * @brief Read only wall mesh geometry data without reading velocity/pressure
     * Useful for reading wall data before agents are created
     */
    void ReadInitialWallData() {
        if (!is_initialized_) return;

        if (params_.coupling.static_mesh && wall_data_read_once_) {
            return;
        }

        // Read at the START of the coupling window. In parallel-explicit this
        // reliably returns initialized wall data, while reading at dt can hit
        // not-yet-written end-of-window values and silently return zeros.
        double read_time = 0.0;
        bool read_any_wall = false;
        for (auto& [mesh_name, mesh] : wall_meshes_) {
            if (mesh.is_initialized) {
                ReadWallMeshData(mesh_name, mesh, read_time);
                read_any_wall = true;
            }
        }
        if (read_any_wall) {
            wall_data_read_once_ = true;
        }
        // Re-read exact static triangles at the same reliable initialization
        // instant used for the legacy wall fields.
        InitializeBoundaryTriangleMesh();
        ComputeAndSetGlobalDomainBounds();
        domain_bounds_computed_ = true;
    }
    /**
     * @brief Fetch fluid material properties provided by the CFD participant
     * @return true if both density and viscosity were received
     */
    bool FetchFluidProperties(double relative_time, double* density_out, double* mu_out) {
        if (!is_initialized_) {
            Initialize();
        }
        if (!participant_) {
            return false;
        }
        return ReadFluidPropertiesInternal(relative_time, density_out, mu_out);
    }

    /**
     * @brief Send output data to preCICE (for two-way coupling)
     *
     * Tissue feedback pathway: optional agent positions and scalar uptake sinks.
     *
     * Only active when two_way_coupling is enabled in params.
     */
    void SendOutputData() {
        if (!params_.coupling.two_way_coupling) {
            return;  // One-way coupling - no data to send
        }

        if (!is_initialized_ || !participant_ || !participant_->isCouplingOngoing()) {
            return;
        }

        // #region agent log
        {
            std::ostringstream payload;
            payload << "{"
                    << "\"iteration\":" << coupling_iteration_ << ","
                    << "\"twoWay\":" << (params_.coupling.two_way_coupling ? "true" : "false") << ","
                    << "\"oxygenSinkEnabled\":" << (params_.coupling.enable_oxygen_sink ? "true" : "false") << ","
                    << "\"oxygenSinkDataNameEnabled\":"
                    << (IsDataNameEnabled(params_.precice.oxygen_sink_data_name) ? "true" : "false")
                    << ",\"no2SinkEnabled\":"
                    << (params_.coupling.enable_no2_sink ? "true" : "false")
                    << ",\"no2SinkDataNameEnabled\":"
                    << (IsDataNameEnabled(params_.precice.no2_sink_data_name) ? "true" : "false")
                    << "}";
            DebugNdjsonLog("baseline", "H2", "bdm_adapter.h:SendOutputData",
                           "send_output_gate", payload.str());
        }
        // #endregion

        // OpenFOAM pattern: Write positions and source terms
        if (IsDataNameEnabled(params_.precice.position_data_name)) {
            WritePositionData();
        }

        WriteOxygenSinkData();
        WriteH2O2SinkData();
        WriteNO2SinkData();
    }

    /**
     * @brief Reset source terms (OpenFOAM pattern: resetSourceTerms)
     *
     * Call at the beginning of each timestep to clear accumulated sources.
     */
    void ResetSourceTerms() {
        oxygen_sink_source_.clear();
        h2o2_sink_source_.clear();
        no2_sink_source_.clear();
    }

    /**
     * @brief Accumulate oxygen sink-rate contribution for a CFD cell
     *
     * @param cell_id CFD cell (preCICE vertex id)
     * @param sink_rate Per-cell oxygen concentration sink-rate contribution [O2-units/s]
     */
    void TransferOxygenSink(int cell_id, double sink_rate) {
        if (sink_rate <= 0.0) {
            return;
        }
        oxygen_sink_source_[cell_id] += sink_rate;
    }

    /**
     * @brief Accumulate NO2- sink-rate contribution for a CFD cell
     */
    void TransferH2O2Sink(int cell_id, double sink_rate) {
        if (sink_rate <= 0.0) {
            return;
        }
        h2o2_sink_source_[cell_id] += sink_rate;
    }

    void TransferNO2Sink(int cell_id, double sink_rate) {
        if (sink_rate <= 0.0) {
            return;
        }
        no2_sink_source_[cell_id] += sink_rate;
    }

    /**
     * @brief Check if two-way coupling is enabled
     */
    bool IsTwoWayCoupling() const {
        return params_.coupling.two_way_coupling;
    }

    /**
     * @brief Get coupling configuration
     */
    const typename AdapterParameters::CouplingConfig& GetCouplingConfig() const {
        return params_.coupling;
    }

    /**
     * @brief Get fluid configuration
     */
    const typename AdapterParameters::FluidConfig& GetFluidConfig() const {
        return params_.fluid;
    }

    /**
     * @brief Get particle physics configuration
     */
    const typename AdapterParameters::ParticlePhysicsConfig& GetParticlePhysicsConfig() const {
        return params_.particle_physics;
    }

    // =========================================================================
    // Configuration Access
    // =========================================================================

    /**
     * @brief Get the coupling configuration
     */
    const Config& GetConfig() const { return config_; }

    /**
     * @brief Get mutable configuration reference
     */
    Config& GetConfig() { return config_; }

    /**
     * @brief Get adapter parameters
     */
    const AdapterParameters& GetParameters() const { return params_; }

    /**
     * @brief Get domain bounding box
     */
    const BoundingBox& GetDomainBounds() const { return domain_bounds_; }

    /**
     * @brief Set domain bounding box
     */
    void SetDomainBounds(const BoundingBox& bounds) { domain_bounds_ = bounds; }

    /**
     * @brief Override domain bounds with explicit values (for degenerate mesh workaround)
     *
     * Use this when preCICE mesh bounds are incorrect (e.g., axis mismatch,
     * partial MPI partition). Logs detailed info about the override.
     *
     * @param x_min Minimum X coordinate
     * @param x_max Maximum X coordinate
     * @param y_min Minimum Y coordinate
     * @param y_max Maximum Y coordinate
     * @param z_min Minimum Z coordinate
     * @param z_max Maximum Z coordinate
     * @param reason Description of why override is needed
     */
    void SetDomainBoundsOverride(double x_min, double x_max,
                                  double y_min, double y_max,
                                  double z_min, double z_max,
                                  const std::string& reason = "manual override") {
        BoundingBox old_bounds = domain_bounds_;
        domain_bounds_.min = {x_min, y_min, z_min};
        domain_bounds_.max = {x_max, y_max, z_max};

        bdm::Log::Warning("PreciceAdapter", "DOMAIN BOUNDS OVERRIDE (", reason, "):");
        bdm::Log::Warning("PreciceAdapter", "  Old bounds: [",
                         old_bounds.min[0], ",", old_bounds.min[1], ",", old_bounds.min[2], "] -> [",
                         old_bounds.max[0], ",", old_bounds.max[1], ",", old_bounds.max[2], "]");
        bdm::Log::Warning("PreciceAdapter", "  New bounds: [",
                         x_min, ",", y_min, ",", z_min, "] -> [",
                         x_max, ",", y_max, ",", z_max, "]");
        bdm::Log::Warning("PreciceAdapter", "  New extents: X=", x_max - x_min,
                         ", Y=", y_max - y_min, ", Z=", z_max - z_min);
    }

    /**
     * @brief Check if current domain bounds appear degenerate
     *
     * @param min_extent Minimum acceptable extent in X or Y (default 0.1)
     * @return true if any X/Y extent is below min_extent
     */
    bool AreBoundsDegenerate(double min_extent = 0.1) const {
        auto size = domain_bounds_.Size();
        return size[0] < min_extent || size[1] < min_extent;
    }

    /**
     * @brief Get mesh dimensions from preCICE for diagnostics
     * @param mesh_name Name of the mesh to query
     * @return Number of dimensions (2 or 3), or -1 if mesh not found
     */
    int GetMeshDimensions(const std::string& mesh_name) const {
        if (!participant_) return -1;
        try {
            return participant_->getMeshDimensions(mesh_name);
        } catch (...) {
            return -1;
        }
    }

    const std::shared_ptr<const BoundaryGeometry>& GetBoundaryGeometry() const {
        return boundary_geometry_;
    }

    bool HasExactBoundaryGeometry() const {
        return boundary_geometry_ && !boundary_geometry_->Empty();
    }

    bool HasClosedBoundaryGeometry() const {
        return boundary_geometry_ && boundary_geometry_->IsClosed();
    }

    /**
     * @brief Get wall meshes for collision detection
     */
    const std::map<std::string, WallMeshData>& GetWallMeshes() const {
        return wall_meshes_;
    }

    /**
     * @brief Check if wall mesh data is available
     */
    bool HasWallMeshData() const {
        for (const auto& [name, mesh] : wall_meshes_) {
            if (mesh.is_initialized && mesh.NumVertices() > 0) {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Structure representing a single wall face for collision detection
     *
     * Contains all data needed for accurate particle-wall interactions:
     * - center: Face center position for distance calculations
     * - normal: OpenFOAM patch.nf() unit normal, outward from fluid domain
     * - area: Face area for proximity filtering
     * - velocity: Face velocity for moving walls (e.g., lid-driven cavity)
     */
    struct WallFace {
        std::array<double, 3> center;   // Face center position
        std::array<double, 3> normal;   // OpenFOAM patch.nf() unit normal, outward from fluid domain
        double area;                     // Face area
        std::array<double, 3> velocity; // Face velocity (for moving walls)

        // Check if this is a moving wall
        bool IsMoving() const {
            return (std::abs(velocity[0]) > 1e-12 ||
                    std::abs(velocity[1]) > 1e-12 ||
                    std::abs(velocity[2]) > 1e-12);
        }
    };

    /**
     * @brief Get all wall faces for collision detection
     *
     * Returns a vector of wall faces with their centers, normals, areas, and velocities.
     * This data comes from OpenFOAM via preCICE WallSurfaceMesh and can be
     * used for accurate particle-wall collision detection in any geometry
     * (not just axis-aligned boxes).
     *
     * For moving walls (e.g., lid-driven cavity), the velocity field contains
     * the boundary condition velocity, enabling correct relative velocity
     * calculations in collision responses.
     *
     * @return Vector of WallFace structs
     */
    std::vector<WallFace> GetWallFaces() const {
        std::vector<WallFace> faces;

        for (const auto& [mesh_name, mesh] : wall_meshes_) {
            if (!mesh.is_initialized) continue;

            // Check that we have geometry, normals, and areas
            size_t num_faces = mesh.areas.size();
            if (num_faces == 0) continue;
            if (mesh.geometry.size() < num_faces * 3) continue;
            if (mesh.normals.size() < num_faces * 3) continue;

            // Check if velocity data is available
            bool has_velocity = (mesh.velocity.size() >= num_faces * 3);

            faces.reserve(faces.size() + num_faces);

            size_t moving_count = 0;
            for (size_t i = 0; i < num_faces; ++i) {
                WallFace face;
                face.center = {
                    mesh.geometry[i * 3],
                    mesh.geometry[i * 3 + 1],
                    mesh.geometry[i * 3 + 2]
                };
                face.normal = {
                    mesh.normals[i * 3],
                    mesh.normals[i * 3 + 1],
                    mesh.normals[i * 3 + 2]
                };
                face.area = mesh.areas[i];

                // Get velocity if available, otherwise zero (stationary wall)
                if (has_velocity) {
                    face.velocity = {
                        mesh.velocity[i * 3],
                        mesh.velocity[i * 3 + 1],
                        mesh.velocity[i * 3 + 2]
                    };
                    if (face.IsMoving()) ++moving_count;
                } else {
                    face.velocity = {0.0, 0.0, 0.0};
                }

                faces.push_back(face);
            }

            if (moving_count > 0) {
                bdm::Log::Info("PreciceAdapter", "Loaded ", num_faces,
                              " wall faces from '", mesh_name, "' (",
                              moving_count, " moving)");
            } else {
                bdm::Log::Info("PreciceAdapter", "Loaded ", num_faces,
                              " wall faces from '", mesh_name, "' (all stationary)");
            }
        }

        return faces;
    }

    /**
     * @brief Get total number of wall faces across all wall meshes
     */
    size_t GetTotalWallFaceCount() const {
        size_t count = 0;
        for (const auto& [name, mesh] : wall_meshes_) {
            if (mesh.is_initialized) {
                count += mesh.areas.size();
            }
        }
        return count;
    }

    /**
     * @brief Get a specific volume mesh
     */
    const MeshData* GetVolumeMesh(const std::string& name) const {
        auto it = volume_meshes_.find(name);
        return (it != volume_meshes_.end()) ? &it->second : nullptr;
    }

    /**
     * @brief Clamp all dynamic agents to domain bounds
     *
     * Ensures agents stay within the fluid domain. Useful for
     * preventing particles from escaping the CFD domain.
     */
    void ClampAgentsToDomain() {
        // Require finite, non-degenerate bounds on every axis before clamping.
        // The previous guard (`!min[0] && !max[0]`) silently disabled clamping
        // whenever the X bounds happened to contain zero — including the
        // degenerate [0,0] produced by a broken wall-mesh exchange, which is
        // precisely when clamping is most needed.
        for (int axis = 0; axis < 3; ++axis) {
            if (!std::isfinite(domain_bounds_.min[axis]) ||
                !std::isfinite(domain_bounds_.max[axis]) ||
                !(domain_bounds_.max[axis] > domain_bounds_.min[axis])) {
                return; // bounds not set or degenerate on this axis
            }
        }

        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;

        rm->ForEachAgent([&](bdm::Agent* agent) {
            auto* cell = dynamic_cast<CellType*>(agent);
            if (!cell || cell->IsStationary()) return;

            auto position = cell->GetPosition();
            bool modified = false;
            double diameter = std::max(0.0, params_.agent.default_diameter);
            if constexpr (detail::has_get_diameter<CellType>::value) {
                const double cell_diameter = static_cast<double>(cell->GetDiameter());
                if (std::isfinite(cell_diameter) && cell_diameter > 0.0) {
                    diameter = cell_diameter;
                }
            }
            const double radius = 0.5 * diameter;

            for (int axis = 0; axis < 3; ++axis) {
                // Keep the full agent sphere inside the domain when possible:
                // center must satisfy [min + r, max - r].
                double lower = domain_bounds_.min[axis] + radius;
                double upper = domain_bounds_.max[axis] - radius;
                if (!(upper > lower)) {
                    // If domain thickness is smaller than diameter on this axis,
                    // fall back to center-only clamping to avoid invalid bounds.
                    lower = domain_bounds_.min[axis];
                    upper = domain_bounds_.max[axis];
                }
                double clamped = std::clamp(position[axis], lower, upper);
                if (std::abs(clamped - position[axis]) > 1e-12) {
                    position[axis] = clamped;
                    modified = true;
                }
            }

            // Force 2D if configured
            if (params_.coupling.force_2d) {
                double z_plane = (domain_bounds_.min[2] + domain_bounds_.max[2]) / 2.0;
                if (std::abs(position[2] - z_plane) > 1e-12) {
                    position[2] = z_plane;
                    modified = true;
                }
            }

            if (modified) {
                cell->SetPosition(position);
            }
        });
    }

    // =========================================================================
    // Compatibility Aliases (for legacy code migration)
    // =========================================================================

    /**
     * @brief Alias for Initialize() - for backward compatibility
     */
    void InitializeCoupling() { Initialize(); }

    /**
     * @brief Alias for ClampAgentsToDomain() - for backward compatibility
     */
    void ClampDynamicAgentsToDomain() { ClampAgentsToDomain(); }

    /**
     * @brief Update dynamic agent-mesh mapping
     *
     * IMPORTANT: Dynamic mapping is handled automatically by preCICE's Just-In-Time (JIT) mapping!
     *
     * How dynamic mapping works:
     * 1. ONE-WAY COUPLING (reading fluid data):
     *    - Each timestep, FetchInputData() calls mapAndReadData() with CURRENT agent positions
     *    - preCICE dynamically interpolates fluid fields at these moving positions
     *    - No pre-registration needed - positions update every call
     *
     * 2. TWO-WAY COUPLING (writing source terms):
     *    - CollectSourceTerms() calls FindNearestMeshVertex() with CURRENT agent positions
     *    - Each particle finds its containing CFD cell dynamically
     *    - Source terms are accumulated to the correct cells based on current positions
     *
     * This method provides optional validation and diagnostics:
     * - Counts agents outside domain bounds (warning)
     * - Tracks mapping statistics for debugging
     * - Can be extended for custom remapping strategies
     */
    void UpdateDynamicMapping() {
        if (!is_initialized_) return;

        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;

        // Validate agent positions are within reasonable bounds
        size_t out_of_bounds = 0;
        size_t dynamic_agents = 0;

        rm->ForEachAgent([&](bdm::Agent* agent) {
            auto* cell = dynamic_cast<CellType*>(agent);
            if (!cell || cell->IsStationary()) return;

            dynamic_agents++;

            // Check if position is within domain bounds
            const auto& pos = cell->GetPosition();
            if (!domain_bounds_.Contains({pos[0], pos[1], pos[2]})) {
                out_of_bounds++;
                if (out_of_bounds <= 5) {
                    bdm::Log::Warning("PreciceAdapter",
                        "Out-of-bounds agent pos=(",
                        pos[0], ",", pos[1], ",", pos[2], "), bounds X:[",
                        domain_bounds_.min[0], ",", domain_bounds_.max[0], "] Y:[",
                        domain_bounds_.min[1], ",", domain_bounds_.max[1], "] Z:[",
                        domain_bounds_.min[2], ",", domain_bounds_.max[2], "]");
                }
            }
        });

        if (out_of_bounds > 0 && dynamic_agents > 0) {
            double fraction = static_cast<double>(out_of_bounds) / dynamic_agents;
            if (fraction > 0.1) {  // More than 10% out of bounds
                bdm::Log::Warning("PreciceAdapter",
                    "Dynamic mapping: ", out_of_bounds, " / ", dynamic_agents,
                    " agents outside domain bounds (", int(fraction * 100), "%)" );
            }
        }

        bdm::Log::Info("PreciceAdapter",
            "Dynamic mapping OK: ", dynamic_agents, " agents tracked");

        // Invalidate nearest-cell cache: agents may have moved, so the cached
        // FVM-cell assignment from the previous window is stale.
        // For stationary simulations (o2Cube) UpdateDynamicMapping is never
        // called (motion_enabled=false in cells.h), so the cache persists for
        // the full run — which is exactly the desired behaviour.
        agent_fvm_cell_cache_.clear();
    }

    // =========================================================================
    // Agent Management
    // =========================================================================

    /**
     * @brief Create agents at mesh vertex positions
     *
     * Creates CellType agents at each vertex of the specified volume mesh.
     * Useful for initializing particles at CFD cell centers.
     *
     * @param mesh_name Name of the volume mesh
     * @param agent_factory Factory function to create configured agents
     * @return Number of agents created
     */
    size_t CreateAgentsAtMeshVertices(
        const std::string& mesh_name,
        std::function<CellType*(const std::array<double, 3>&)> agent_factory)
    {
        auto it = volume_meshes_.find(mesh_name);
        if (it == volume_meshes_.end() || !it->second.is_initialized) {
            bdm::Log::Warning("PreciceAdapter", "Mesh '", mesh_name, "' not found or not initialized");
            return 0;
        }

        const auto& mesh = it->second;
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();

        size_t count = 0;
        for (size_t i = 0; i < mesh.NumVertices(); ++i) {
            auto pos = mesh.GetVertex(i);
            auto* agent = agent_factory(pos);
            if (agent) {
                rm->AddAgent(agent);
                ++count;
            }
        }

        bdm::Log::Info("PreciceAdapter", "Created ", count, " agents from mesh '", mesh_name, "'");
        return count;
    }

    /**
     * @brief Create agents from CSV file positions
     *
     * @param filepath Path to CSV file with x,y,z columns
     * @param agent_factory Factory function to create configured agents
     * @return Number of agents created
     */
    /**
     * @brief Load agent positions from CSV file
     * Expected format: x,y,z (with header row)
     */
    std::vector<std::array<double, 3>> LoadPositionsFromCSV(const std::string& filepath) {
        std::vector<std::array<double, 3>> positions;

        bdm::Log::Info("PreciceAdapter", "Attempting to load positions from: ", filepath);
        std::ifstream file(filepath);

        if (!file.is_open()) {
            bdm::Log::Warning("PreciceAdapter", "Could not open CSV file: ", filepath);
            return positions;
        }

        bdm::Log::Info("PreciceAdapter", "CSV file opened successfully");
        std::string line;
        bool first_line = true;
        int line_num = 0;

        while (std::getline(file, line)) {
            line_num++;

            // Skip empty lines and comments
            if (line.empty() || line[0] == '#') {
                bdm::Log::Info("PreciceAdapter", "Line ", line_num, " skipped (comment/empty): ", line);
                continue;
            }

            // Skip header line
            if (first_line) {
                bdm::Log::Info("PreciceAdapter", "Line ", line_num, " skipped (header): ", line);
                first_line = false;
                continue;
            }

            bdm::Log::Info("PreciceAdapter", "Parsing line ", line_num, ": ", line);
            std::istringstream iss(line);
            std::string token;
            std::array<double, 3> pos;
            int idx = 0;

            while (std::getline(iss, token, ',') && idx < 3) {
                try {
                    pos[idx] = std::stod(token);
                    bdm::Log::Info("PreciceAdapter", "  Parsed pos[", idx, "] = ", pos[idx]);
                    idx++;
                } catch (...) {
                    bdm::Log::Warning("PreciceAdapter", "Failed to parse value: ", token);
                    break;
                }
            }

            if (idx == 3) {
                positions.push_back(pos);
                bdm::Log::Info("PreciceAdapter", "Added position: (", pos[0], ", ", pos[1], ", ", pos[2], ")");
            } else {
                bdm::Log::Warning("PreciceAdapter", "Line ", line_num, " incomplete (", idx, " values)");
            }
        }

        bdm::Log::Info("PreciceAdapter", "Loaded ", positions.size(), " positions from CSV");
        return positions;
    }

    size_t CreateAgentsFromCSV(
        const std::string& filepath,
        std::function<CellType*(const std::array<double, 3>&)> agent_factory)
    {
        auto positions = LoadPositionsFromCSV(filepath);
        if (positions.empty()) {
            bdm::Log::Warning("PreciceAdapter", "No positions loaded from '", filepath, "'");
            return 0;
        }

        // Optionally remove positions that are adjacent to walls.
        // Effective distance can be specified directly, or derived from
        // (wall_clearance_layers * wall_clearance_layer_thickness).
        double wall_clearance_distance = params_.agent.wall_clearance_distance;
        if (params_.agent.wall_clearance_layers > 0 &&
            params_.agent.wall_clearance_layer_thickness > 0.0) {
            wall_clearance_distance = std::max(
                wall_clearance_distance,
                static_cast<double>(params_.agent.wall_clearance_layers) *
                    params_.agent.wall_clearance_layer_thickness);
        }
        if (params_.agent.enable_wall_clearance && wall_clearance_distance > 0.0) {
            const size_t before = positions.size();
            positions = FilterWallAdjacentPositions(
                positions, wall_meshes_, wall_clearance_distance);
            const size_t removed = before - positions.size();
            bdm::Log::Info("PreciceAdapter",
                "Wall-clearance filter (d=", wall_clearance_distance, " m): ",
                "kept ", positions.size(), " / ", before,
                " positions (removed ", removed, " wall-adjacent)");
        }

        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        size_t count = 0;

        for (const auto& pos : positions) {
            auto* agent = agent_factory(pos);
            if (agent) {
                rm->AddAgent(agent);
                ++count;
            }
        }

        bdm::Log::Info("PreciceAdapter", "Created ", count, " agents from '", filepath, "'");
        return count;
    }

    // =========================================================================
    // Inlet/Outlet Agent Management (for complex geometries)
    // =========================================================================

    std::vector<size_t> SelectInletFaceIndices(const MeshData& mesh) const {
        std::vector<size_t> all_indices(mesh.NumVertices());
        std::iota(all_indices.begin(), all_indices.end(), size_t{0});

        if (!params_.agent.inlet_keep_largest_component || mesh.NumVertices() < 2) {
            return all_indices;
        }

        const size_t n = mesh.NumVertices();
        std::vector<double> nearest_distances;
        nearest_distances.reserve(n);

        for (size_t i = 0; i < n; ++i) {
            const auto pi = mesh.GetVertex(i);
            double nearest_d2 = std::numeric_limits<double>::max();
            for (size_t j = 0; j < n; ++j) {
                if (i == j) continue;
                const auto pj = mesh.GetVertex(j);
                const double dx = pi[0] - pj[0];
                const double dy = pi[1] - pj[1];
                const double dz = pi[2] - pj[2];
                const double d2 = dx * dx + dy * dy + dz * dz;
                nearest_d2 = std::min(nearest_d2, d2);
            }
            if (nearest_d2 < std::numeric_limits<double>::max()) {
                nearest_distances.push_back(std::sqrt(nearest_d2));
            }
        }

        if (nearest_distances.empty()) {
            return all_indices;
        }

        std::sort(nearest_distances.begin(), nearest_distances.end());
        const double median_nn = nearest_distances[nearest_distances.size() / 2];
        double radius = params_.agent.inlet_component_radius;
        if (radius <= 0.0) {
            radius = std::max(3.0 * median_nn, 1.0e-12);
        }
        const double radius2 = radius * radius;

        std::vector<int> component(n, -1);
        std::vector<size_t> component_sizes;
        int current_component = 0;

        for (size_t seed = 0; seed < n; ++seed) {
            if (component[seed] != -1) continue;

            component_sizes.push_back(0);
            std::vector<size_t> stack{seed};
            component[seed] = current_component;

            while (!stack.empty()) {
                const size_t i = stack.back();
                stack.pop_back();
                ++component_sizes.back();
                const auto pi = mesh.GetVertex(i);

                for (size_t j = 0; j < n; ++j) {
                    if (component[j] != -1) continue;
                    const auto pj = mesh.GetVertex(j);
                    const double dx = pi[0] - pj[0];
                    const double dy = pi[1] - pj[1];
                    const double dz = pi[2] - pj[2];
                    const double d2 = dx * dx + dy * dy + dz * dz;
                    if (d2 <= radius2) {
                        component[j] = current_component;
                        stack.push_back(j);
                    }
                }
            }

            ++current_component;
        }

        if (component_sizes.size() <= 1) {
            return all_indices;
        }

        const auto largest_it = std::max_element(component_sizes.begin(), component_sizes.end());
        const int largest_component = static_cast<int>(std::distance(component_sizes.begin(), largest_it));

        std::vector<size_t> selected;
        selected.reserve(*largest_it);
        for (size_t i = 0; i < n; ++i) {
            if (component[i] == largest_component) {
                selected.push_back(i);
            }
        }

        bdm::Log::Warning("PreciceAdapter",
            "Inlet component filter kept largest connected face cluster: ",
            selected.size(), "/", n, " faces (radius=", radius,
            " m, components=", component_sizes.size(), "). Removed ",
            n - selected.size(), " disconnected inlet faces.");

        return selected;
    }

    std::vector<std::array<double, 3>> KeepLargestPositionComponent(
        const std::vector<std::array<double, 3>>& positions,
        const std::string& label) const {
        if (!params_.agent.inlet_keep_largest_component || positions.size() < 2) {
            return positions;
        }

        const size_t n = positions.size();
        std::vector<double> nearest_distances;
        nearest_distances.reserve(n);

        for (size_t i = 0; i < n; ++i) {
            double nearest_d2 = std::numeric_limits<double>::max();
            for (size_t j = 0; j < n; ++j) {
                if (i == j) continue;
                const double dx = positions[i][0] - positions[j][0];
                const double dy = positions[i][1] - positions[j][1];
                const double dz = positions[i][2] - positions[j][2];
                nearest_d2 = std::min(nearest_d2, dx * dx + dy * dy + dz * dz);
            }
            if (nearest_d2 < std::numeric_limits<double>::max()) {
                nearest_distances.push_back(std::sqrt(nearest_d2));
            }
        }

        if (nearest_distances.empty()) {
            return positions;
        }

        std::sort(nearest_distances.begin(), nearest_distances.end());
        const double median_nn = nearest_distances[nearest_distances.size() / 2];
        double radius = params_.agent.inlet_component_radius;
        if (radius <= 0.0) {
            radius = std::max(3.0 * median_nn, 1.0e-12);
        }
        const double radius2 = radius * radius;

        std::vector<int> component(n, -1);
        std::vector<size_t> component_sizes;
        int current_component = 0;

        for (size_t seed = 0; seed < n; ++seed) {
            if (component[seed] != -1) continue;

            component_sizes.push_back(0);
            std::vector<size_t> stack{seed};
            component[seed] = current_component;

            while (!stack.empty()) {
                const size_t i = stack.back();
                stack.pop_back();
                ++component_sizes.back();

                for (size_t j = 0; j < n; ++j) {
                    if (component[j] != -1) continue;
                    const double dx = positions[i][0] - positions[j][0];
                    const double dy = positions[i][1] - positions[j][1];
                    const double dz = positions[i][2] - positions[j][2];
                    const double d2 = dx * dx + dy * dy + dz * dz;
                    if (d2 <= radius2) {
                        component[j] = current_component;
                        stack.push_back(j);
                    }
                }
            }

            ++current_component;
        }

        if (component_sizes.size() <= 1) {
            return positions;
        }

        const auto largest_it = std::max_element(component_sizes.begin(), component_sizes.end());
        const int largest_component = static_cast<int>(std::distance(component_sizes.begin(), largest_it));

        std::vector<std::array<double, 3>> selected;
        selected.reserve(*largest_it);
        for (size_t i = 0; i < n; ++i) {
            if (component[i] == largest_component) {
                selected.push_back(positions[i]);
            }
        }

        bdm::Log::Warning("PreciceAdapter",
            "Inlet component filter kept largest connected ", label,
            " cluster: ", selected.size(), "/", n,
            " positions (radius=", radius, " m, components=",
            component_sizes.size(), "). Removed ", n - selected.size(),
            " disconnected mapped positions.");

        return selected;
    }

    /**
     * @brief Create agents at inlet volume cell centers
     *
     * For complex geometries (aneurysm, venturi, etc.) the OpenFOAM adapter
     * can export the owner volume-cell center associated with each inlet
     * boundary face. This method uses those exact centers when available and
     * only falls back to a guarded nearest-neighbor search for older configs.
     *
     * Wall-adjacent volume cells (cells that also touch the wall boundary)
     * are omitted via the wall-clearance filter.
     *
     * @param mesh_name Name of the inlet surface mesh (default: params_.mesh.inlet_mesh)
     * @param agent_factory Factory function to create configured agents
     * @param offset_into_domain Legacy parameter (kept for API compat; ignored when
     *        volume mesh mapping is available)
     * @return Number of agents created
     */
    size_t CreateAgentsAtInletMesh(
        const std::string& mesh_name,
        std::function<CellType*(const std::array<double, 3>&)> agent_factory,
        double offset_into_domain = 0.0,
        size_t max_agents = 0)
    {
        // Try inlet_meshes_ first, then fall back to volume_meshes_ for compatibility
        const MeshData* mesh_ptr = nullptr;

        auto inlet_it = inlet_meshes_.find(mesh_name);
        if (inlet_it != inlet_meshes_.end() && inlet_it->second.is_initialized) {
            mesh_ptr = &inlet_it->second;
        } else {
            auto vol_it = volume_meshes_.find(mesh_name);
            if (vol_it != volume_meshes_.end() && vol_it->second.is_initialized) {
                mesh_ptr = &vol_it->second;
            }
        }

        if (!mesh_ptr) {
            bdm::Log::Warning("PreciceAdapter", "Inlet mesh '", mesh_name, "' not found or not initialized");
            return 0;
        }

        const auto& mesh = *mesh_ptr;
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();

        std::vector<std::array<double, 3>> candidates;
        const std::vector<size_t> inlet_face_indices = SelectInletFaceIndices(mesh);

        const auto& ib = mesh.bounds;
        const double inlet_dx = ib.max[0] - ib.min[0];
        const double inlet_dy = ib.max[1] - ib.min[1];
        const double inlet_dz = ib.max[2] - ib.min[2];
        const double inlet_diag =
            std::sqrt(inlet_dx * inlet_dx + inlet_dy * inlet_dy + inlet_dz * inlet_dz);
        double max_map_distance = 2.0 * inlet_diag;
        if (params_.agent.inlet_max_cell_distance > 0.0) {
            max_map_distance = params_.agent.inlet_max_cell_distance;
        }
        const double max_map_d2 = max_map_distance * max_map_distance;

        if (mesh.HasMappedCellCenters()) {
            candidates.reserve(inlet_face_indices.size());
            size_t rejected_count = 0;
            double max_accepted_d = 0.0;
            double max_rejected_d = 0.0;

            for (const size_t i : inlet_face_indices) {
                const auto face_center = mesh.GetVertex(i);
                const auto pos = mesh.GetMappedCellCenter(i);
                const double dx = face_center[0] - pos[0];
                const double dy = face_center[1] - pos[1];
                const double dz = face_center[2] - pos[2];
                const double d2 = dx * dx + dy * dy + dz * dz;

                if (d2 > max_map_d2) {
                    ++rejected_count;
                    max_rejected_d = std::max(max_rejected_d, std::sqrt(d2));
                    continue;
                }

                max_accepted_d = std::max(max_accepted_d, std::sqrt(d2));
                const auto duplicate = std::find(candidates.begin(), candidates.end(), pos);
                if (duplicate == candidates.end()) {
                    candidates.push_back(pos);
                }
            }

            if (rejected_count > 0) {
                bdm::Log::Warning("PreciceAdapter",
                    "Exact inlet owner-cell mapping rejected ", rejected_count, " of ",
                    inlet_face_indices.size(), " entries because the owner cell center was farther than ",
                    max_map_distance, " m from its inlet face center (max accepted = ",
                    max_accepted_d, " m, worst rejected = ", max_rejected_d,
                    " m). Check OpenFOAM/preCICE vertex-data ordering if this warning persists.");
            }

            bdm::Log::Info("PreciceAdapter",
                "Using exact OpenFOAM owner cell centers for inlet seeding: ",
                inlet_face_indices.size(), " inlet faces -> ", candidates.size(),
                " unique adjacent volume cells");
        } else {
            // ---------------------------------------------------------------
            // Fallback: map inlet face centers to nearest FVMmesh cell centers
            // ---------------------------------------------------------------
            // This path is kept for older configurations that do not exchange
            // InletCellCenters. It is intentionally guarded by a distance
            // threshold because spatial nearest-neighbor is ambiguous in
            // folded vascular geometries.
            const MeshData* vol_mesh_ptr = nullptr;
            {
                auto vm_it = volume_meshes_.find(params_.mesh.primary_volume_mesh);
                if (vm_it != volume_meshes_.end() && vm_it->second.is_initialized &&
                    vm_it->second.NumVertices() > 0) {
                    vol_mesh_ptr = &vm_it->second;
                }
            }

            if (!vol_mesh_ptr) {
                bdm::Log::Warning("PreciceAdapter",
                    "Volume mesh '", params_.mesh.primary_volume_mesh,
                    "' not available for inlet cell-center mapping; "
                    "cannot create agents at inlet");
                return 0;
            }

            const auto& vol_mesh = *vol_mesh_ptr;

            bdm::Log::Info("PreciceAdapter",
                "Exact inlet owner cell centers are unavailable; mapping ",
                inlet_face_indices.size(), " inlet face centers to nearest volume cell centers from '",
                vol_mesh.name, "' (", vol_mesh.NumVertices(), " cells)");

            const size_t inlet_faces = inlet_face_indices.size();
            const size_t volume_cells = vol_mesh.NumVertices();
            std::vector<size_t> best_indices(inlet_faces, std::numeric_limits<size_t>::max());
            std::vector<double> best_dist2(inlet_faces, std::numeric_limits<double>::max());

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for (long long fi_signed = 0; fi_signed < static_cast<long long>(inlet_faces); ++fi_signed) {
                const size_t fi = static_cast<size_t>(fi_signed);
                const auto fc = mesh.GetVertex(inlet_face_indices[fi]);

                size_t local_best_idx = 0;
                double local_best_d2 = std::numeric_limits<double>::max();

                for (size_t vi = 0; vi < volume_cells; ++vi) {
                    const auto vc = vol_mesh.GetVertex(vi);
                    const double dx = fc[0] - vc[0];
                    const double dy = fc[1] - vc[1];
                    const double dz = fc[2] - vc[2];
                    const double d2 = dx * dx + dy * dy + dz * dz;
                    if (d2 < local_best_d2) {
                        local_best_d2 = d2;
                        local_best_idx = vi;
                    }
                }

                best_indices[fi] = local_best_idx;
                best_dist2[fi] = local_best_d2;
            }

            std::vector<bool> selected(volume_cells, false);
            size_t mapped_count = 0;
            size_t rejected_count = 0;
            double max_accepted_d = 0.0;
            double max_rejected_d = 0.0;

            for (size_t fi = 0; fi < inlet_faces; ++fi) {
                const size_t best_idx = best_indices[fi];
                const double best_d2 = best_dist2[fi];

                if (best_d2 > max_map_d2) {
                    ++rejected_count;
                    max_rejected_d = std::max(max_rejected_d, std::sqrt(best_d2));
                    continue;
                }

                max_accepted_d = std::max(max_accepted_d, std::sqrt(best_d2));
                if (!selected[best_idx]) {
                    selected[best_idx] = true;
                    ++mapped_count;
                }
            }

            if (rejected_count > 0) {
                bdm::Log::Warning("PreciceAdapter",
                    "Inlet-to-cell fallback mapping rejected ", rejected_count, " of ",
                    inlet_faces, " inlet face centers because their nearest ",
                    "FVMmesh cell was farther than ", max_map_distance,
                    " m (max accepted = ", max_accepted_d,
                    " m, worst rejected = ", max_rejected_d,
                    " m, inlet bbox diag = ", inlet_diag, " m). ",
                    "Use InletCellCenters exchange for exact OpenFOAM owner-cell seeding.");
            }

            candidates.reserve(mapped_count);
            for (size_t vi = 0; vi < vol_mesh.NumVertices(); ++vi) {
                if (selected[vi]) {
                    candidates.push_back(vol_mesh.GetVertex(vi));
                }
            }

            bdm::Log::Info("PreciceAdapter",
                "Inlet fallback cell-center mapping: ", inlet_faces,
                " face centers -> ", candidates.size(),
                " unique volume cells");
        }

        candidates = KeepLargestPositionComponent(candidates, "mapped inlet cell-center");

        if (!candidates.empty() && !inlet_face_indices.empty()) {
            std::vector<std::array<double, 3>> inlet_adjacent_candidates;
            inlet_adjacent_candidates.reserve(candidates.size());
            size_t rejected_count = 0;
            double max_accepted_d = 0.0;
            double max_rejected_d = 0.0;

            for (const auto& pos : candidates) {
                double nearest_d2 = std::numeric_limits<double>::max();
                for (const size_t face_index : inlet_face_indices) {
                    const auto face_center = mesh.GetVertex(face_index);
                    const double dx = pos[0] - face_center[0];
                    const double dy = pos[1] - face_center[1];
                    const double dz = pos[2] - face_center[2];
                    nearest_d2 = std::min(nearest_d2, dx * dx + dy * dy + dz * dz);
                }

                if (nearest_d2 <= max_map_d2) {
                    inlet_adjacent_candidates.push_back(pos);
                    max_accepted_d = std::max(max_accepted_d, std::sqrt(nearest_d2));
                } else {
                    ++rejected_count;
                    max_rejected_d = std::max(max_rejected_d, std::sqrt(nearest_d2));
                }
            }

            if (rejected_count > 0) {
                bdm::Log::Warning("PreciceAdapter",
                    "Final inlet adjacency filter rejected ", rejected_count, " of ",
                    candidates.size(), " candidate cell centers because they were farther than ",
                    max_map_distance, " m from the selected inlet face cloud (max accepted = ",
                    max_accepted_d, " m, worst rejected = ", max_rejected_d, " m).");
            }
            candidates = std::move(inlet_adjacent_candidates);
        }

        // Keep only owner volume-cell centers near the physical inlet center.
        // Surface meshes are generally nonuniform, so an ordinary mean or
        // median of face/owner centers can visibly shift the seed disk.  When
        // OpenFOAM face areas are available, calculate the center and robust
        // radius in the owner-cell layer with face-area weights.  This also
        // keeps the center correctly offset inside an arbitrarily oriented
        // inlet instead of placing it on the boundary surface.
        const double center_radius_fraction =
            params_.agent.inlet_center_radius_fraction;
        if (!candidates.empty() && center_radius_fraction > 0.0 &&
            center_radius_fraction < 1.0) {
            std::array<double, 3> inlet_center = {{0.0, 0.0, 0.0}};
            double effective_inlet_radius = 0.0;
            bool used_area_weighting = false;

            if (mesh.HasMappedCellCenters() && mesh.HasFaceAreas()) {
                double total_area = 0.0;
                for (const size_t face_index : inlet_face_indices) {
                    const double area = mesh.face_areas[face_index];
                    if (!std::isfinite(area) || area <= 0.0) continue;

                    const auto face_center = mesh.GetVertex(face_index);
                    const auto owner_center = mesh.GetMappedCellCenter(face_index);
                    const double mdx = face_center[0] - owner_center[0];
                    const double mdy = face_center[1] - owner_center[1];
                    const double mdz = face_center[2] - owner_center[2];
                    if (mdx * mdx + mdy * mdy + mdz * mdz > max_map_d2) continue;

                    inlet_center[0] += area * owner_center[0];
                    inlet_center[1] += area * owner_center[1];
                    inlet_center[2] += area * owner_center[2];
                    total_area += area;
                }

                if (total_area > 0.0) {
                    inlet_center[0] /= total_area;
                    inlet_center[1] /= total_area;
                    inlet_center[2] /= total_area;

                    std::vector<std::pair<double, double>> distance_areas;
                    distance_areas.reserve(inlet_face_indices.size());
                    double radius_area = 0.0;
                    for (const size_t face_index : inlet_face_indices) {
                        const double area = mesh.face_areas[face_index];
                        if (!std::isfinite(area) || area <= 0.0) continue;

                        const auto face_center = mesh.GetVertex(face_index);
                        const auto owner_center = mesh.GetMappedCellCenter(face_index);
                        const double mdx = face_center[0] - owner_center[0];
                        const double mdy = face_center[1] - owner_center[1];
                        const double mdz = face_center[2] - owner_center[2];
                        if (mdx * mdx + mdy * mdy + mdz * mdz > max_map_d2) continue;

                        const double dx = owner_center[0] - inlet_center[0];
                        const double dy = owner_center[1] - inlet_center[1];
                        const double dz = owner_center[2] - inlet_center[2];
                        distance_areas.emplace_back(
                            std::sqrt(dx * dx + dy * dy + dz * dz), area);
                        radius_area += area;
                    }

                    std::sort(distance_areas.begin(), distance_areas.end(),
                              [](const auto& lhs, const auto& rhs) {
                                  return lhs.first < rhs.first;
                              });
                    const double target_area = 0.95 * radius_area;
                    double cumulative_area = 0.0;
                    for (const auto& distance_area : distance_areas) {
                        cumulative_area += distance_area.second;
                        effective_inlet_radius = distance_area.first;
                        if (cumulative_area >= target_area) break;
                    }
                    used_area_weighting = !distance_areas.empty() &&
                                          effective_inlet_radius > 0.0;
                }
            }

            if (!used_area_weighting) {
                const auto median_coordinate = [&](const int axis) {
                    std::vector<double> values;
                    values.reserve(candidates.size());
                    for (const auto& pos : candidates) values.push_back(pos[axis]);
                    const size_t middle = values.size() / 2;
                    std::nth_element(values.begin(), values.begin() + middle,
                                     values.end());
                    if (values.size() % 2 == 1) return values[middle];
                    const double upper = values[middle];
                    std::nth_element(values.begin(), values.begin() + middle - 1,
                                     values.begin() + middle);
                    return 0.5 * (values[middle - 1] + upper);
                };
                inlet_center = {{median_coordinate(0), median_coordinate(1),
                                 median_coordinate(2)}};

                std::vector<double> radial_distances;
                radial_distances.reserve(candidates.size());
                for (const auto& pos : candidates) {
                    const double dx = pos[0] - inlet_center[0];
                    const double dy = pos[1] - inlet_center[1];
                    const double dz = pos[2] - inlet_center[2];
                    radial_distances.push_back(
                        std::sqrt(dx * dx + dy * dy + dz * dz));
                }
                const size_t radius_index = static_cast<size_t>(std::floor(
                    0.95 * static_cast<double>(radial_distances.size() - 1)));
                std::nth_element(radial_distances.begin(),
                                 radial_distances.begin() + radius_index,
                                 radial_distances.end());
                effective_inlet_radius = radial_distances[radius_index];
                bdm::Log::Warning(
                    "PreciceAdapter",
                    "Inlet face areas unavailable; using median-based inlet center fallback");
            }

            const double selected_radius =
                center_radius_fraction * effective_inlet_radius;
            const double selected_radius2 = selected_radius * selected_radius;

            std::vector<std::array<double, 3>> central_candidates;
            central_candidates.reserve(candidates.size());
            for (const auto& pos : candidates) {
                const double dx = pos[0] - inlet_center[0];
                const double dy = pos[1] - inlet_center[1];
                const double dz = pos[2] - inlet_center[2];
                if (dx * dx + dy * dy + dz * dz <= selected_radius2) {
                    central_candidates.push_back(pos);
                }
            }

            bdm::Log::Info(
                "PreciceAdapter", "Central inlet volume-cell filter: center=(",
                inlet_center[0], ", ", inlet_center[1], ", ",
                inlet_center[2], "), effective radius=",
                effective_inlet_radius, " m, selected radius=",
                selected_radius, " m (fraction=", center_radius_fraction,
                ", center method=",
                used_area_weighting ? "face-area-weighted" : "median fallback",
                "): kept ", central_candidates.size(), " / ",
                candidates.size(), " owner cell centers");
            if (!central_candidates.empty()) {
                candidates = std::move(central_candidates);
            } else {
                bdm::Log::Warning(
                    "PreciceAdapter",
                    "Central inlet filter selected no owner cells; retaining the unfiltered inlet set");
            }
        } else if (center_radius_fraction <= 0.0) {
            bdm::Log::Warning(
                "PreciceAdapter",
                "Ignoring non-positive agent.inlet_center_radius_fraction=",
                center_radius_fraction, " (expected 0 < fraction <= 1)");
        }

        // Filter out volume cells that are adjacent to walls.
        // These are cells at the periphery of the inlet cross-section
        // whose cell centers are too close to wall face centers.
        double wall_clearance_distance = params_.agent.wall_clearance_distance;
        if (params_.agent.wall_clearance_layers > 0 &&
            params_.agent.wall_clearance_layer_thickness > 0.0) {
            wall_clearance_distance = std::max(
                wall_clearance_distance,
                static_cast<double>(params_.agent.wall_clearance_layers) *
                    params_.agent.wall_clearance_layer_thickness);
        }
        if (params_.agent.enable_wall_clearance && wall_clearance_distance > 0.0) {
            const size_t before = candidates.size();
            candidates = FilterWallAdjacentPositions(
                candidates, wall_meshes_, wall_clearance_distance);
            const size_t removed = before - candidates.size();
            bdm::Log::Info("PreciceAdapter",
                "Wall-clearance filter (d=", wall_clearance_distance, " m): ",
                "kept ", candidates.size(), " / ", before,
                " inlet volume cells (removed ", removed, " wall-adjacent)");
        }

        // Optionally reduce a 3-D inlet cross-section to one or more central
        // lines. Each configured axis is the normal of a retained middle
        // slice. For example, axes y;z form an orthogonal cross: y=midpoint
        // gives the z-directed arm and z=midpoint gives the y-directed arm.
        // The result is the UNION of the slices, not their intersection.
        // Preserve original mapped owner-cell coordinates so every agent is
        // still seeded at a real CFD volume-cell center.
        std::vector<std::string> midplane_axes =
            params_.agent.inlet_midplane_axes;
        if (midplane_axes.empty() &&
            params_.agent.inlet_midplane_axis != "none") {
            midplane_axes.push_back(params_.agent.inlet_midplane_axis);
        }
        if (!candidates.empty() && !midplane_axes.empty()) {
            const auto unsliced_candidates = candidates;
            std::vector<std::array<double, 3>> middle_slice_candidates;
            middle_slice_candidates.reserve(
                unsliced_candidates.size() * midplane_axes.size());

            size_t valid_axis_count = 0;
            for (const auto& axis_name : midplane_axes) {
                int axis = -1;
                if (axis_name == "x" || axis_name == "X") axis = 0;
                if (axis_name == "y" || axis_name == "Y") axis = 1;
                if (axis_name == "z" || axis_name == "Z") axis = 2;

                if (axis < 0) {
                    bdm::Log::Warning(
                        "PreciceAdapter",
                        "Ignoring invalid inlet midplane axis '", axis_name,
                        "' (expected x, y, or z)");
                    continue;
                }
                ++valid_axis_count;

                const double midpoint = 0.5 * (ib.min[axis] + ib.max[axis]);
                double nearest_distance = std::numeric_limits<double>::max();
                for (const auto& pos : unsliced_candidates) {
                    nearest_distance = std::min(
                        nearest_distance, std::abs(pos[axis] - midpoint));
                }

                const double axis_extent = std::abs(ib.max[axis] - ib.min[axis]);
                const double tolerance = std::max(1e-12, axis_extent * 1e-9);
                double slice_half_width =
                    params_.agent.inlet_midplane_tolerance;
                if (slice_half_width <= 0.0) {
                    // A surface candidate cloud has roughly sqrt(N) samples
                    // across either transverse direction. Use that spacing as
                    // an automatic one-layer slice width.
                    slice_half_width = axis_extent /
                        std::sqrt(static_cast<double>(unsliced_candidates.size()));
                }
                slice_half_width = std::max(slice_half_width, nearest_distance);

                size_t arm_count = 0;
                for (const auto& pos : unsliced_candidates) {
                    if (std::abs(pos[axis] - midpoint) >
                        slice_half_width + tolerance) {
                        continue;
                    }
                    if (std::find(middle_slice_candidates.begin(),
                                  middle_slice_candidates.end(), pos) ==
                        middle_slice_candidates.end()) {
                        middle_slice_candidates.push_back(pos);
                    }
                    ++arm_count;
                }

                bdm::Log::Info(
                    "PreciceAdapter", "Inlet middle-slice arm (axis=",
                    axis_name, ", coordinate=", midpoint, ", half-width=",
                    slice_half_width, " m): selected ", arm_count, " / ",
                    unsliced_candidates.size(), " candidates before cross-arm deduplication");
            }

            if (valid_axis_count > 0 && !middle_slice_candidates.empty()) {
                bdm::Log::Info(
                    "PreciceAdapter", "Inlet cross/multi-slice union: kept ",
                    middle_slice_candidates.size(), " unique / ",
                    unsliced_candidates.size(), " inlet candidates across ",
                    valid_axis_count, " arm(s)");
                candidates = std::move(middle_slice_candidates);
            }
        }

        double occupied_tolerance = params_.agent.inlet_occupied_tolerance;
        if (occupied_tolerance == 0.0) {
            occupied_tolerance = 0.51 * params_.agent.default_diameter;
        }
        if (params_.agent.inlet_filter_occupied &&
            occupied_tolerance > 0.0 && rm) {
            const double occupied_tol2 = occupied_tolerance * occupied_tolerance;
            std::vector<std::array<double, 3>> free_candidates;
            free_candidates.reserve(candidates.size());

            for (const auto& pos : candidates) {
                bool occupied = false;
                rm->ForEachAgent([&](bdm::Agent* existing) {
                    auto* cell = dynamic_cast<CellType*>(existing);
                    if (!cell) return true;
                    const auto& p = cell->GetPosition();
                    const double dx = p[0] - pos[0];
                    const double dy = p[1] - pos[1];
                    const double dz = p[2] - pos[2];
                    if (dx * dx + dy * dy + dz * dz <= occupied_tol2) {
                        occupied = true;
                        return false;
                    }
                    return true;
                });
                if (!occupied) {
                    free_candidates.push_back(pos);
                }
            }

            const size_t skipped = candidates.size() - free_candidates.size();
            if (skipped > 0) {
                bdm::Log::Info("PreciceAdapter",
                    "Inlet occupancy filter (d=", occupied_tolerance,
                    " m): skipped ", skipped,
                    " already occupied mapped volume cell centers");
            }
            candidates = std::move(free_candidates);
        } else if (!params_.agent.inlet_filter_occupied) {
            bdm::Log::Info(
                "PreciceAdapter",
                "Inlet occupancy filter disabled: reusing the same validated "
                "owner cell centers for a fixed-size injection batch");
        }

        if (!candidates.empty()) {
            std::array<double, 3> min_pos = candidates.front();
            std::array<double, 3> max_pos = candidates.front();
            for (const auto& pos : candidates) {
                for (int d = 0; d < 3; ++d) {
                    min_pos[d] = std::min(min_pos[d], pos[d]);
                    max_pos[d] = std::max(max_pos[d], pos[d]);
                }
            }
            bdm::Log::Info("PreciceAdapter",
                "Final inlet seeding candidates: ", candidates.size(),
                " cell centers, bounds min=(", min_pos[0], ", ", min_pos[1],
                ", ", min_pos[2], ") max=(", max_pos[0], ", ", max_pos[1],
                ", ", max_pos[2], ")");
        }

        // Optional cap for spawn calls: create up to max_agents.
        // Random mode samples without replacement.  Deterministic mode uses
        // orientation-independent farthest-point sampling, with stable
        // lexicographic tie-breaking.  It therefore returns the exact same
        // mapped owner-cell set on every injection while retaining good
        // coverage of the validated inlet region.
        std::vector<std::array<double, 3>> selected_candidates;
        const bool use_cap = (max_agents > 0);
        if (use_cap && !candidates.empty()) {
            selected_candidates = candidates;
            const size_t target = std::min(max_agents, candidates.size());
            if (params_.agent.inlet_random_sampling) {
                if (target < selected_candidates.size()) {
                    static thread_local std::mt19937 rng(std::random_device{}());
                    std::shuffle(selected_candidates.begin(), selected_candidates.end(), rng);
                    selected_candidates.resize(target);
                }
            } else if (target < selected_candidates.size()) {
                std::sort(selected_candidates.begin(), selected_candidates.end(),
                          [](const auto& lhs, const auto& rhs) {
                              if (lhs[0] != rhs[0]) return lhs[0] < rhs[0];
                              if (lhs[1] != rhs[1]) return lhs[1] < rhs[1];
                              return lhs[2] < rhs[2];
                          });

                std::array<double, 3> center = {{0.0, 0.0, 0.0}};
                for (const auto& pos : selected_candidates) {
                    center[0] += pos[0];
                    center[1] += pos[1];
                    center[2] += pos[2];
                }
                const double inv_count =
                    1.0 / static_cast<double>(selected_candidates.size());
                center[0] *= inv_count;
                center[1] *= inv_count;
                center[2] *= inv_count;

                auto distance2 = [](const auto& a, const auto& b) {
                    const double dx = a[0] - b[0];
                    const double dy = a[1] - b[1];
                    const double dz = a[2] - b[2];
                    return dx * dx + dy * dy + dz * dz;
                };

                size_t next = 0;
                double nearest_center_d2 =
                    distance2(selected_candidates.front(), center);
                for (size_t i = 1; i < selected_candidates.size(); ++i) {
                    const double d2 = distance2(selected_candidates[i], center);
                    if (d2 < nearest_center_d2) {
                        nearest_center_d2 = d2;
                        next = i;
                    }
                }

                std::vector<bool> chosen(selected_candidates.size(), false);
                std::vector<double> nearest_selected_d2(
                    selected_candidates.size(),
                    std::numeric_limits<double>::max());
                std::vector<std::array<double, 3>> deterministic_selection;
                deterministic_selection.reserve(target);

                for (size_t selection = 0; selection < target; ++selection) {
                    chosen[next] = true;
                    deterministic_selection.push_back(selected_candidates[next]);

                    size_t farthest = 0;
                    double farthest_d2 = -1.0;
                    for (size_t i = 0; i < selected_candidates.size(); ++i) {
                        if (chosen[i]) continue;
                        nearest_selected_d2[i] = std::min(
                            nearest_selected_d2[i],
                            distance2(selected_candidates[i],
                                      selected_candidates[next]));
                        // Strict comparison intentionally preserves the
                        // lexicographically first point on exact ties.
                        if (nearest_selected_d2[i] > farthest_d2) {
                            farthest_d2 = nearest_selected_d2[i];
                            farthest = i;
                        }
                    }
                    next = farthest;
                }
                selected_candidates = std::move(deterministic_selection);
                bdm::Log::Info("PreciceAdapter",
                    "Inlet random sampling disabled (agent.inlet_random_sampling=false): "
                    "deterministically selected the same ", selected_candidates.size(),
                    " mapped owner cells from ", candidates.size(),
                    " available candidates");
            } else {
                bdm::Log::Info("PreciceAdapter",
                    "Inlet random sampling disabled (agent.inlet_random_sampling=false): "
                    "using all ", selected_candidates.size(),
                    " mapped owner cells deterministically");
            }
        } else {
            selected_candidates = candidates;
        }

        // Enforce minimum spacing among agents generated in this batch, too,
        // unless exact one-agent-per-selected-cell injection was requested.
        // Reuse inlet_occupied_tolerance as the user-defined spacing control.
        if (params_.agent.inlet_filter_occupied &&
            occupied_tolerance > 0.0 && !selected_candidates.empty()) {
            const double occupied_tol2 = occupied_tolerance * occupied_tolerance;
            std::vector<std::array<double, 3>> spaced_candidates;
            spaced_candidates.reserve(selected_candidates.size());

            for (const auto& pos : selected_candidates) {
                bool too_close = false;
                for (const auto& accepted : spaced_candidates) {
                    const double dx = accepted[0] - pos[0];
                    const double dy = accepted[1] - pos[1];
                    const double dz = accepted[2] - pos[2];
                    if (dx * dx + dy * dy + dz * dz <= occupied_tol2) {
                        too_close = true;
                        break;
                    }
                }
                if (!too_close) {
                    spaced_candidates.push_back(pos);
                }
            }

            const size_t removed_due_spacing =
                selected_candidates.size() - spaced_candidates.size();
            if (removed_due_spacing > 0) {
                bdm::Log::Info("PreciceAdapter",
                    "Inlet spacing filter (d=", occupied_tolerance,
                    " m): removed ", removed_due_spacing,
                    " candidates that were too close to newly generated neighbors");
            }
            selected_candidates = std::move(spaced_candidates);
        }

        // Create agents at the validated cell center positions
        size_t count = 0;
        for (const auto& pos : selected_candidates) {
            auto* agent = agent_factory(pos);
            if (agent) {
                rm->AddAgent(agent);
                ++count;
            }
        }

        if (use_cap) {
            bdm::Log::Info("PreciceAdapter", "Created ", count,
                " capped inlet agents from '", mesh_name, "' (requested=", max_agents,
                ", available=", candidates.size(), ")");
        } else {
            bdm::Log::Info("PreciceAdapter", "Created ", count,
                " agents at inlet volume cell centers from '", mesh_name, "'");
        }

        // Log first few positions for diagnostics
        size_t diag_count = std::min(count, size_t(5));
        for (size_t i = 0; i < diag_count && i < selected_candidates.size(); ++i) {
            bdm::Log::Info("PreciceAdapter", "  Agent[", i, "] cell center: (",
                selected_candidates[i][0], ", ", selected_candidates[i][1], ", ", selected_candidates[i][2], ")");
        }

        return count;
    }

    /**
     * @brief Create agents at inlet using default mesh name from params
     */
    size_t CreateAgentsAtInletMesh(
        std::function<CellType*(const std::array<double, 3>&)> agent_factory,
        double offset_into_domain = 0.0,
        size_t max_agents = 0)
    {
        return CreateAgentsAtInletMesh(params_.mesh.inlet_mesh, agent_factory, offset_into_domain, max_agents);
    }

    /**
     * @brief Get inlet mesh positions for manual agent creation
     *
     * Returns face centers from the inlet mesh. User can then create agents
     * with custom logic (e.g., random sampling, specific patterns).
     *
     * @param mesh_name Name of the inlet surface mesh
     * @return Vector of 3D positions (face centers)
     */
    std::vector<std::array<double, 3>> GetInletMeshPositions(const std::string& mesh_name = "") const {
        const std::string& target = mesh_name.empty() ? params_.mesh.inlet_mesh : mesh_name;
        std::vector<std::array<double, 3>> positions;

        const MeshData* mesh_ptr = nullptr;
        auto inlet_it = inlet_meshes_.find(target);
        if (inlet_it != inlet_meshes_.end() && inlet_it->second.is_initialized) {
            mesh_ptr = &inlet_it->second;
        }

        if (!mesh_ptr || !mesh_ptr->is_initialized) {
            bdm::Log::Warning("PreciceAdapter", "Inlet mesh '", target, "' not available");
            return positions;
        }

        positions.reserve(mesh_ptr->NumVertices());
        for (size_t i = 0; i < mesh_ptr->NumVertices(); ++i) {
            positions.push_back(mesh_ptr->GetVertex(i));
        }
        return positions;
    }

    /**
     * @brief Get outlet mesh positions for proximity checking
     *
     * @param mesh_name Name of the outlet surface mesh
     * @return Vector of 3D positions (face centers)
     */
    std::vector<std::array<double, 3>> GetOutletMeshPositions(const std::string& mesh_name = "") const {
        const std::string& target = mesh_name.empty() ? params_.mesh.outlet_mesh : mesh_name;
        std::vector<std::array<double, 3>> positions;

        const MeshData* mesh_ptr = nullptr;
        auto outlet_it = outlet_meshes_.find(target);
        if (outlet_it != outlet_meshes_.end() && outlet_it->second.is_initialized) {
            mesh_ptr = &outlet_it->second;
        }

        if (!mesh_ptr || !mesh_ptr->is_initialized) {
            bdm::Log::Warning("PreciceAdapter", "Outlet mesh '", target, "' not available");
            return positions;
        }

        positions.reserve(mesh_ptr->NumVertices());
        for (size_t i = 0; i < mesh_ptr->NumVertices(); ++i) {
            positions.push_back(mesh_ptr->GetVertex(i));
        }
        return positions;
    }

    /**
     * @brief Check if an agent is near the outlet and should be removed
     *
     * Computes minimum distance from agent position to all outlet face centers.
     *
     * @param position Agent position
     * @param threshold Distance threshold for removal [m]
     * @param mesh_name Optional specific outlet mesh name
     * @return true if agent should be removed (within threshold distance of outlet)
     */
    bool IsNearOutlet(const std::array<double, 3>& position, double threshold,
                      const std::string& mesh_name = "") const {
        const std::string& target = mesh_name.empty() ? params_.mesh.outlet_mesh : mesh_name;

        const MeshData* mesh_ptr = nullptr;
        auto outlet_it = outlet_meshes_.find(target);
        if (outlet_it != outlet_meshes_.end() && outlet_it->second.is_initialized) {
            mesh_ptr = &outlet_it->second;
        }

        if (!mesh_ptr || !mesh_ptr->is_initialized) {
            return false;  // No outlet mesh - can't remove
        }

        // Check distance to each outlet face center
        double threshold_sq = threshold * threshold;
        for (size_t i = 0; i < mesh_ptr->NumVertices(); ++i) {
            auto face_center = mesh_ptr->GetVertex(i);
            double dx = position[0] - face_center[0];
            double dy = position[1] - face_center[1];
            double dz = position[2] - face_center[2];
            double dist_sq = dx*dx + dy*dy + dz*dz;

            if (dist_sq <= threshold_sq) {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Remove agents that are near the outlet
     *
     * Iterates through all dynamic agents and removes those within the
     * specified distance of any outlet face center.
     *
     * @param threshold Distance threshold for removal [m]
     * @param mesh_name Optional specific outlet mesh name
     * @return Number of agents removed
     */
    size_t RemoveAgentsNearOutlet(double threshold, const std::string& mesh_name = "") {
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return 0;

        const std::string& target = mesh_name.empty() ? params_.mesh.outlet_mesh : mesh_name;

        // Collect UIDs of agents to remove (can't modify during iteration)
        std::vector<bdm::AgentUid> to_remove;

        rm->ForEachAgent([&](bdm::Agent* agent) {
            auto* cell = dynamic_cast<CellType*>(agent);
            if (!cell || cell->IsStationary()) return;

            const auto& pos = cell->GetPosition();
            std::array<double, 3> pos_arr = {pos[0], pos[1], pos[2]};

            if (IsNearOutlet(pos_arr, threshold, target)) {
                to_remove.push_back(cell->GetUid());
            }
        });

        // Remove collected agents
        for (const auto& uid : to_remove) {
            rm->RemoveAgent(uid);
        }

        if (!to_remove.empty()) {
            bdm::Log::Info("PreciceAdapter", "Removed ", to_remove.size(),
                          " agents near outlet '", target, "' (threshold=", threshold, " m)");
        }

        return to_remove.size();
    }

    /**
     * @brief Get inlet mesh data (for advanced use)
     */
    const std::map<std::string, MeshData>& GetInletMeshes() const {
        return inlet_meshes_;
    }

    /**
     * @brief Get outlet mesh data (for advanced use)
     */
    const std::map<std::string, MeshData>& GetOutletMeshes() const {
        return outlet_meshes_;
    }

    /**
     * @brief Check if inlet mesh is available
     */
    bool HasInletMesh(const std::string& mesh_name = "") const {
        const std::string& target = mesh_name.empty() ? params_.mesh.inlet_mesh : mesh_name;
        auto it = inlet_meshes_.find(target);
        return it != inlet_meshes_.end() && it->second.is_initialized && it->second.NumVertices() > 0;
    }

    /**
     * @brief Check if outlet mesh is available
     */
    bool HasOutletMesh(const std::string& mesh_name = "") const {
        const std::string& target = mesh_name.empty() ? params_.mesh.outlet_mesh : mesh_name;
        auto it = outlet_meshes_.find(target);
        return it != outlet_meshes_.end() && it->second.is_initialized && it->second.NumVertices() > 0;
    }

    // =========================================================================
    // Two-Way Coupling: Scalar Sink Collection (Public API)
    // =========================================================================

    /**
     * @brief Collect cellular uptake sinks at the nearest finite-volume cells.
     */
    void CollectSourceTerms() {
        if (!params_.coupling.two_way_coupling) return;

        ResetSourceTerms();

        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;

        // Get the FVM mesh for particle-to-cell mapping
        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto mesh_it = volume_meshes_.find(mesh_name);
        if (mesh_it == volume_meshes_.end() || !mesh_it->second.is_initialized) {
            bdm::Log::Warning("PreciceAdapter", "Cannot collect source terms: mesh not initialized");
            return;
        }

        const auto& mesh = mesh_it->second;
        const bool collect_oxygen_sink =
            params_.coupling.enable_oxygen_sink &&
            IsDataNameEnabled(params_.precice.oxygen_sink_data_name);
        const bool collect_h2o2_sink =
            params_.coupling.enable_h2o2_sink &&
            IsDataNameEnabled(params_.precice.h2o2_sink_data_name);
        const bool collect_no2_sink =
            params_.coupling.enable_no2_sink &&
            IsDataNameEnabled(params_.precice.no2_sink_data_name);

        // Accumulate scalar uptake from all tissue agents.
        rm->ForEachAgent([&](bdm::Agent* agent) {
            auto* cell = dynamic_cast<CellType*>(agent);
            if (!cell) return;

            const auto& position = cell->GetPosition();

            // Find nearest mesh vertex (= OpenFOAM cell center).
            // Cache hit: O(1) for stationary agents (all of o2Cube) and for
            // mobile agents whose cell has not moved since last
            // UpdateDynamicMapping() call.  Cache miss: new cell (post-
            // division) or first call → brute-force once, then cached.
            int nearest_vertex = -1;
            {
                const auto uid = cell->GetUid();
                auto cache_it = agent_fvm_cell_cache_.find(uid);
                if (cache_it != agent_fvm_cell_cache_.end()) {
                    nearest_vertex = cache_it->second;  // O(1) cache hit
                } else {
                    nearest_vertex = FindNearestMeshVertex(mesh, position);
                    if (nearest_vertex >= 0) {
                        agent_fvm_cell_cache_[uid] = nearest_vertex;
                    }
                }
            }

            if (nearest_vertex >= 0) {
                if (collect_oxygen_sink) {
                    double uptake = 0.0;

                    // Michaelis-Menten kinetics: uptake = Vmax * c / (Km + c)
                    // Preferred for cancer spheroid O2 consumption; saturates at high [O2]
                    // and avoids the linear model's unrealistic near-zero uptake at hypoxia.
                    //
                    // IMPORTANT: oxygen_uptake_ on the cell is treated as a cached
                    // diagnostic/output value, not as authoritative input. Reusing
                    // the cached value here would freeze O2Sink across coupling
                    // iterations and produce zero residual vectors in preCICE.
                    if (params_.coupling.oxygen_uptake_km > 0.0
                            && params_.coupling.oxygen_uptake_vmax > 0.0) {
                        double oxygen = 0.0;
                        if constexpr (has_oxygen_getter_v<CellType>) {
                            oxygen = cell->GetOxygen();
                        } else if constexpr (has_temperature_v<CellType>) {
                            oxygen = cell->GetTemperature();
                        }
                        oxygen = std::max(0.0, oxygen);
                        uptake = params_.coupling.oxygen_uptake_vmax * oxygen
                                 / (params_.coupling.oxygen_uptake_km + oxygen);
                    }

                    // Fallback: first-order (linear) uptake model: uptake = k * c
                    if (uptake <= 0.0 && params_.coupling.oxygen_uptake_rate > 0.0) {
                        double oxygen = 0.0;
                        if constexpr (has_oxygen_getter_v<CellType>) {
                            oxygen = cell->GetOxygen();
                        } else if constexpr (has_temperature_v<CellType>) {
                            oxygen = cell->GetTemperature();
                        }

                        oxygen = std::max(0.0, oxygen);
                        uptake = params_.coupling.oxygen_uptake_rate * oxygen;
                    }

                    if (params_.coupling.oxygen_uptake_max > 0.0) {
                        uptake = std::min(uptake, params_.coupling.oxygen_uptake_max);
                    }
                    uptake = std::max(0.0, uptake);

                    if constexpr (has_oxygen_uptake_setter_v<CellType>) {
                        cell->SetOxygenUptake(uptake);
                    }

                    TransferOxygenSink(nearest_vertex, uptake);
                }

                if (collect_no2_sink) {
                    double no2 = 0.0;
                    if constexpr (has_no2_getter_v<CellType>) {
                        no2 = std::max(0.0, cell->GetNO2());
                    }

                    double uptake = params_.coupling.no2_uptake_rate * no2;
                    if (params_.coupling.no2_uptake_max > 0.0) {
                        uptake = std::min(uptake, params_.coupling.no2_uptake_max);
                    }
                    TransferNO2Sink(nearest_vertex, std::max(0.0, uptake));
                }

                if (collect_h2o2_sink) {
                    double h2o2 = 0.0;
                    if constexpr (has_h2o2_getter_v<CellType>) {
                        h2o2 = std::max(0.0, cell->GetH2O2());
                    }
                    double uptake = params_.coupling.h2o2_uptake_rate * h2o2;
                    if (params_.coupling.h2o2_uptake_max > 0.0) {
                        uptake = std::min(uptake, params_.coupling.h2o2_uptake_max);
                    }
                    uptake = std::max(0.0, uptake);
                    if constexpr (has_h2o2_uptake_setter_v<CellType>) {
                        cell->SetH2O2Uptake(uptake);
                    }
                    TransferH2O2Sink(nearest_vertex, uptake);
                }
            }
        });

        if (collect_oxygen_sink && (coupling_iteration_ < 5 || coupling_iteration_ % 20 == 0)) {
            double total_sink = 0.0;
            double max_sink = 0.0;
            for (const auto& [cell_id, sink] : oxygen_sink_source_) {
                static_cast<void>(cell_id);
                total_sink += sink;
                max_sink = std::max(max_sink, sink);
            }
            bdm::Log::Info("PreciceAdapter", "O2 sink aggregation: activeCells=",
                           oxygen_sink_source_.size(), " total=", total_sink, " max=", max_sink);

            // #region agent log
            std::ostringstream payload;
            payload << "{"
                    << "\"iteration\":" << coupling_iteration_ << ","
                    << "\"activeCells\":" << oxygen_sink_source_.size() << ","
                    << "\"totalSink\":" << total_sink << ","
                    << "\"maxSink\":" << max_sink
                    << "}";
            DebugNdjsonLog("baseline", "H1", "bdm_adapter.h:CollectSourceTerms",
                           "oxygen_sink_aggregation", payload.str());
            // #endregion
        }

        if (collect_no2_sink && (coupling_iteration_ < 5 || coupling_iteration_ % 20 == 0)) {
            double total_sink = 0.0;
            double max_sink = 0.0;
            for (const auto& [cell_id, sink] : no2_sink_source_) {
                static_cast<void>(cell_id);
                total_sink += sink;
                max_sink = std::max(max_sink, sink);
            }
            bdm::Log::Info("PreciceAdapter", "NO2 sink aggregation: activeCells=",
                           no2_sink_source_.size(), " total=", total_sink,
                           " max=", max_sink);
        }

    }

protected:
    // =========================================================================
    // Internal Implementation
    // =========================================================================

    // #region agent log
    void DebugNdjsonLog(const std::string& run_id,
                        const std::string& hypothesis_id,
                        const std::string& location,
                        const std::string& message,
                        const std::string& data_json) const {
        static_cast<void>(run_id);
        static_cast<void>(hypothesis_id);
        static_cast<void>(location);
        static_cast<void>(message);
        static_cast<void>(data_json);
    }
    // #endregion

    void SetupProfiling() {
        if (params_.precice.enable_profiling) {
            setenv(params_.precice.profiling_variable.c_str(),
                   params_.precice.profiling_value.c_str(), 1);
            setenv(params_.precice.profiling_mode_variable.c_str(),
                   params_.precice.profiling_mode_value.c_str(), 1);
        }
    }

    void SetWideAccessRegion(const std::string& mesh_name, int dimensions) {
        // Always use 3D bounding box - preCICE expects [xmin, ymin, zmin, xmax, ymax, zmax]
        // Even for 2D meshes, we need the full 3D box to capture all vertices
        constexpr double kInitialAccessExtent = 1.0e6;
        // preCICE expects the layout [xmin, xmax, ymin, ymax, zmin, zmax]
        std::vector<double> bounding_box = {
            -kInitialAccessExtent,  kInitialAccessExtent,  // X range
            -kInitialAccessExtent,  kInitialAccessExtent,  // Y range
            -kInitialAccessExtent,  kInitialAccessExtent   // Z range
        };

        bdm::Log::Info("PreciceAdapter", "Setting access region for '", mesh_name,
                      "' with dimensions=", dimensions, ", using 3D box: [",
                      bounding_box[0], ",", bounding_box[2], ",", bounding_box[4], "] to [",
                      bounding_box[1], ",", bounding_box[3], ",", bounding_box[5], "]");

        try {
            participant_->setMeshAccessRegion(mesh_name,
                precice::span<const double>(bounding_box.data(), bounding_box.size()));
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter", "Failed to set access region for '",
                             mesh_name, "': ", e.what());
        }
    }

    void InitializeVolumeMeshes() {
        for (const auto& mesh_name : params_.mesh.default_volume_meshes) {
            MeshData mesh;
            mesh.name = mesh_name;

            // Get mesh vertices from preCICE
            try {
                int n_vertices = participant_->getMeshVertexSize(mesh_name);
                if (n_vertices > 0) {
                    mesh.vertex_ids.resize(n_vertices);
                    mesh.coordinates.resize(n_vertices * 3);

                    participant_->getMeshVertexIDsAndCoordinates(
                        mesh_name,
                        precice::span<int>(mesh.vertex_ids.data(), mesh.vertex_ids.size()),
                        precice::span<double>(mesh.coordinates.data(), mesh.coordinates.size())
                    );

                    mesh.ComputeBounds();
                    mesh.is_initialized = true;

                    // Update domain bounds
                    if (volume_meshes_.empty()) {
                        domain_bounds_ = mesh.bounds;
                    } else {
                        domain_bounds_.min[0] = std::min(domain_bounds_.min[0], mesh.bounds.min[0]);
                        domain_bounds_.min[1] = std::min(domain_bounds_.min[1], mesh.bounds.min[1]);
                        domain_bounds_.min[2] = std::min(domain_bounds_.min[2], mesh.bounds.min[2]);
                        domain_bounds_.max[0] = std::max(domain_bounds_.max[0], mesh.bounds.max[0]);
                        domain_bounds_.max[1] = std::max(domain_bounds_.max[1], mesh.bounds.max[1]);
                        domain_bounds_.max[2] = std::max(domain_bounds_.max[2], mesh.bounds.max[2]);
                    }

                    bdm::Log::Info("PreciceAdapter", "Initialized volume mesh '", mesh_name,
                                  "' with ", n_vertices, " vertices");
                }
            } catch (const std::exception& e) {
                bdm::Log::Warning("PreciceAdapter", "Failed to initialize mesh '", mesh_name, "': ", e.what());
            }

            volume_meshes_[mesh_name] = std::move(mesh);
        }
    }

    void InitializeWallMeshes() {
        for (const auto& mesh_name : params_.mesh.default_wall_meshes) {
            WallMeshData mesh;
            mesh.name = mesh_name;

            // Get wall mesh vertices from preCICE
            try {
                int n_vertices = participant_->getMeshVertexSize(mesh_name);
                if (n_vertices > 0) {
                    mesh.vertex_ids.resize(n_vertices);
                    mesh.coordinates.resize(n_vertices * 3);

                    participant_->getMeshVertexIDsAndCoordinates(
                        mesh_name,
                        precice::span<int>(mesh.vertex_ids.data(), mesh.vertex_ids.size()),
                        precice::span<double>(mesh.coordinates.data(), mesh.coordinates.size())
                    );

                    mesh.ComputeBounds();
                    mesh.is_initialized = true;

                    // Seed face-center geometry with the received mesh vertex
                    // coordinates. For a static wall these *are* the face
                    // centres (OpenFOAM registers WallSurfaceMesh vertices at
                    // face centres in Interface.C). This guarantees wall
                    // collision physics and domain-bounds queries remain
                    // correct even when the WallGeometry data field is not
                    // actively populated (e.g. restart-from-steady-state
                    // runs where the initial data exchange silently returns
                    // zeros). For moving walls (lid-driven cavity),
                    // ReadWallMeshData overwrites mesh.geometry with the
                    // current face centres each step.
                    mesh.geometry = mesh.coordinates;

                    bdm::Log::Info("PreciceAdapter", "Initialized wall mesh '", mesh_name,
                                  "' with ", n_vertices, " vertices (geometry seeded from coordinates)");
                }
            } catch (const std::exception& e) {
                bdm::Log::Warning("PreciceAdapter", "Failed to initialize wall mesh '", mesh_name, "': ", e.what());
            }

            wall_meshes_[mesh_name] = std::move(mesh);
        }
    }

    void InitializeBoundaryTriangleMesh() {
        const std::string& mesh_name = params_.mesh.boundary_triangle_mesh;
        if (mesh_name.empty()) return;

        BoundaryTriangleMeshData mesh;
        mesh.name = mesh_name;
        try {
            const int n_vertices = participant_->getMeshVertexSize(mesh_name);
            if (n_vertices <= 0) {
                bdm::Log::Warning("PreciceAdapter",
                    "Exact boundary mesh '", mesh_name, "' is empty");
                return;
            }

            mesh.vertex_ids.resize(n_vertices);
            mesh.coordinates.resize(static_cast<std::size_t>(n_vertices) * 3);
            participant_->getMeshVertexIDsAndCoordinates(
                mesh_name,
                precice::span<int>(mesh.vertex_ids.data(), mesh.vertex_ids.size()),
                precice::span<double>(mesh.coordinates.data(), mesh.coordinates.size()));
            mesh.ComputeBounds();
            mesh.is_initialized = true;

            const auto read_field = [&](const std::string& data_name,
                                        int expected_dim,
                                        std::vector<double>* values,
                                        bool required) {
                try {
                    const int dim =
                        participant_->getDataDimensions(mesh_name, data_name);
                    if (dim != expected_dim) {
                        throw std::runtime_error("dimension mismatch");
                    }
                    values->assign(
                        static_cast<std::size_t>(n_vertices) * expected_dim, 0.0);
                    participant_->readData(
                        mesh_name, data_name,
                        precice::span<const int>(
                            mesh.vertex_ids.data(), mesh.vertex_ids.size()),
                        0.0,
                        precice::span<double>(values->data(), values->size()));
                } catch (const std::exception& e) {
                    values->clear();
                    if (required) {
                        throw std::runtime_error(
                            "required boundary data '" + data_name
                            + "' unavailable: " + e.what());
                    }
                }
            };

            read_field(params_.precice.boundary_vertex0_data_name, 3,
                       &mesh.vertex0, true);
            read_field(params_.precice.boundary_vertex1_data_name, 3,
                       &mesh.vertex1, true);
            read_field(params_.precice.boundary_vertex2_data_name, 3,
                       &mesh.vertex2, true);
            read_field(params_.precice.boundary_kind_data_name, 1,
                       &mesh.kinds, true);
            read_field(params_.precice.boundary_normals_data_name, 3,
                       &mesh.normals, false);
            read_field(params_.precice.boundary_velocity_data_name, 3,
                       &mesh.velocity, false);

            if (mesh.normals.size() != mesh.vertex0.size()) {
                mesh.normals.assign(mesh.vertex0.size(), 0.0);
            }
            if (mesh.velocity.size() != mesh.vertex0.size()) {
                mesh.velocity.assign(mesh.vertex0.size(), 0.0);
            }
            if (!mesh.HasCompleteRecords()) {
                throw std::runtime_error("incomplete triangle records");
            }

            std::vector<BoundaryTriangle> triangles;
            triangles.reserve(mesh.NumVertices());
            for (std::size_t i = 0; i < mesh.NumVertices(); ++i) {
                const std::size_t j = i * 3;
                BoundaryTriangle tri;
                tri.a = {mesh.vertex0[j], mesh.vertex0[j+1], mesh.vertex0[j+2]};
                tri.b = {mesh.vertex1[j], mesh.vertex1[j+1], mesh.vertex1[j+2]};
                tri.c = {mesh.vertex2[j], mesh.vertex2[j+1], mesh.vertex2[j+2]};
                tri.normal = {
                    mesh.normals[j], mesh.normals[j+1], mesh.normals[j+2]};
                tri.velocity = {
                    mesh.velocity[j], mesh.velocity[j+1], mesh.velocity[j+2]};
                tri.kind = mesh.kinds[i] >= 0.5
                    ? BoundaryKind::kWall : BoundaryKind::kPassThrough;
                triangles.push_back(tri);
            }

            auto geometry = std::make_shared<BoundaryGeometry>();
            if (!geometry->Build(std::move(triangles))) {
                throw std::runtime_error("all boundary triangles are degenerate");
            }
            boundary_triangle_mesh_ = std::move(mesh);
            boundary_geometry_ = std::move(geometry);
            bdm::Log::Info("PreciceAdapter", "Initialized exact boundary mesh '",
                mesh_name, "' with ", boundary_geometry_->Size(), " triangles; closed=",
                boundary_geometry_->IsClosed() ? "true" : "false");
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter",
                "Failed to initialize exact boundary mesh '", mesh_name,
                "': ", e.what());
            boundary_geometry_.reset();
        }
    }

    /**
     * @brief Initialize inlet and outlet surface meshes
     *
     * These meshes contain face centers from OpenFOAM inlet/outlet patches.
     * Used for:
     * - Inlet: Agent generation at inlet face positions
     * - Outlet: Agent removal when approaching outlet
     */
    void InitializeInletOutletMeshes() {
        // Initialize inlet meshes
        for (const auto& mesh_name : params_.mesh.default_inlet_meshes) {
            MeshData mesh;
            mesh.name = mesh_name;

            try {
                int n_vertices = participant_->getMeshVertexSize(mesh_name);
                if (n_vertices > 0) {
                    mesh.vertex_ids.resize(n_vertices);
                    mesh.coordinates.resize(n_vertices * 3);

                    participant_->getMeshVertexIDsAndCoordinates(
                        mesh_name,
                        precice::span<int>(mesh.vertex_ids.data(), mesh.vertex_ids.size()),
                        precice::span<double>(mesh.coordinates.data(), mesh.coordinates.size())
                    );

                    mesh.ComputeBounds();
                    mesh.is_initialized = true;
                    ReadInletCellCenterData(mesh_name, mesh, 0.0);
                    ReadInletAreaData(mesh_name, mesh, 0.0);

                    bdm::Log::Info("PreciceAdapter", "Initialized inlet mesh '", mesh_name,
                                  "' with ", n_vertices, " face centers",
                                  mesh.HasMappedCellCenters() ? " and exact owner cell centers" : "",
                                  mesh.HasFaceAreas() ? " and face areas" : "");
                }
            } catch (const std::exception& e) {
                bdm::Log::Warning("PreciceAdapter", "Failed to initialize inlet mesh '", mesh_name, "': ", e.what());
            }

            inlet_meshes_[mesh_name] = std::move(mesh);
        }

        // Initialize outlet meshes
        for (const auto& mesh_name : params_.mesh.default_outlet_meshes) {
            MeshData mesh;
            mesh.name = mesh_name;

            try {
                int n_vertices = participant_->getMeshVertexSize(mesh_name);
                if (n_vertices > 0) {
                    mesh.vertex_ids.resize(n_vertices);
                    mesh.coordinates.resize(n_vertices * 3);

                    participant_->getMeshVertexIDsAndCoordinates(
                        mesh_name,
                        precice::span<int>(mesh.vertex_ids.data(), mesh.vertex_ids.size()),
                        precice::span<double>(mesh.coordinates.data(), mesh.coordinates.size())
                    );

                    mesh.ComputeBounds();
                    mesh.is_initialized = true;

                    bdm::Log::Info("PreciceAdapter", "Initialized outlet mesh '", mesh_name,
                                  "' with ", n_vertices, " face centers");
                }
            } catch (const std::exception& e) {
                bdm::Log::Warning("PreciceAdapter", "Failed to initialize outlet mesh '", mesh_name, "': ", e.what());
            }

            outlet_meshes_[mesh_name] = std::move(mesh);
        }

        // Log summary
        size_t total_inlet = 0, total_outlet = 0;
        for (const auto& [name, mesh] : inlet_meshes_) {
            if (mesh.is_initialized) total_inlet += mesh.NumVertices();
        }
        for (const auto& [name, mesh] : outlet_meshes_) {
            if (mesh.is_initialized) total_outlet += mesh.NumVertices();
        }
        if (total_inlet > 0 || total_outlet > 0) {
            bdm::Log::Info("PreciceAdapter", "Inlet/Outlet meshes: ", total_inlet, " inlet faces, ", total_outlet, " outlet faces");
        }
    }

    static bool IsDataNameEnabled(const std::string& name) {
        if (name.empty()) {
            return false;
        }
        std::string lowered = name;
        std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        return lowered != "none";
    }

    static std::string ToLower(std::string value) {
        std::transform(value.begin(), value.end(), value.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        return value;
    }

    bool IsSteadyFlowMode() const {
        return ToLower(params_.coupling.flow_mode) == "steady";
    }

    int NormalizeReadInterval(int interval) const {
        return interval < 0 ? 1 : interval;
    }

    bool ShouldReadWithInterval(int configured_interval, bool already_read_once) const {
        const int interval = NormalizeReadInterval(configured_interval);
        if (interval == 0) {
            return !already_read_once;
        }
        if (fetch_input_call_count_ == 0) {
            return true;
        }
        return ((fetch_input_call_count_ - 1) % static_cast<std::size_t>(interval)) == 0;
    }

    bool ShouldReadFieldData() const {
        return ShouldReadWithInterval(params_.coupling.field_read_interval,
                                      field_data_read_once_);
    }

    bool ShouldReadWallData() const {
        if (params_.coupling.static_mesh) {
            return !wall_data_read_once_;
        }
        return ShouldReadWithInterval(params_.coupling.wall_read_interval,
                                      wall_data_read_once_);
    }

    bool ShouldUseCachedFieldResampling() const {
        const int interval = NormalizeReadInterval(params_.coupling.field_read_interval);
        return interval == 0 || interval > 1;
    }

    bool HasCachedFieldSnapshot() const {
        return cached_velocity_valid_ ||
               cached_pressure_valid_ ||
               cached_pressure_gradient_valid_ ||
               cached_fluid_acceleration_valid_ ||
               cached_fluid_vorticity_valid_ ||
               cached_temperature_valid_ ||
               cached_oxygen_valid_;
    }

    void CacheCurrentFieldSnapshot(double relative_time) {
        if (!participant_) {
            return;
        }

        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto mesh_it = volume_meshes_.find(mesh_name);
        if (mesh_it == volume_meshes_.end() || !mesh_it->second.is_initialized) {
            return;
        }

        const auto& mesh = mesh_it->second;
        const std::size_t n_vertices = mesh.vertex_ids.size();
        if (n_vertices == 0) {
            return;
        }

        auto ids_span =
            precice::span<const int>(mesh.vertex_ids.data(), mesh.vertex_ids.size());

        bool any_cached = false;

        auto read_vector_snapshot = [&](bool available, const std::string& data_name,
                                        std::vector<double>& target, bool& valid) {
            if (!available || !IsDataNameEnabled(data_name)) {
                target.clear();
                valid = false;
                return;
            }

            std::vector<double> values(n_vertices * 3, 0.0);
            auto values_span = precice::span<double>(values.data(), values.size());
            try {
                participant_->readData(mesh_name, data_name, ids_span, relative_time,
                                       values_span);
                target = std::move(values);
                valid = true;
                any_cached = true;
            } catch (const std::exception& e) {
                if (!valid) {
                    bdm::Log::Warning("PreciceAdapter",
                                      "Failed to cache vector field '", data_name,
                                      "' for frozen re-sampling: ", e.what());
                }
            }
        };

        auto read_scalar_snapshot = [&](bool available, const std::string& data_name,
                                        std::vector<double>& target, bool& valid) {
            if (!available || !IsDataNameEnabled(data_name)) {
                target.clear();
                valid = false;
                return;
            }

            std::vector<double> values(n_vertices, 0.0);
            auto values_span = precice::span<double>(values.data(), values.size());
            try {
                participant_->readData(mesh_name, data_name, ids_span, relative_time,
                                       values_span);
                target = std::move(values);
                valid = true;
                any_cached = true;
            } catch (const std::exception& e) {
                if (!valid) {
                    bdm::Log::Warning("PreciceAdapter",
                                      "Failed to cache scalar field '", data_name,
                                      "' for frozen re-sampling: ", e.what());
                }
            }
        };

        read_vector_snapshot(velocity_available_, params_.precice.velocity_data_name,
                             cached_velocity_field_, cached_velocity_valid_);
        read_scalar_snapshot(pressure_available_, params_.precice.pressure_data_name,
                             cached_pressure_field_, cached_pressure_valid_);
        read_vector_snapshot(pressure_gradient_available_, params_.precice.pressure_gradient_data_name,
                             cached_pressure_gradient_field_, cached_pressure_gradient_valid_);
        read_vector_snapshot(fluid_acceleration_available_, params_.precice.fluid_acceleration_data_name,
                             cached_fluid_acceleration_field_, cached_fluid_acceleration_valid_);
        read_vector_snapshot(fluid_vorticity_available_, params_.precice.fluid_vorticity_data_name,
                             cached_fluid_vorticity_field_, cached_fluid_vorticity_valid_);
        read_scalar_snapshot(temperature_available_, params_.precice.temperature_data_name,
                             cached_temperature_field_, cached_temperature_valid_);
        read_scalar_snapshot(oxygen_available_, params_.precice.oxygen_data_name,
                             cached_oxygen_field_, cached_oxygen_valid_);

        if (any_cached) {
            cached_field_vertex_count_ = n_vertices;
            cached_field_notice_emitted_ = false;
        }
    }

    bool ApplyCachedFieldDataToAgents() {
        if (!HasCachedFieldSnapshot()) {
            return false;
        }

        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto mesh_it = volume_meshes_.find(mesh_name);
        if (mesh_it == volume_meshes_.end() || !mesh_it->second.is_initialized) {
            return false;
        }

        const auto& mesh = mesh_it->second;
        if (mesh.vertex_ids.empty()) {
            return false;
        }

        if (cached_field_vertex_count_ > 0 &&
            cached_field_vertex_count_ != mesh.vertex_ids.size()) {
            if (!cached_field_notice_emitted_) {
                bdm::Log::Warning("PreciceAdapter",
                                  "Cannot apply cached fields: mesh vertex count changed "
                                  "from ",
                                  cached_field_vertex_count_, " to ",
                                  mesh.vertex_ids.size());
                cached_field_notice_emitted_ = true;
            }
            return false;
        }

        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) {
            return false;
        }

        bool applied = false;
        rm->ForEachAgent([&](bdm::Agent* agent) {
            auto* cell = dynamic_cast<CellType*>(agent);
            if (!cell) {
                return;
            }

            const auto& p = cell->GetPosition();
            const int nearest_index = FindNearestMeshVertexIndex(mesh, p);
            if (nearest_index < 0) {
                return;
            }

            const std::size_t idx = static_cast<std::size_t>(nearest_index);

            if (cached_velocity_valid_ &&
                cached_velocity_field_.size() >= (idx + 1) * 3 &&
                !cell->IsStationary()) {
                bdm::Real3 vel = {cached_velocity_field_[idx * 3],
                                  cached_velocity_field_[idx * 3 + 1],
                                  cached_velocity_field_[idx * 3 + 2]};
                if (params_.coupling.force_2d) {
                    vel[2] = 0.0;
                }
                cell->SetFluidVelocity(vel);
                applied = true;
            }

            if (cached_pressure_valid_ &&
                cached_pressure_field_.size() > idx) {
                cell->SetPressure(cached_pressure_field_[idx]);
                applied = true;
            }

            if (cached_pressure_gradient_valid_ &&
                cached_pressure_gradient_field_.size() >= (idx + 1) * 3) {
                std::array<double, 3> grad = {
                    cached_pressure_gradient_field_[idx * 3],
                    cached_pressure_gradient_field_[idx * 3 + 1],
                    cached_pressure_gradient_field_[idx * 3 + 2]
                };
                cell->SetPressureGradient(grad);
                applied = true;
            }

            if (cached_fluid_acceleration_valid_ &&
                cached_fluid_acceleration_field_.size() >= (idx + 1) * 3) {
                bdm::Real3 accel = {cached_fluid_acceleration_field_[idx * 3],
                                    cached_fluid_acceleration_field_[idx * 3 + 1],
                                    cached_fluid_acceleration_field_[idx * 3 + 2]};
                if constexpr (has_fluid_acceleration_v<CellType>) {
                    cell->SetFluidAcceleration(accel);
                }
                applied = true;
            }

            if (cached_fluid_vorticity_valid_ &&
                cached_fluid_vorticity_field_.size() >= (idx + 1) * 3) {
                bdm::Real3 vort = {cached_fluid_vorticity_field_[idx * 3],
                                   cached_fluid_vorticity_field_[idx * 3 + 1],
                                   cached_fluid_vorticity_field_[idx * 3 + 2]};
                if constexpr (has_fluid_vorticity_v<CellType>) {
                    cell->SetFluidVorticity(vort);
                }
                applied = true;
            }

            if (cached_temperature_valid_ &&
                cached_temperature_field_.size() > idx) {
                if constexpr (has_temperature_v<CellType>) {
                    cell->SetTemperature(cached_temperature_field_[idx]);
                }
                applied = true;
            }

            if (cached_oxygen_valid_ &&
                cached_oxygen_field_.size() > idx) {
                if constexpr (has_oxygen_setter_v<CellType>) {
                    cell->SetOxygen(cached_oxygen_field_[idx]);
                } else if constexpr (has_temperature_v<CellType>) {
                    cell->SetTemperature(cached_oxygen_field_[idx]);
                }
                applied = true;
            }
        });

        if (!applied && !cached_field_notice_emitted_) {
            bdm::Log::Warning("PreciceAdapter",
                              "Cached field snapshot exists but could not be applied to agents");
            cached_field_notice_emitted_ = true;
        }

        return applied;
    }

    void SetupDataAvailability() {
        // Check what data is available on the primary volume mesh
        const auto& mesh_name = params_.mesh.primary_volume_mesh;

        try {
            velocity_available_ = params_.coupling.enable_velocity &&
                                  IsDataNameEnabled(params_.precice.velocity_data_name) &&
                                  participant_->getDataDimensions(mesh_name, params_.precice.velocity_data_name) == 3;
        } catch (...) { velocity_available_ = false; }

        try {
            pressure_available_ = params_.coupling.enable_pressure &&
                                  IsDataNameEnabled(params_.precice.pressure_data_name) &&
                                  participant_->getDataDimensions(mesh_name, params_.precice.pressure_data_name) == 1;
        } catch (...) { pressure_available_ = false; }

        try {
            pressure_gradient_available_ = params_.coupling.enable_pressure_gradient &&
                                           IsDataNameEnabled(params_.precice.pressure_gradient_data_name) &&
                                           participant_->getDataDimensions(mesh_name, params_.precice.pressure_gradient_data_name) == 3;
        } catch (...) { pressure_gradient_available_ = false; }

        try {
            fluid_acceleration_available_ = IsDataNameEnabled(params_.precice.fluid_acceleration_data_name) &&
                                            participant_->getDataDimensions(mesh_name, params_.precice.fluid_acceleration_data_name) == 3;
        } catch (...) { fluid_acceleration_available_ = false; }

        try {
            fluid_vorticity_available_ = IsDataNameEnabled(params_.precice.fluid_vorticity_data_name) &&
                                         participant_->getDataDimensions(mesh_name, params_.precice.fluid_vorticity_data_name) == 3;
        } catch (...) { fluid_vorticity_available_ = false; }

        if (IsDataNameEnabled(params_.precice.fluid_density_data_name)) {
            try {
                fluid_density_available_ = participant_->getDataDimensions(mesh_name, params_.precice.fluid_density_data_name) == 1;
            } catch (...) {
                fluid_density_available_ = false;
            }
        } else {
            fluid_density_available_ = false;
        }

        if (IsDataNameEnabled(params_.precice.fluid_dynamic_viscosity_data_name)) {
            try {
                fluid_dynamic_viscosity_available_ = participant_->getDataDimensions(mesh_name, params_.precice.fluid_dynamic_viscosity_data_name) == 1;
            } catch (...) {
                fluid_dynamic_viscosity_available_ = false;
            }
        } else {
            fluid_dynamic_viscosity_available_ = false;
        }

        try {
            temperature_available_ = params_.coupling.enable_temperature &&
                                     IsDataNameEnabled(params_.precice.temperature_data_name) &&
                                     participant_->getDataDimensions(mesh_name, params_.precice.temperature_data_name) == 1;
        } catch (...) { temperature_available_ = false; }

        try {
            oxygen_available_ = params_.coupling.enable_oxygen &&
                                IsDataNameEnabled(params_.precice.oxygen_data_name) &&
                                participant_->getDataDimensions(mesh_name, params_.precice.oxygen_data_name) == 1;
        } catch (...) { oxygen_available_ = false; }

        try {
            tke_available_ = params_.coupling.enable_turbulence &&
                             IsDataNameEnabled(params_.precice.tke_data_name) &&
                             participant_->getDataDimensions(mesh_name, params_.precice.tke_data_name) == 1;
        } catch (...) { tke_available_ = false; }

        try {
            epsilon_available_ = params_.coupling.enable_turbulence &&
                                 IsDataNameEnabled(params_.precice.epsilon_data_name) &&
                                 participant_->getDataDimensions(mesh_name, params_.precice.epsilon_data_name) == 1;
        } catch (...) { epsilon_available_ = false; }

        try {
            omega_available_ = params_.coupling.enable_turbulence &&
                               IsDataNameEnabled(params_.precice.omega_data_name) &&
                               participant_->getDataDimensions(mesh_name, params_.precice.omega_data_name) == 1;
        } catch (...) { omega_available_ = false; }

        try {
            nut_available_ = params_.coupling.enable_turbulence &&
                             IsDataNameEnabled(params_.precice.nut_data_name) &&
                             participant_->getDataDimensions(mesh_name, params_.precice.nut_data_name) == 1;
        } catch (...) { nut_available_ = false; }

        bdm::Log::Info("PreciceAdapter", "Data availability - V:", velocity_available_,
                      " P:", pressure_available_, " gradP:", pressure_gradient_available_,
                      " DUcDt:", fluid_acceleration_available_, " curlUc:", fluid_vorticity_available_,
                      " T:", temperature_available_, " O2:", oxygen_available_,
                      " rho:", fluid_density_available_, " mu:", fluid_dynamic_viscosity_available_,
                      " k:", tke_available_, " eps:", epsilon_available_,
                      " omg:", omega_available_, " nut:", nut_available_);
    }

    // -------------------------------------------------------------------------
    // Data Reading Methods
    // -------------------------------------------------------------------------

    // Populates read_coords_buf_ and read_uid_buf_ from the current agent set.
    // Called once per FetchInputData() before any Read*Data calls.
    // If called by a Read*Data function directly (standalone use), builds the
    // buffer lazily and sets read_buf_valid_ = true.
    void RebuildReadBuffer() {
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        const size_t n = rm ? rm->GetNumAgents() : 0;
        read_coords_buf_.clear();
        read_uid_buf_.clear();
        read_coords_buf_.reserve(3 * n);
        read_uid_buf_.reserve(n);
        if (!rm) { read_buf_valid_ = false; return; }
        rm->ForEachAgent([&](bdm::Agent* agent) {
            auto* cell = dynamic_cast<CellType*>(agent);
            if (!cell) return;
            const auto& p = cell->GetPosition();
            read_coords_buf_.push_back(p[0]);
            read_coords_buf_.push_back(p[1]);
            read_coords_buf_.push_back(p[2]);
            read_uid_buf_.push_back(cell->GetUid());
        });
        read_buf_valid_ = true;
    }

    void ReadVelocityData(double relative_time) {
        if (!velocity_available_) return;

        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;

        // Use shared read-buffer built by FetchInputData; fallback to rebuild if called standalone.
        if (!read_buf_valid_) RebuildReadBuffer();
        const auto& coords = read_coords_buf_;
        const auto& uids   = read_uid_buf_;

        // Guard: if no agents exist the JIT read would use an empty coordinate set.
        // This SHOULD NOT happen after correct initialisation (agents are created at inlet
        // cell centres before the first FetchInputData call).  Emit a warning so the
        // problem is immediately visible rather than silently producing zero fields.
        if (uids.empty()) {
            bdm::Log::Warning("PreciceAdapter",
                              "ReadVelocityData: coordinate set is EMPTY (no agents exist). "
                              "mapAndReadData will not be called.  Check that agents are "
                              "created (inlet seeding) before the first FetchInputData call.");
            return;
        }

        const bool is_first_jit_read = (fetch_input_call_count_ == 1);

        if (is_first_jit_read) {
            // ---- First JIT velocity read: emit full diagnostics ----
            // Compute min/max of the query coordinate set.
            std::array<double, 3> cmin = {coords[0], coords[1], coords[2]};
            std::array<double, 3> cmax = {coords[0], coords[1], coords[2]};
            for (std::size_t ci = 0; ci + 2 < coords.size(); ci += 3) {
                cmin[0] = std::min(cmin[0], coords[ci]);
                cmin[1] = std::min(cmin[1], coords[ci + 1]);
                cmin[2] = std::min(cmin[2], coords[ci + 2]);
                cmax[0] = std::max(cmax[0], coords[ci]);
                cmax[1] = std::max(cmax[1], coords[ci + 1]);
                cmax[2] = std::max(cmax[2], coords[ci + 2]);
            }
            bdm::Log::Info("PreciceAdapter",
                           "FIRST JIT READ: mapAndReadData(\"", mesh_name, "\", \"",
                           params_.precice.velocity_data_name, "\")"
                           " coords=", uids.size(),
                           " agents (empty=", uids.empty() ? "true" : "false", ")"
                           " bounds min=(", cmin[0], ",", cmin[1], ",", cmin[2],
                           ") max=(", cmax[0], ",", cmax[1], ",", cmax[2], ")"
                           " relative_time=", relative_time);
            // Log first 5 query coordinates.
            const std::size_t diag_n = std::min(uids.size(), std::size_t(5));
            for (std::size_t di = 0; di < diag_n; ++di) {
                bdm::Log::Info("PreciceAdapter",
                               "  query[", di, "]=(",
                               coords[di * 3], ",",
                               coords[di * 3 + 1], ",",
                               coords[di * 3 + 2], ")");
            }
            // Access region reminder.
            bdm::Log::Info("PreciceAdapter",
                           "  setMeshAccessRegion was \u00b11e6 (all FVMmesh vertices included). "
                           "\"Mapping distance not available due to empty partition\" "
                           "emitted during initialize() is expected for JIT mapping and is benign.");
        }

        // Read velocity using just-in-time mapping
        std::vector<double> values(uids.size() * 3, 0.0);
        auto coords_span = precice::span<const double>(coords.data(), coords.size());
        auto values_span = precice::span<double>(values.data(), values.size());

        try {
            participant_->mapAndReadData(mesh_name, params_.precice.velocity_data_name,
                                        coords_span, relative_time, values_span);

            if (is_first_jit_read) {
                // Spot-check returned values for the first 5 agents.
                const std::size_t diag_n = std::min(uids.size(), std::size_t(5));
                for (std::size_t di = 0; di < diag_n; ++di) {
                    bdm::Log::Info("PreciceAdapter",
                                   "  result[", di, "]=(",
                                   values[di * 3], ",",
                                   values[di * 3 + 1], ",",
                                   values[di * 3 + 2], ")");
                }
            }

            // Apply velocity to agents using the exact same ordering as in coords/uids
            auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
            size_t i = 0;
            int log_first = 0;
            for (const auto& uid : uids) {
                if (i + 2 >= values.size()) {
                    break;
                }

                auto* cell = dynamic_cast<CellType*>(rm->GetAgent(uid));
                if (!cell || cell->IsStationary()) {
                    i += 3;
                    continue;
                }

                bdm::Real3 vel = {values[i], values[i+1], values[i+2]};

                // For 2D simulations: zero out Z-velocity component
                if (params_.coupling.force_2d) {
                    vel[2] = 0.0;
                }

                // Store fluid velocity for drag integration; do NOT overwrite particle kinematic velocity.
                cell->SetFluidVelocity(vel);

                // Diagnostic: log first few mapped velocities with positions
                if (log_first < 3) {
                    const auto& p = cell->GetPosition();
                    bdm::Log::Debug("PreciceAdapter",
                        "VEL MAP agent[", i/3, "] pos=(",
                        p[0], ",", p[1], ",", p[2],
                        ") vel=(", vel[0], ",", vel[1], ",", vel[2], ")");
                    log_first++;
                }

                i += 3;
            }

        } catch (const std::exception& e) {
            bdm::Log::Error("PreciceAdapter", "Error reading velocity data: ", e.what());
        }
    }

    void ReadPressureData(double relative_time) {
        if (!pressure_available_) return;

        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;

        // Use shared read-buffer built by FetchInputData; fallback to rebuild if called standalone.
        if (!read_buf_valid_) RebuildReadBuffer();
        const auto& coords = read_coords_buf_;
        const auto& uids   = read_uid_buf_;

        if (uids.empty()) return;

        std::vector<double> values(uids.size(), 0.0);
        auto coords_span = precice::span<const double>(coords.data(), coords.size());
        auto values_span = precice::span<double>(values.data(), values.size());

        try {
            participant_->mapAndReadData(mesh_name, params_.precice.pressure_data_name,
                                        coords_span, relative_time, values_span);

            size_t i = 0;
            for (const auto& uid : uids) {
                auto* cell = dynamic_cast<CellType*>(rm->GetAgent(uid));
                if (cell) {
                    cell->SetPressure(values[i]);
                }
                ++i;
            }
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter", "Failed to read pressure: ", e.what());
        }
    }

    bool ReadFluidPropertiesInternal(double relative_time, double* density_out, double* mu_out) {
        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto mesh_it = volume_meshes_.find(mesh_name);
        if (mesh_it == volume_meshes_.end() || !mesh_it->second.is_initialized) {
            bdm::Log::Warning("PreciceAdapter", "Cannot read fluid properties: mesh '", mesh_name, "' not initialized");
            return false;
        }

        const auto& mesh = mesh_it->second;
        if (mesh.vertex_ids.empty()) {
            bdm::Log::Warning("PreciceAdapter", "Cannot read fluid properties: no vertices on mesh '", mesh_name, "'");
            return false;
        }

        bool ok = true;

        if (fluid_density_available_) {
            std::vector<double> rho_vals(1, 0.0);
            auto id_span = precice::span<const int>(mesh.vertex_ids.data(), 1);
            auto rho_span = precice::span<double>(rho_vals.data(), rho_vals.size());
            try {
                participant_->readData(mesh_name, params_.precice.fluid_density_data_name,
                                      id_span, relative_time, rho_span);
                if (density_out) { *density_out = rho_vals[0]; }
            } catch (const std::exception& e) {
                ok = false;
                bdm::Log::Warning("PreciceAdapter", "Failed to read fluid density: ", e.what());
            }
        } else if (IsDataNameEnabled(params_.precice.fluid_density_data_name)) {
            ok = false;
            bdm::Log::Warning("PreciceAdapter", "Fluid density data '", params_.precice.fluid_density_data_name,
                             "' not available on mesh '", mesh_name, "'");
        }

        if (fluid_dynamic_viscosity_available_) {
            std::vector<double> mu_vals(1, 0.0);
            auto id_span = precice::span<const int>(mesh.vertex_ids.data(), 1);
            auto mu_span = precice::span<double>(mu_vals.data(), mu_vals.size());
            try {
                participant_->readData(mesh_name, params_.precice.fluid_dynamic_viscosity_data_name,
                                      id_span, relative_time, mu_span);
                if (mu_out) { *mu_out = mu_vals[0]; }
            } catch (const std::exception& e) {
                ok = false;
                bdm::Log::Warning("PreciceAdapter", "Failed to read fluid viscosity: ", e.what());
            }
        } else if (IsDataNameEnabled(params_.precice.fluid_dynamic_viscosity_data_name)) {
            ok = false;
            bdm::Log::Warning("PreciceAdapter", "Fluid viscosity data '", params_.precice.fluid_dynamic_viscosity_data_name,
                             "' not available on mesh '", mesh_name, "'");
        }

        return ok;
    }

    void ReadPressureGradientData(double relative_time) {
        if (!pressure_gradient_available_) return;

        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;

        // Use shared read-buffer built by FetchInputData; fallback to rebuild if called standalone.
        if (!read_buf_valid_) RebuildReadBuffer();
        const auto& coords = read_coords_buf_;
        const auto& uids   = read_uid_buf_;

        if (uids.empty()) return;

        std::vector<double> values(uids.size() * 3, 0.0);
        auto coords_span = precice::span<const double>(coords.data(), coords.size());
        auto values_span = precice::span<double>(values.data(), values.size());

        try {
            participant_->mapAndReadData(mesh_name, params_.precice.pressure_gradient_data_name,
                                        coords_span, relative_time, values_span);

            size_t i = 0;
            for (const auto& uid : uids) {
                auto* cell = dynamic_cast<CellType*>(rm->GetAgent(uid));
                if (cell) {
                    std::array<double, 3> grad = {values[i*3], values[i*3+1], values[i*3+2]};
                    cell->SetPressureGradient(grad);
                }
                ++i;
            }
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter", "Failed to read pressure gradient: ", e.what());
        }
    }

    void ReadTemperatureData(double relative_time) {
        if (!temperature_available_) return;

        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;

        // Use shared read-buffer built by FetchInputData; fallback to rebuild if called standalone.
        if (!read_buf_valid_) RebuildReadBuffer();
        const auto& coords = read_coords_buf_;
        const auto& uids   = read_uid_buf_;

        if (uids.empty()) return;

        std::vector<double> values(uids.size(), 0.0);
        auto coords_span = precice::span<const double>(coords.data(), coords.size());
        auto values_span = precice::span<double>(values.data(), values.size());

        try {
            participant_->mapAndReadData(mesh_name, params_.precice.temperature_data_name,
                                        coords_span, relative_time, values_span);

            size_t i = 0;
            double min_t = std::numeric_limits<double>::max();
            double max_t = std::numeric_limits<double>::lowest();
            double sum_t = 0.0;
            size_t count_t = 0;
            for (const auto& uid : uids) {
                auto* cell = dynamic_cast<CellType*>(rm->GetAgent(uid));
                if (cell) {
                    if constexpr (has_temperature_v<CellType>) {
                        const double t = values[i];
                        cell->SetTemperature(t);
                        min_t = std::min(min_t, t);
                        max_t = std::max(max_t, t);
                        sum_t += t;
                        ++count_t;
                    }
                }
                ++i;
            }

            if (count_t > 0 && (coupling_iteration_ < 5 || coupling_iteration_ % 20 == 0)) {
                const double avg_t = sum_t / static_cast<double>(count_t);
                bdm::Log::Info("PreciceAdapter", "Temperature read stats: count=", count_t,
                               " min=", min_t, " max=", max_t, " avg=", avg_t);
            }
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter", "Failed to read temperature: ", e.what());
        }
    }

    void ReadOxygenData(double relative_time) {
        if (!oxygen_available_) return;

        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;

        // Use shared read-buffer built by FetchInputData; fallback to rebuild if called standalone.
        if (!read_buf_valid_) RebuildReadBuffer();
        const auto& coords = read_coords_buf_;
        const auto& uids   = read_uid_buf_;

        if (uids.empty()) return;

        std::vector<double> values(uids.size(), 0.0);
        auto coords_span = precice::span<const double>(coords.data(), coords.size());
        auto values_span = precice::span<double>(values.data(), values.size());

        try {
            participant_->mapAndReadData(mesh_name, params_.precice.oxygen_data_name,
                                        coords_span, relative_time, values_span);

            size_t i = 0;
            double min_o2 = std::numeric_limits<double>::max();
            double max_o2 = std::numeric_limits<double>::lowest();
            double sum_o2 = 0.0;
            size_t count_o2 = 0;
            for (const auto& uid : uids) {
                auto* cell = dynamic_cast<CellType*>(rm->GetAgent(uid));
                if (!cell) {
                    ++i;
                    continue;
                }

                const double o2 = values[i];
                if constexpr (has_oxygen_setter_v<CellType>) {
                    cell->SetOxygen(o2);
                } else if constexpr (has_temperature_v<CellType>) {
                    // Backward compatibility: oxygen can be piped through temperature storage.
                    cell->SetTemperature(o2);
                }
                min_o2 = std::min(min_o2, o2);
                max_o2 = std::max(max_o2, o2);
                sum_o2 += o2;
                ++count_o2;
                ++i;
            }

            if (count_o2 > 0 && (coupling_iteration_ < 5 || coupling_iteration_ % 20 == 0)) {
                const double avg_o2 = sum_o2 / static_cast<double>(count_o2);
                bdm::Log::Info("PreciceAdapter", "O2 read stats: count=", count_o2,
                               " min=", min_o2, " max=", max_o2, " avg=", avg_o2);

                // #region agent log
                std::ostringstream payload;
                payload << "{"
                        << "\"iteration\":" << coupling_iteration_ << ","
                        << "\"count\":" << count_o2 << ","
                        << "\"min\":" << min_o2 << ","
                        << "\"max\":" << max_o2 << ","
                        << "\"avg\":" << avg_o2
                        << "}";
                DebugNdjsonLog("baseline", "H5", "bdm_adapter.h:ReadOxygenData",
                               "oxygen_read_stats", payload.str());
                // #endregion
            }
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter", "Failed to read oxygen data: ", e.what());
        }
    }

    void ReadH2O2Data(double relative_time) {
        if (!IsDataNameEnabled(params_.precice.h2o2_data_name)) return;
        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;
        if (!read_buf_valid_) RebuildReadBuffer();
        const auto& coords = read_coords_buf_;
        const auto& uids = read_uid_buf_;
        if (uids.empty()) return;

        std::vector<double> values(uids.size(), 0.0);
        auto coords_span = precice::span<const double>(coords.data(), coords.size());
        auto values_span = precice::span<double>(values.data(), values.size());
        try {
            participant_->mapAndReadData(mesh_name, params_.precice.h2o2_data_name,
                                        coords_span, relative_time, values_span);
            size_t i = 0;
            for (const auto& uid : uids) {
                auto* cell = dynamic_cast<CellType*>(rm->GetAgent(uid));
                if (cell) {
                    if constexpr (has_h2o2_setter_v<CellType>) {
                        cell->SetH2O2(values[i]);
                    }
                }
                ++i;
            }
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter", "Failed to read H2O2 data: ", e.what());
        }
    }

    void ReadNO2Data(double relative_time) {
        if (!IsDataNameEnabled(params_.precice.no2_data_name)) return;
        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;
        if (!read_buf_valid_) RebuildReadBuffer();
        const auto& coords = read_coords_buf_;
        const auto& uids = read_uid_buf_;
        if (uids.empty()) return;

        std::vector<double> values(uids.size(), 0.0);
        auto coords_span = precice::span<const double>(coords.data(), coords.size());
        auto values_span = precice::span<double>(values.data(), values.size());
        try {
            participant_->mapAndReadData(mesh_name, params_.precice.no2_data_name,
                                        coords_span, relative_time, values_span);
            size_t i = 0;
            for (const auto& uid : uids) {
                auto* cell = dynamic_cast<CellType*>(rm->GetAgent(uid));
                if (cell) {
                    if constexpr (has_no2_setter_v<CellType>) {
                        cell->SetNO2(values[i]);
                    }
                }
                ++i;
            }
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter", "Failed to read NO2 data: ", e.what());
        }
    }

    void ReadTKEData(double relative_time) {
        if (!tke_available_) return;

        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;

        // Use shared read-buffer built by FetchInputData; fallback to rebuild if called standalone.
        if (!read_buf_valid_) RebuildReadBuffer();
        const auto& coords = read_coords_buf_;
        const auto& uids   = read_uid_buf_;

        if (uids.empty()) return;

        std::vector<double> values(uids.size(), 0.0);
        auto coords_span = precice::span<const double>(coords.data(), coords.size());
        auto values_span = precice::span<double>(values.data(), values.size());

        try {
            participant_->mapAndReadData(mesh_name, params_.precice.tke_data_name,
                                        coords_span, relative_time, values_span);
            size_t i = 0;
            for (const auto& uid : uids) {
                auto* cell = dynamic_cast<CellType*>(rm->GetAgent(uid));
                if (cell) {
                    if constexpr (has_tke_v<CellType>) {
                        cell->SetTKE(values[i]);
                    }
                }
                ++i;
            }
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter", "Failed to read TKE data: ", e.what());
        }
    }

    void ReadEpsilonData(double relative_time) {
        if (!epsilon_available_) return;

        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;

        // Use shared read-buffer built by FetchInputData; fallback to rebuild if called standalone.
        if (!read_buf_valid_) RebuildReadBuffer();
        const auto& coords = read_coords_buf_;
        const auto& uids   = read_uid_buf_;

        if (uids.empty()) return;

        std::vector<double> values(uids.size(), 0.0);
        auto coords_span = precice::span<const double>(coords.data(), coords.size());
        auto values_span = precice::span<double>(values.data(), values.size());

        try {
            participant_->mapAndReadData(mesh_name, params_.precice.epsilon_data_name,
                                        coords_span, relative_time, values_span);
            size_t i = 0;
            for (const auto& uid : uids) {
                auto* cell = dynamic_cast<CellType*>(rm->GetAgent(uid));
                if (cell) {
                    if constexpr (has_epsilon_v<CellType>) {
                        cell->SetEpsilon(values[i]);
                    }
                }
                ++i;
            }
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter", "Failed to read epsilon data: ", e.what());
        }
    }

    void ReadOmegaData(double relative_time) {
        if (!omega_available_) return;

        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;

        // Use shared read-buffer built by FetchInputData; fallback to rebuild if called standalone.
        if (!read_buf_valid_) RebuildReadBuffer();
        const auto& coords = read_coords_buf_;
        const auto& uids   = read_uid_buf_;

        if (uids.empty()) return;

        std::vector<double> values(uids.size(), 0.0);
        auto coords_span = precice::span<const double>(coords.data(), coords.size());
        auto values_span = precice::span<double>(values.data(), values.size());

        try {
            participant_->mapAndReadData(mesh_name, params_.precice.omega_data_name,
                                        coords_span, relative_time, values_span);
            size_t i = 0;
            for (const auto& uid : uids) {
                auto* cell = dynamic_cast<CellType*>(rm->GetAgent(uid));
                if (cell) {
                    if constexpr (has_omega_v<CellType>) {
                        cell->SetOmega(values[i]);
                    }
                }
                ++i;
            }
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter", "Failed to read omega data: ", e.what());
        }
    }

    void ReadNutData(double relative_time) {
        if (!nut_available_) return;

        const auto& mesh_name = params_.mesh.primary_volume_mesh;
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;

        // Use shared read-buffer built by FetchInputData; fallback to rebuild if called standalone.
        if (!read_buf_valid_) RebuildReadBuffer();
        const auto& coords = read_coords_buf_;
        const auto& uids   = read_uid_buf_;

        if (uids.empty()) return;

        std::vector<double> values(uids.size(), 0.0);
        auto coords_span = precice::span<const double>(coords.data(), coords.size());
        auto values_span = precice::span<double>(values.data(), values.size());

        try {
            participant_->mapAndReadData(mesh_name, params_.precice.nut_data_name,
                                        coords_span, relative_time, values_span);
            size_t i = 0;
            for (const auto& uid : uids) {
                auto* cell = dynamic_cast<CellType*>(rm->GetAgent(uid));
                if (cell) {
                    if constexpr (has_nut_v<CellType>) {
                        cell->SetNut(values[i]);
                    }
                }
                ++i;
            }
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter", "Failed to read nut data: ", e.what());
        }
    }

    void ReadWallMeshData(const std::string& mesh_name, WallMeshData& mesh, double relative_time) {
        if (!mesh.is_initialized || mesh.vertex_ids.empty()) return;

        try {
            // Read wall geometry.
            //
            // For static walls `mesh.geometry` has already been seeded from
            // `mesh.coordinates` in `InitializeWallMeshes`. For moving walls
            // (lid-driven cavity and similar) the CFD side continuously
            // rewrites the WallGeometry data field with the current face
            // centres; we therefore read into a temporary buffer and only
            // adopt it when the values actually look like a populated mesh
            // (not a silently-zeroed buffer). This prevents accidentally
            // clobbering the coordinate-based fallback when the WallGeometry
            // data is not exchanged – observed during restart-from-steady
            // runs on the aneurysm case.
            int geom_dim = participant_->getDataDimensions(mesh_name, params_.precice.wall_geometry_data_name);
            if (geom_dim == 3) {
                const size_t n_entries = mesh.vertex_ids.size() * 3;
                std::vector<double> tmp_geom(n_entries, 0.0);
                auto ids_span = precice::span<const int>(mesh.vertex_ids.data(), mesh.vertex_ids.size());
                auto geom_span = precice::span<double>(tmp_geom.data(), tmp_geom.size());
                participant_->readData(mesh_name, params_.precice.wall_geometry_data_name,
                                      ids_span, relative_time, geom_span);

                bool has_nonzero = false;
                for (double v : tmp_geom) {
                    if (std::abs(v) > 1e-30) { has_nonzero = true; break; }
                }

                if (has_nonzero) {
                    mesh.geometry = std::move(tmp_geom);
                } else if (mesh.geometry.size() != n_entries) {
                    // No usable fallback yet – keep a zero buffer rather than
                    // leaving inconsistent sizes for downstream consumers.
                    mesh.geometry = std::move(tmp_geom);
                }
                // else: retain the coordinates-seeded fallback.
            }

            // Read wall normals. Like WallGeometry, we read into a temporary
            // buffer and only adopt it when it actually carries data. This
            // keeps a previously-read (or seeded) normals buffer valid on
            // subsequent silent-zero reads instead of wiping wall collisions.
            int norm_dim = participant_->getDataDimensions(mesh_name, params_.precice.wall_normals_data_name);
            if (norm_dim == 3) {
                const size_t n_entries = mesh.vertex_ids.size() * 3;
                std::vector<double> tmp_norm(n_entries, 0.0);
                auto ids_span = precice::span<const int>(mesh.vertex_ids.data(), mesh.vertex_ids.size());
                auto norm_span = precice::span<double>(tmp_norm.data(), tmp_norm.size());
                participant_->readData(mesh_name, params_.precice.wall_normals_data_name,
                                      ids_span, relative_time, norm_span);

                bool has_nonzero = false;
                for (double v : tmp_norm) {
                    if (std::abs(v) > 1e-30) { has_nonzero = true; break; }
                }

                if (has_nonzero) {
                    mesh.normals = std::move(tmp_norm);
                }
                // else: retain previously populated normals.
            }

            // Read wall areas. Same defensive strategy: do not overwrite a
            // previously valid buffer with silent zeros.
            int area_dim = participant_->getDataDimensions(mesh_name, params_.precice.wall_areas_data_name);
            if (area_dim == 1) {
                const size_t n_entries = mesh.vertex_ids.size();
                std::vector<double> tmp_area(n_entries, 0.0);
                auto ids_span = precice::span<const int>(mesh.vertex_ids.data(), mesh.vertex_ids.size());
                auto area_span = precice::span<double>(tmp_area.data(), tmp_area.size());
                participant_->readData(mesh_name, params_.precice.wall_areas_data_name,
                                      ids_span, relative_time, area_span);

                bool has_nonzero = false;
                for (double v : tmp_area) {
                    if (std::abs(v) > 1e-30) { has_nonzero = true; break; }
                }

                if (has_nonzero) {
                    mesh.areas = std::move(tmp_area);
                }
                // else: retain previously populated areas.
            }

            // Read wall velocity (for moving walls like lid-driven cavity)
            // This is optional - if not configured in preCICE, walls are stationary
            if (IsDataNameEnabled(params_.precice.wall_velocity_data_name)) {
                try {
                    int vel_dim = participant_->getDataDimensions(mesh_name, params_.precice.wall_velocity_data_name);
                    if (vel_dim == 3) {
                        mesh.velocity.resize(mesh.vertex_ids.size() * 3);
                        auto ids_span = precice::span<const int>(mesh.vertex_ids.data(), mesh.vertex_ids.size());
                        auto vel_span = precice::span<double>(mesh.velocity.data(), mesh.velocity.size());
                        participant_->readData(mesh_name, params_.precice.wall_velocity_data_name,
                                              ids_span, relative_time, vel_span);

                        // Count moving faces for logging
                        size_t moving = 0;
                        for (size_t i = 0; i < mesh.vertex_ids.size(); ++i) {
                            if (std::abs(mesh.velocity[i*3]) > 1e-12 ||
                                std::abs(mesh.velocity[i*3+1]) > 1e-12 ||
                                std::abs(mesh.velocity[i*3+2]) > 1e-12) {
                                ++moving;
                            }
                        }
                        if (moving > 0) {
                            bdm::Log::Info("PreciceAdapter", "Read wall velocity for '", mesh_name,
                                          "': ", moving, "/", mesh.vertex_ids.size(), " moving faces");
                        }
                    }
                } catch (...) {
                    // WallVelocity not configured - walls are stationary (default)
                    mesh.velocity.clear();
                }
            } else {
                mesh.velocity.clear();
            }
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter", "Failed to read wall data for '", mesh_name, "': ", e.what());
        }
    }

    void ReadInletCellCenterData(const std::string& mesh_name, MeshData& mesh, double relative_time) {
        if (!mesh.is_initialized || mesh.vertex_ids.empty()) return;
        if (!IsDataNameEnabled(params_.precice.inlet_cell_centers_data_name)) return;

        try {
            const int data_dim = participant_->getDataDimensions(
                mesh_name, params_.precice.inlet_cell_centers_data_name);
            if (data_dim != 3) {
                bdm::Log::Warning("PreciceAdapter",
                    "Inlet cell-center data '", params_.precice.inlet_cell_centers_data_name,
                    "' on '", mesh_name, "' has dimension ", data_dim,
                    " (expected 3). Falling back to nearest-volume-cell mapping.");
                return;
            }

            std::vector<double> centers(mesh.vertex_ids.size() * 3, 0.0);
            auto ids_span = precice::span<const int>(mesh.vertex_ids.data(), mesh.vertex_ids.size());
            auto centers_span = precice::span<double>(centers.data(), centers.size());
            participant_->readData(mesh_name, params_.precice.inlet_cell_centers_data_name,
                                   ids_span, relative_time, centers_span);

            bool has_finite_nonzero = false;
            for (const double v : centers) {
                if (std::isfinite(v) && std::abs(v) > 1e-30) {
                    has_finite_nonzero = true;
                    break;
                }
            }

            if (!has_finite_nonzero) {
                bdm::Log::Warning("PreciceAdapter",
                    "Inlet cell-center data '", params_.precice.inlet_cell_centers_data_name,
                    "' on '", mesh_name,
                    "' was empty/zero at initialization. Falling back to nearest-volume-cell mapping.");
                return;
            }

            mesh.mapped_cell_centers = std::move(centers);
            bdm::Log::Info("PreciceAdapter",
                "Read exact inlet owner cell centers from '", mesh_name,
                "' using data '", params_.precice.inlet_cell_centers_data_name,
                "' (", mesh.vertex_ids.size(), " entries)");
        } catch (const std::exception& e) {
            bdm::Log::Info("PreciceAdapter",
                "Exact inlet owner cell-center data not available on '", mesh_name,
                "' (", e.what(), "). Falling back to nearest-volume-cell mapping.");
        }
    }

    void ReadInletAreaData(const std::string& mesh_name, MeshData& mesh, double relative_time) {
        if (!mesh.is_initialized || mesh.vertex_ids.empty()) return;
        if (!IsDataNameEnabled(params_.precice.inlet_areas_data_name)) return;

        try {
            const int data_dim = participant_->getDataDimensions(
                mesh_name, params_.precice.inlet_areas_data_name);
            if (data_dim != 1) {
                bdm::Log::Warning("PreciceAdapter",
                    "Inlet area data '", params_.precice.inlet_areas_data_name,
                    "' on '", mesh_name, "' has dimension ", data_dim,
                    " (expected 1). Inlet centering will use its unweighted fallback.");
                return;
            }

            std::vector<double> areas(mesh.vertex_ids.size(), 0.0);
            auto ids_span = precice::span<const int>(mesh.vertex_ids.data(), mesh.vertex_ids.size());
            auto areas_span = precice::span<double>(areas.data(), areas.size());
            participant_->readData(mesh_name, params_.precice.inlet_areas_data_name,
                                   ids_span, relative_time, areas_span);

            double total_area = 0.0;
            for (const double area : areas) {
                if (std::isfinite(area) && area > 0.0) {
                    total_area += area;
                }
            }
            if (!(total_area > 0.0)) {
                bdm::Log::Warning("PreciceAdapter",
                    "Inlet area data '", params_.precice.inlet_areas_data_name,
                    "' on '", mesh_name,
                    "' was empty or non-positive. Inlet centering will use its unweighted fallback.");
                return;
            }

            mesh.face_areas = std::move(areas);
            bdm::Log::Info("PreciceAdapter",
                "Read ", mesh.face_areas.size(), " inlet face areas from '",
                mesh_name, "' (total area=", total_area, " m^2)");
        } catch (const std::exception& e) {
            bdm::Log::Info("PreciceAdapter",
                "Inlet face-area data not available on '", mesh_name,
                "' (", e.what(), "). Inlet centering will use its unweighted fallback.");
        }
    }

    /**
     * @brief Write agent positions to preCICE for two-way coupling
     *
     * Collects positions of all dynamic agents and writes them to ABMmesh
     * for the CFD solver to use (e.g., for particle presence in fluid).
     */
    void WritePositionData() {
        if (!params_.coupling.two_way_coupling) return;
        if (!IsDataNameEnabled(params_.precice.position_data_name)) return;

        const auto& mesh_name = params_.mesh.particle_mesh;  // ABMmesh
        auto* rm = bdm::Simulation::GetActive()->GetResourceManager();
        if (!rm) return;

        // Check if ABMmesh exists and position data is available
        try {
            if (!participant_->requiresMeshConnectivityFor(mesh_name)) {
                // Mesh doesn't require connectivity - we can write positions
            }
        } catch (const std::exception& e) {
            // Mesh might not be configured for writing - skip silently
            return;
        }

        // Collect agent positions
        std::vector<double> positions;
        std::vector<int> vertex_ids;
        positions.reserve(3 * 1024);
        vertex_ids.reserve(1024);

        int vertex_id = 0;
        rm->ForEachAgent([&](bdm::Agent* agent) {
            auto* cell = dynamic_cast<CellType*>(agent);
            if (!cell || cell->IsStationary()) return;

            const auto& p = cell->GetPosition();
            positions.push_back(p[0]);
            positions.push_back(p[1]);
            positions.push_back(p[2]);
            vertex_ids.push_back(vertex_id++);
        });

        if (vertex_ids.empty()) return;

        try {
            // Set mesh vertices (positions define the ABMmesh)
            auto pos_span = precice::span<const double>(positions.data(), positions.size());
            auto ids_span = precice::span<int>(vertex_ids.data(), vertex_ids.size());

            // For two-way coupling, we need to set vertex positions on ABMmesh
            // This requires the mesh to be defined with setMeshVertices
            participant_->setMeshVertices(mesh_name, pos_span, ids_span);

            auto write_pos_span = precice::span<const double>(positions.data(), positions.size());
            participant_->writeData(mesh_name, params_.precice.position_data_name,
                                   precice::span<const int>(vertex_ids.data(), vertex_ids.size()),
                                   write_pos_span);

            bdm::Log::Info("PreciceAdapter", "Wrote ", vertex_ids.size(), " agent positions to '", mesh_name, "'");
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter", "Failed to write position data: ", e.what());
        }
    }

    void WriteScalarSinkData(const std::map<int, double>& sink_source,
                             std::vector<int>& previous_sink_cells,
                             const std::string& sink_data_name,
                             const std::string& log_label) {
        const auto& mesh_name = params_.mesh.primary_volume_mesh;

        try {
            std::vector<int> cell_ids;
            std::vector<double> sink_values;
            std::vector<int> active_cells;
            cell_ids.reserve(sink_source.size() + previous_sink_cells.size());
            sink_values.reserve(sink_source.size() + previous_sink_cells.size());
            active_cells.reserve(sink_source.size());
            double total_sink = 0.0;
            double max_sink = 0.0;

            for (const auto& [cell_id, sink] : sink_source) {
                if (sink <= 0.0) {
                    continue;
                }
                cell_ids.push_back(cell_id);
                sink_values.push_back(sink);
                active_cells.push_back(cell_id);
                total_sink += sink;
                max_sink = std::max(max_sink, sink);
            }

            // Reset cells that were active in the previous step but are inactive now.
            for (const int cell_id : previous_sink_cells) {
                if (sink_source.find(cell_id) != sink_source.end()) {
                    continue;
                }
                cell_ids.push_back(cell_id);
                sink_values.push_back(0.0);
            }

            if (!cell_ids.empty()) {
                auto ids_span = precice::span<const int>(cell_ids.data(), cell_ids.size());
                auto sink_span = precice::span<const double>(sink_values.data(), sink_values.size());
                participant_->writeData(mesh_name, sink_data_name, ids_span, sink_span);
            }

            std::sort(active_cells.begin(), active_cells.end());
            active_cells.erase(std::unique(active_cells.begin(), active_cells.end()), active_cells.end());
            previous_sink_cells.swap(active_cells);

            if (coupling_iteration_ < 5 || coupling_iteration_ % 20 == 0) {
                const size_t cleared_cells = cell_ids.size() >= previous_sink_cells.size()
                    ? (cell_ids.size() - previous_sink_cells.size())
                    : 0;
                bdm::Log::Info("PreciceAdapter", log_label, " sink write: activeCells=",
                               previous_sink_cells.size(), " updatedCells=",
                               cell_ids.size(), " clearedCells=", cleared_cells,
                               " total=", total_sink, " max=", max_sink);

                // #region agent log
                std::ostringstream payload;
                payload << "{"
                        << "\"iteration\":" << coupling_iteration_ << ","
                        << "\"activeCells\":" << previous_sink_cells.size() << ","
                        << "\"updatedCells\":" << cell_ids.size() << ","
                        << "\"clearedCells\":" << cleared_cells << ","
                        << "\"dataName\":\"" << sink_data_name << "\","
                        << "\"totalSink\":" << total_sink << ","
                        << "\"maxSink\":" << max_sink
                        << "}";
                DebugNdjsonLog("baseline", "H2", "bdm_adapter.h:WriteScalarSinkData",
                               "scalar_sink_write_stats", payload.str());
                // #endregion
            }
        } catch (const std::exception& e) {
            bdm::Log::Warning("PreciceAdapter", "Failed to write ", log_label,
                              " sink '", sink_data_name, "': ", e.what());
        }
    }

    /**
     * @brief Write H2O2 sink data to preCICE for two-way scalar coupling
     */
    void WriteOxygenSinkData() {
        if (!params_.coupling.two_way_coupling) return;
        if (!params_.coupling.enable_oxygen_sink) return;
        if (!IsDataNameEnabled(params_.precice.oxygen_sink_data_name)) return;

        WriteScalarSinkData(oxygen_sink_source_, previous_oxygen_sink_cells_,
                            params_.precice.oxygen_sink_data_name, "O2");
    }

    void WriteH2O2SinkData() {
        if (!params_.coupling.two_way_coupling) return;
        if (!params_.coupling.enable_h2o2_sink) return;
        if (!IsDataNameEnabled(params_.precice.h2o2_sink_data_name)) return;

        WriteScalarSinkData(h2o2_sink_source_, previous_h2o2_sink_cells_,
                            params_.precice.h2o2_sink_data_name, "H2O2");
    }

    /**
     * @brief Write NO2- sink data to preCICE for two-way scalar coupling
     */
    void WriteNO2SinkData() {
        if (!params_.coupling.two_way_coupling) return;
        if (!params_.coupling.enable_no2_sink) return;
        if (!IsDataNameEnabled(params_.precice.no2_sink_data_name)) return;

        WriteScalarSinkData(no2_sink_source_, previous_no2_sink_cells_,
                            params_.precice.no2_sink_data_name, "NO2");
    }

    // =========================================================================
    // Member Variables
    // =========================================================================

    AdapterParameters params_;
    Config config_;

    std::unique_ptr<precice::Participant> participant_;

    std::map<std::string, MeshData> volume_meshes_;
    std::map<std::string, WallMeshData> wall_meshes_;
    BoundaryTriangleMeshData boundary_triangle_mesh_;
    std::shared_ptr<const BoundaryGeometry> boundary_geometry_;

    // Inlet/Outlet surface meshes for agent generation and removal
    // These store face centers from OpenFOAM inlet/outlet patches
    std::map<std::string, MeshData> inlet_meshes_;
    std::map<std::string, MeshData> outlet_meshes_;
    std::map<std::string, size_t> inlet_seed_cursor_;

    BoundingBox domain_bounds_;

    double dt_ = 0.0;
    bool is_initialized_ = false;
    int coupling_iteration_ = 0;

    std::size_t fetch_input_call_count_ = 0;
    bool field_data_read_once_ = false;
    bool wall_data_read_once_ = false;
    bool domain_bounds_computed_ = false;
    bool field_freeze_notice_emitted_ = false;
    bool wall_freeze_notice_emitted_ = false;

    // Data availability flags
    bool velocity_available_ = false;
    bool pressure_available_ = false;
    bool pressure_gradient_available_ = false;
    bool fluid_acceleration_available_ = false;
    bool fluid_vorticity_available_ = false;
    bool temperature_available_ = false;
    bool oxygen_available_ = false;
    bool fluid_density_available_ = false;
    bool fluid_dynamic_viscosity_available_ = false;
    bool tke_available_ = false;
    bool epsilon_available_ = false;
    bool omega_available_ = false;
    bool nut_available_ = false;

    // Cached CFD field snapshot on FVMmesh vertices.
    // Used for standalone/frozen-field re-sampling at moving agent positions.
    std::vector<double> cached_velocity_field_;
    std::vector<double> cached_pressure_field_;
    std::vector<double> cached_pressure_gradient_field_;
    std::vector<double> cached_fluid_acceleration_field_;
    std::vector<double> cached_fluid_vorticity_field_;
    std::vector<double> cached_temperature_field_;
    std::vector<double> cached_oxygen_field_;
    std::size_t cached_field_vertex_count_ = 0;
    bool cached_velocity_valid_ = false;
    bool cached_pressure_valid_ = false;
    bool cached_pressure_gradient_valid_ = false;
    bool cached_fluid_acceleration_valid_ = false;
    bool cached_fluid_vorticity_valid_ = false;
    bool cached_temperature_valid_ = false;
    bool cached_oxygen_valid_ = false;
    bool cached_field_notice_emitted_ = false;

    // Two-way scalar sink source terms.
    // O2 sink source [concentration/s per cell]
    std::map<int, double> oxygen_sink_source_;
    // Cell ids that had non-zero oxygen sink in the previous write.
    std::vector<int> previous_oxygen_sink_cells_;
    std::map<int, double> h2o2_sink_source_;
    std::vector<int> previous_h2o2_sink_cells_;
    // NO2- sink source [concentration/s per cell].
    std::map<int, double> no2_sink_source_;
    // Cell ids that had non-zero NO2- sink in the previous write.
    std::vector<int> previous_no2_sink_cells_;

    // Shared read-buffer: agent positions + UIDs, built once per FetchInputData
    // call and reused by all Read*Data functions.  Capacity is preserved across
    // calls (vectors only grow, never shrink) so steady-state == zero allocation.
    std::vector<double>        read_coords_buf_;
    std::vector<bdm::AgentUid> read_uid_buf_;
    bool                       read_buf_valid_ = false;

    // Nearest-FVM-cell cache: AgentUid → preCICE vertex ID.
    // Populated lazily in CollectSourceTerms(); cleared by
    // UpdateDynamicMapping() (only called when motion_enabled=true).
    // Stationary simulations (o2Cube) keep the cache for the entire run,
    // paying the O(N_agents × N_OF_cells) brute-force cost only once per
    // new agent.  New agents from division get a cache-miss on first touch.
    std::unordered_map<bdm::AgentUid, int> agent_fvm_cell_cache_;

    /**
     * @brief Find nearest mesh vertex to given position (particle-to-cell mapping)
     *
     * OpenFOAM equivalent: mesh.findNearestCell(position)
     * This finds which CFD cell contains (or is nearest to) the particle.
     *
     * Nearest-cell assignment:
     *   c(i) = argmin_c ||x_i - x_c||
     * where mesh vertices are the exported OpenFOAM cell centers x_c.
     *
     * @param mesh The volume mesh (FVMmesh with cell centers)
     * @param position Particle position
     * @return Vertex ID (= CFD cell index) or -1 if not found
     */
    template<typename PosType>
    int FindNearestMeshVertexIndex(const MeshData& mesh, const PosType& position) const {
        if (mesh.vertex_ids.empty()) return -1;

        int nearest_index = -1;
        double min_dist_sq = std::numeric_limits<double>::max();

        // Simple brute-force nearest neighbor search.
        // For the current simple case (O(10^2-10^3) CFD cells), this is cheap
        // and robust.
        for (size_t i = 0; i < mesh.NumVertices(); ++i) {
            auto vertex = mesh.GetVertex(i);
            double dx = position[0] - vertex[0];
            double dy = position[1] - vertex[1];
            double dz = position[2] - vertex[2];
            double dist_sq = dx*dx + dy*dy + dz*dz;

            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                nearest_index = static_cast<int>(i);
            }
        }

        return nearest_index;
    }

    template<typename PosType>
    int FindNearestMeshVertex(const MeshData& mesh, const PosType& position) const {
        const int nearest_index = FindNearestMeshVertexIndex(mesh, position);
        if (nearest_index < 0) return -1;
        return mesh.vertex_ids[static_cast<std::size_t>(nearest_index)];
    }
};

} // namespace bdm_adapter

#endif // BDM_ADAPTER_H_
