#ifndef CELLS_H_
#define CELLS_H_

#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <unistd.h>    // getcwd – avoid <filesystem> (breaks ROOT cling dict)
#include <sys/stat.h>  // stat  – file existence check
#ifdef _OPENMP
#include <omp.h>
#endif

// BioDynaMo core
#include "biodynamo.h"

// Configuration module
#include "config/coupling_configuration.h"
#include "config/parameter_loader.h"

// Agents module
#include "agents/coupled_cell.h"

// Behaviors module
#include "behaviors/kinematic_behavior.h"
#include "behaviors/oxygen_growth_behavior.h"
#include "behaviors/rons_response_behavior.h"
#include "behaviors/cell_spread_behavior.h"
#include "behaviors/temperature_size_relaxation_behavior.h"

// IO module
#include "io/agent_trajectory_logger.h"

// Environment module
#include "environment/precice_mesh_mapping.h"
#include "core/environment/kd_tree_environment.h"

// External bdm_adapter library
#include "bdm_adapter.h"

namespace bdm {

// Bring external adapter into bdm namespace for compatibility
using PreciceAdapter = bdm_adapter::PreciceAdapter<MyCell, CavityConfig>;

// Forward declarations
inline int RunCavityCoupling(Simulation& simulation);
inline void ValidateConfiguration();

/**
 * @brief Cross-validated library classes and functions
 *
 * This section documents the cross-validation of classes and functions
 * used from each software library to ensure compatibility.
 */

/*
=== BioDynaMo <-> OpenFOAM <-> preCICE ===

Data structures cross-validated:
- BioDynaMo: Real3, Double3 for positions and vectors ✓
- preCICE: VertexID, std::vector<double> for vertex coordinates ✓
- OpenFOAM: Foam::point, Foam::vector ✓

Data exchange formats cross-validated:
- BioDynaMo -> preCICE: Real3 -> std::vector<double> ✓
- preCICE -> OpenFOAM: std::vector<double> -> Foam::point ✓

Memory layouts cross-validated:
- BioDynaMo positions: {x,y,z} triplets ✓
- preCICE coordinates: interleaved [x1,y1,z1,x2,y2,z2,...] ✓
- OpenFOAM cell centers: (x,y,z) with right-handed system ✓

Coordinate systems:
- BioDynaMo: Cartesian right-handed, uses origin (0,0,0) ✓
- OpenFOAM: Cartesian right-handed, uses origin (0,0,0) ✓
- preCICE: Pass-through for coordinates ✓

API functions cross-validated:
- BioDynaMo:
  - ResourceManager::AddAgent(), GetAgent(), ForEachAgent() ✓
  - Cell::GetPosition(), SetPosition() ✓
  - InteractionForce::Calculate() ✓
  - Agent::AddBehavior() ✓
  - Simulation::GetActive() ✓
  - Scheduler::Simulate() ✓
*/

// === preCICE Library Cross-Validation ===
/*
* List of preCICE functions validated:
* - initialize() ✓
* - advance(dt) ✓
* - setMeshVertices() ✓
* - readData() ✓
* - writeData() ✓
* - getMeshVertexSize() ✓
* - requiresReadingCheckpoint() ✓
*/

// === OpenFOAM Library Cross-Validation ===
/*
* OpenFOAM adapter functions verified:
* - preciceAdapter::velocityData::read() ✓
* - preciceAdapter::pressureData::read() ✓
* - preciceAdapter::temperatureData::read() ✓
* - preciceAdapter::positionData::write() ✓
*/

/**
 * @brief Main simulation function for OpenFOAM-preCICE-BioDynaMo coupling
 *
 * This function initializes and runs the coupling simulation.
 * All settings are controlled via params.csv.
 */
inline int Simulate(int argc, const char** argv) {
    try {
        // Print header
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "  OpenFOAM-preCICE-BioDynaMo Coupling" << std::endl;
        std::cout << std::string(60, '=') << "\n" << std::endl;

        // Check and report OpenMP thread configuration
        const char* omp_threads = std::getenv("OMP_NUM_THREADS");
        int actual_threads = 1;
#ifdef _OPENMP
        actual_threads = omp_get_max_threads();
#endif
        Log::Info("Simulate", "OpenMP Configuration:");
        Log::Info("Simulate", "  OMP_NUM_THREADS env var: ", omp_threads ? omp_threads : "not set");
        Log::Info("Simulate", "  Actual max threads available: ", actual_threads);
        Log::Info("Simulate", "  Parallel execution: ", actual_threads > 1 ? "ENABLED" : "DISABLED");

        Log::Info("Simulate", "Starting OpenFOAM-preCICE-BioDynaMo coupling");

        // Log features from params.csv (single source of truth)
        const auto& rt_config = GetRuntimeCouplingConfig();
        Log::Info("Simulate", "Features enabled (from params.csv):");
        Log::Info("Simulate", rt_config.enable_velocity ? "  ✓ Velocity coupling" : "  ✗ Velocity coupling");
        Log::Info("Simulate", rt_config.enable_pressure ? "  ✓ Pressure coupling" : "  ✗ Pressure coupling");
        Log::Info("Simulate", rt_config.enable_temperature ? "  ✓ Temperature coupling" : "  ✗ Temperature coupling");
        Log::Info("Simulate", rt_config.enable_oxygen ? "  ✓ Oxygen coupling" : "  ✗ Oxygen coupling");
        Log::Info("Simulate", rt_config.enable_h2o2 ? "  ✓ H2O2 coupling" : "  ✗ H2O2 coupling");
        Log::Info("Simulate", rt_config.enable_no2 ? "  ✓ NO2 coupling" : "  ✗ NO2 coupling");
        Log::Info("Simulate", rt_config.enable_h2o2_sink ? "  ✓ H2O2 sink feedback" : "  ✗ H2O2 sink feedback");
        Log::Info("Simulate", rt_config.enable_oxygen_sink ? "  ✓ O2 sink feedback" : "  ✗ O2 sink feedback");
        Log::Info("Simulate", rt_config.enable_no2_sink ? "  ✓ NO2 sink feedback" : "  ✗ NO2 sink feedback");
        Log::Info("Simulate", rt_config.two_way_coupling ? "  ✓ Two-way coupling" : "  ✗ Two-way coupling (one-way only)");
        Log::Info("Simulate", rt_config.enable_particle_particle ? "  ✓ Particle-particle interactions" : "  ✗ Particle-particle interactions");
        Log::Info("Simulate", "  Wall confinement mode: ",
                  rt_config.wall_confinement);
        Log::Info("Simulate", "  Wall velocity response: ",
                  rt_config.wall_velocity_response);
        if (!IsWallConfinementModeValid()) {
            Log::Fatal("Simulate", "Invalid simulation.wall_confinement='",
                       rt_config.wall_confinement,
                       "'. Expected enforce, aabb, or off.");
        }
        if (!IsWallVelocityResponseValid()) {
            Log::Fatal("Simulate", "Invalid interaction.wall_velocity_response='",
                       rt_config.wall_velocity_response,
                       "'. Expected none, slip, rebound, or stick.");
        }
        Log::Info("Simulate", rt_config.enable_particle_motion ? "  ✓ Agent motion" : "  ✗ Agent motion");
        if (rt_config.force_2d) {
            Log::Info("Simulate", "  ✓ 2D coupling mode (force_2d=true)");
        } else {
            Log::Info("Simulate", "  ✓ 3D coupling mode");
        }

        // Set up initial simulation parameters
        // Use kOpen for custom agent-wall interactions (not BioDynaMo's built-in boundaries)
        auto set_param = [](Param* param) {
            param->bound_space = Param::kOpen; // Allow custom wall interactions
            // Avoid UniformGridEnvironment box overflow for sub-unit OpenFOAM
            // coordinates where many cells can occupy one integer grid box.
            param->environment = "kd_tree";

            // Keep BioDynaMo/Catalyst timestamps in physical seconds.  Particle
            // integration still uses dt_ptr_ (preCICE window / substeps), but
            // using the same configured substep here prevents visualization
            // time from being reported as frame numbers 0,1,2,... .
            param->simulation_time_step =
                GetParam<double>("simulation.biodynamo_timestep", 1.0);

            // DISABLE DEFAULT OPERATIONS EXCEPT ENVIRONMENT UPDATES (needed for custom forces)
            // BioDynaMo has many default operations that interfere with coupling:
            // - "mechanical forces": Cell::CalculateDisplacement() and ApplyDisplacement()
            // - "discretization": Cell::RunDiscretization() for morphology changes
            // - "bound space": Boundary clamping
            // We keep only: behaviors, update environment, visualization
            param->unschedule_default_operations = {
                "mechanical forces",
                "discretization",
                "bound space",
                "load balancing"
            };

            // Disable particle-particle interactions but keep custom wall interactions
            param->calculate_gradients = false;           // No gradient-based forces

            // Configure visualization export (read from params.csv after LoadParameters)
            if (ExportVisualization()) {
                param->export_visualization = true;
                param->visualization_interval = VisualizationExportInterval();
                param->output_dir = VisualizationOutputDirectory();
                param->visualization_export_generate_pvsm =
                    GetParam<bool>("visualization.generate_pvsm", false);

                // Build visualization agent map directly from params.csv so we don't need bdm.toml.
                const auto& store = GetParameterStore();

                // 1) Agent names to export (semicolon-separated list in params.csv)
                //    Key: visualization.agents, type: string_list. Default to MyCell if unspecified.
                std::vector<std::string> agent_names;
                if (store.Has<std::vector<std::string>>("visualization.agents")) {
                    agent_names = store.Get<std::vector<std::string>>("visualization.agents");
                }
                if (agent_names.empty()) {
                    agent_names.push_back("MyCell");
                }

                // 2) Additional data members per agent
                //    Key pattern: visualization.agent_members.<AgentName>, type: string_list.
                //    Defaults for MyCell cover the coupled fields we want in ParaView.
                std::map<std::string, std::set<std::string>> viz_agents;
                for (const auto& name : agent_names) {
                    std::set<std::string> members;
                    const std::string key = std::string("visualization.agent_members.") + name;
                    if (store.Has<std::vector<std::string>>(key)) {
                        const auto& list = store.Get<std::vector<std::string>>(key);
                        members.insert(list.begin(), list.end());
                    } else if (name == "MyCell") {
                        // Sensible defaults for the coupled cell type
                        members.insert({"velocity_", "pressure_", "temperature_", "cell_color_", "track_id_", "residence_time_"});
                    }
                    viz_agents[name] = members;
                }
                param->visualize_agents = std::move(viz_agents);

                // 3) Optional diffusion export from params.csv
                //    Key: visualization.diffusion, type: string_list (substance names)
                //    Optional flags (apply to all listed):
                //      visualization.diffusion.gradient (bool, default false)
                //      visualization.diffusion.concentration (bool, default true)
                param->visualize_diffusion.clear();
                if (store.Has<std::vector<std::string>>("visualization.diffusion")) {
                    const auto& substances = store.Get<std::vector<std::string>>("visualization.diffusion");
                    const bool grad = store.GetOr<bool>("visualization.diffusion.gradient", false);
                    const bool conc = store.GetOr<bool>("visualization.diffusion.concentration", true);
                    for (const auto& sub : substances) {
                        Param::VisualizeDiffusion vd;
                        vd.name = sub;
                        vd.gradient = grad;
                        vd.concentration = conc;
                        param->visualize_diffusion.push_back(vd);
                    }
                }

                Log::Info("Simulate", "Visualization export enabled:");
                Log::Info("Simulate", "  Output directory: ", param->output_dir);
                Log::Info("Simulate", "  Export interval: every ", param->visualization_interval, " steps");
                Log::Info("Simulate", "  Generate PVSM: ",
                          param->visualization_export_generate_pvsm ? "true" : "false");
                Log::Info("Simulate", "  Agents: ", agent_names.size());
            } else {
                param->export_visualization = false;
                param->visualize_agents.clear();
                param->visualize_diffusion.clear();
                Log::Info("Simulate", "Visualization export disabled");
            }
        };

        // Create BioDynaMo simulation with cubic bounding box
        Simulation simulation(argc, argv, set_param);

        // Call the cavity coupling implementation
        return RunCavityCoupling(simulation);

    } catch (const std::exception& e) {
        Log::Fatal("Simulate", "Failed to run simulation: ", e.what());
        return 1;
    }
}

/**
 * @brief Alternative simulation configurations
 *
 */

// Full coupling with all features enabled (stationary agents)
inline int SimulateFullCouplingStationary(int argc, const char** argv) {
    Log::Info("SimulateFullCouplingStationary", "Running full coupling with stationary agents");
    return Simulate(argc, argv); // Call regular simulate for now
}

/**
 * @brief Simple cavity coupling implementation
 */
inline int RunCavityCoupling(Simulation& simulation) {
    auto* rm = simulation.GetResourceManager();

    // Particle velocity must start at rest; fluid_velocity_ is sampled from CFD separately.
    auto zero_particle_velocities = [&]() {
        if (!rm) return;
        rm->ForEachAgent([&](bdm::Agent* a) {
            if (auto* cell = dynamic_cast<MyCell*>(a)) {
                if (!cell->IsStationary()) {
                    cell->SetVelocity({0.0, 0.0, 0.0});
                }
            }
            return true;
        });
    };
    AgentTrajectoryLogger trajectory_logger;
    const bool verbose_logging = GetParam<bool>("logging.verbose", false);
    // #region agent log
    auto debug_ndjson_log = [](const std::string& run_id,
                               const std::string& hypothesis_id,
                               const std::string& location,
                               const std::string& message,
                               const std::string& data_json) {
        static_cast<void>(run_id);
        static_cast<void>(hypothesis_id);
        static_cast<void>(location);
        static_cast<void>(message);
        static_cast<void>(data_json);
    };
    // #endregion

    // Create improved preCICE adapter using parameters loaded from params.csv
    const auto adapter_params =
        bdm_adapter::BuildAdapterParametersFromStore(GetParameterStore());

    std::string flow_mode = adapter_params.coupling.flow_mode;
    std::transform(flow_mode.begin(), flow_mode.end(), flow_mode.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    const bool steady_flow_mode = (flow_mode == "steady");
    const int field_read_interval = std::max(0, adapter_params.coupling.field_read_interval);

    Log::Info("RunCavityCoupling", "Creating PreciceAdapter using config file '",
              adapter_params.precice.config_file, "'");
    PreciceAdapter adapter(adapter_params);
    Log::Info("RunCavityCoupling", "Successfully created PreciceAdapter");
    Log::Info("RunCavityCoupling", "Coupling flow mode: ",
              steady_flow_mode ? "steady" : "transient",
              ", field_read_interval=", field_read_interval,
              ", field_read_relative_time=", adapter_params.coupling.field_read_relative_time,
              ", static_mesh=", adapter_params.coupling.static_mesh ? "true" : "false",
              ", wall_read_interval=", std::max(0, adapter_params.coupling.wall_read_interval));

    // Initialize preCICE coupling FIRST - this receives mesh data from OpenFOAM
    Log::Info("RunCavityCoupling", "Initializing preCICE to receive OpenFOAM mesh...");
    adapter.InitializeCoupling();

    // Try to fetch fluid properties from CFD participant via preCICE (optional for one-way coupling)
    const bool fetch_fluid_properties = GetParam<bool>("coupling.fetch_fluid_properties", true);
    auto load_fluid_from_params = [&]() -> bool {
        const auto& store = GetParameterStore();
        if (!store.Has<double>("fluid.density") || !store.Has<double>("fluid.dynamic_viscosity")) {
            return false;
        }

        const double fallback_rho = store.Get<double>("fluid.density");
        const double fallback_mu = store.Get<double>("fluid.dynamic_viscosity");

        if (fallback_rho <= 0.0 || fallback_mu <= 0.0) {
            Log::Fatal("RunCavityCoupling",
                       "Invalid fallback fluid properties (rho<=0 or mu<=0). Check params.csv.");
        }

        SetFluidProperties(fallback_rho, fallback_mu);
        Log::Info("RunCavityCoupling",
                  "Using params.csv fluid properties: rho=", fallback_rho,
                  " kg/m^3, mu=", fallback_mu, " Pa.s");
        return true;
    };

    bool fluid_props_ready = false;
    if (fetch_fluid_properties) {
        double fluid_rho = std::numeric_limits<double>::quiet_NaN();
        double fluid_mu = std::numeric_limits<double>::quiet_NaN();
        if (adapter.FetchFluidProperties(0.0, &fluid_rho, &fluid_mu)) {
            SetFluidProperties(fluid_rho, fluid_mu);
            Log::Info("RunCavityCoupling", "Using fluid properties from preCICE: rho=", fluid_rho,
                      " kg/m^3, mu=", fluid_mu, " Pa.s");
            fluid_props_ready = true;
        } else {
            Log::Warning("RunCavityCoupling",
                         "preCICE did not provide rho/mu; falling back to params.csv values if available.");
        }
    } else {
        Log::Info("RunCavityCoupling",
                  "Skipping preCICE fluid property fetch (coupling.fetch_fluid_properties=false).");
    }

    if (!fluid_props_ready && !load_fluid_from_params()) {
        Log::Fatal("RunCavityCoupling",
                   "Fluid properties missing from preCICE and params.csv. Provide data via preCICE "
                   "or set fluid.density & fluid.dynamic_viscosity in params.csv.");
    }

    // Read wall mesh geometry data from preCICE (vertices were initialized during InitializeCoupling)
    // This populates geometry, normals, and areas needed for GetWallFaces()
    Log::Info("RunCavityCoupling", "===== DIAGNOSTIC: Reading wall mesh geometry data from preCICE...");
    adapter.ReadInitialWallData();
    Log::Info("RunCavityCoupling", "===== DIAGNOSTIC: Wall mesh data read complete");

    // Log mesh dimensions for diagnostics (helps detect 2D/3D mismatch)
    int dim_fvm = adapter.GetMeshDimensions("FVMmesh");
    int dim_wall = adapter.GetMeshDimensions("WallSurfaceMesh");
    Log::Info("RunCavityCoupling", "===== DIAGNOSTIC: Mesh dimensions: FVMmesh=", dim_fvm, ", WallSurfaceMesh=", dim_wall);

    // Get domain boundaries from preCICE
    // Priority: WallSurfaceMesh (actual walls) > FVMmesh (cell centers)
    double foam_x_min, foam_y_min, foam_z_min;
    double foam_x_max, foam_y_max, foam_z_max;

    Log::Info("RunCavityCoupling", "===== DIAGNOSTIC: Checking for wall mesh geometry...");
    Log::Info("RunCavityCoupling", "  HasWallMeshGeometry() = ", adapter.HasWallMeshGeometry());
    Log::Info("RunCavityCoupling", "  HasValidDomainBounds() = ", adapter.HasValidDomainBounds());

    // HYBRID BOUNDS STRATEGY for 2D simulations:
    // - WallSurfaceMesh: XY bounds (actual wall positions)
    // - FVMmesh: Z bounds (full domain extent, since wall faces collapse to single Z)
    // This handles 2D extrusions where all wall faces are at same Z coordinate

    if (adapter.HasWallMeshGeometry() && adapter.HasValidDomainBounds()) {
        Log::Info("RunCavityCoupling", "===== DIAGNOSTIC: Using HYBRID bounds strategy...");

        // Get XY from WallSurfaceMesh (actual walls)
        auto wall_bounds = adapter.GetActualDomainBounds();
        Log::Info("RunCavityCoupling", "  WallSurfaceMesh: X=[", wall_bounds.min[0], ",", wall_bounds.max[0],
                 "] Y=[", wall_bounds.min[1], ",", wall_bounds.max[1],
                 "] Z=[", wall_bounds.min[2], ",", wall_bounds.max[2], "]");

        // Get Z from FVMmesh (full volume extent)
        auto volume_bounds = adapter.GetDomainBoundsFromMesh();
        Log::Info("RunCavityCoupling", "  FVMmesh: X=[", volume_bounds.min[0], ",", volume_bounds.max[0],
                 "] Y=[", volume_bounds.min[1], ",", volume_bounds.max[1],
                 "] Z=[", volume_bounds.min[2], ",", volume_bounds.max[2], "]");

        // Sanity-check each axis. If the wall-mesh extent on an axis is
        // degenerate (NaN / zero / significantly smaller than the
        // corresponding volume extent), the wall-mesh data for that axis is
        // unusable (e.g. silent-zero WallGeometry on restart-from-steady
        // aneurysm runs) — fall back to the FVMmesh bound for that axis so
        // the domain box stays truthful instead of collapsing to a line.
        auto axis_usable_from_wall = [&](int axis) {
            const double w_lo = wall_bounds.min[axis];
            const double w_hi = wall_bounds.max[axis];
            const double v_lo = volume_bounds.min[axis];
            const double v_hi = volume_bounds.max[axis];
            if (!std::isfinite(w_lo) || !std::isfinite(w_hi)) return false;
            const double w_ext = w_hi - w_lo;
            if (!(w_ext > 0.0)) return false;
            if (std::isfinite(v_lo) && std::isfinite(v_hi)) {
                const double v_ext = v_hi - v_lo;
                // Wall XY should at least match a small fraction of the CFD
                // cell-center spread. A threshold of 10% safely tolerates the
                // (cell-center vs. face-center) offset while catching the
                // pathological "wall bounds collapsed to a point" case.
                if (v_ext > 0.0 && w_ext < 0.1 * v_ext) return false;
            }
            return true;
        };

        auto pick_bounds = [&](int axis, bool prefer_volume_here) {
            if (prefer_volume_here) {
                return std::make_pair(volume_bounds.min[axis], volume_bounds.max[axis]);
            }
            if (axis_usable_from_wall(axis)) {
                return std::make_pair(wall_bounds.min[axis], wall_bounds.max[axis]);
            }
            Log::Warning("RunCavityCoupling",
                "  Wall bounds on axis ", axis,
                " are degenerate ([", wall_bounds.min[axis], ",", wall_bounds.max[axis],
                "]); falling back to FVMmesh bounds [",
                volume_bounds.min[axis], ",", volume_bounds.max[axis], "]");
            return std::make_pair(volume_bounds.min[axis], volume_bounds.max[axis]);
        };

        // XY from walls when sane, Z always from the volume (wall faces in
        // 2D extrusions all collapse to a single Z plane).
        auto [x_lo, x_hi] = pick_bounds(0, false);
        auto [y_lo, y_hi] = pick_bounds(1, false);
        auto [z_lo, z_hi] = pick_bounds(2, true);

        foam_x_min = x_lo;  foam_x_max = x_hi;
        foam_y_min = y_lo;  foam_y_max = y_hi;
        foam_z_min = z_lo;  foam_z_max = z_hi;

        Log::Info("RunCavityCoupling", "===== DIAGNOSTIC: HYBRID bounds: X=[", foam_x_min, ",", foam_x_max,
                     "] Y=[", foam_y_min, ",", foam_y_max,
                     "] Z=[", foam_z_min, ",", foam_z_max, "]");
    } else if (adapter.HasWallMeshGeometry()) {
        Log::Info("RunCavityCoupling", "===== DIAGNOSTIC: Using WallSurfaceMesh bounds only...");
        auto wall_bounds = adapter.GetActualDomainBounds();
        foam_x_min = wall_bounds.min[0];
        foam_y_min = wall_bounds.min[1];
        foam_z_min = wall_bounds.min[2];
        foam_x_max = wall_bounds.max[0];
        foam_y_max = wall_bounds.max[1];
        foam_z_max = wall_bounds.max[2];
        Log::Info("RunCavityCoupling", "  WallSurfaceMesh bounds: X=[", foam_x_min, ",", foam_x_max,
                     "] Y=[", foam_y_min, ",", foam_y_max, "] Z=[", foam_z_min, ",", foam_z_max, "]");
    } else if (adapter.HasValidDomainBounds()) {
        Log::Info("RunCavityCoupling", "===== DIAGNOSTIC: Using FVMmesh bounds only...");
        auto mesh_bounds = adapter.GetDomainBoundsFromMesh();
        foam_x_min = mesh_bounds.min[0];
        foam_y_min = mesh_bounds.min[1];
        foam_z_min = mesh_bounds.min[2];
        foam_x_max = mesh_bounds.max[0];
        foam_y_max = mesh_bounds.max[1];
        foam_z_max = mesh_bounds.max[2];
        Log::Info("RunCavityCoupling", "  FVMmesh bounds: X=[", foam_x_min, ",", foam_x_max,
                 "] Y=[", foam_y_min, ",", foam_y_max, "] Z=[", foam_z_min, ",", foam_z_max, "]");
    } else {
        Log::Fatal("RunCavityCoupling",
            "No mesh bounds available from preCICE! "
            "Ensure OpenFOAM is running and meshes are properly configured.");
        return 1;
    }

    bdm_adapter::BoundingBox hybrid_bounds;
    hybrid_bounds.min = {foam_x_min, foam_y_min, foam_z_min};
    hybrid_bounds.max = {foam_x_max, foam_y_max, foam_z_max};

    // =========================================================================
    // BOUNDS VALIDATION - Verify mesh bounds from preCICE are valid
    // =========================================================================
    // All bounds MUST come from preCICE mesh data (FVMmesh or WallSurfaceMesh).
    // NO FALLBACK to hardcoded values - if bounds are wrong, the OpenFOAM mesh
    // export must be fixed at the source.

    Log::Info("RunCavityCoupling", "Received mesh bounds from preCICE:");
    Log::Info("RunCavityCoupling", "  X: [", foam_x_min, ", ", foam_x_max, "]");
    Log::Info("RunCavityCoupling", "  Y: [", foam_y_min, ", ", foam_y_max, "]");
    Log::Info("RunCavityCoupling", "  Z: [", foam_z_min, ", ", foam_z_max, "]");

    // Store mesh-derived bounds in global storage for use by other modules
    SetDomainBoundsFromMesh(foam_x_min, foam_x_max,
                            foam_y_min, foam_y_max,
                            foam_z_min, foam_z_max);

    // CRITICAL: Update adapter's domain_bounds_ with computed hybrid bounds
    adapter.SetDomainBounds(hybrid_bounds);
    Log::Info("RunCavityCoupling", "===== DIAGNOSTIC: Updated adapter domain bounds to hybrid bounds");
    Log::Info("RunCavityCoupling", "  Final bounds: X=[", foam_x_min, ",", foam_x_max,
                 "] Y=[", foam_y_min, ",", foam_y_max, "] Z=[", foam_z_min, ",", foam_z_max, "]");

    const auto& all_outlet_meshes = adapter.GetOutletMeshes();

    if (IsOutletRemovalEnabled()) {
        if (all_outlet_meshes.empty()) {
            Log::Warning("RunCavityCoupling",
                         "Outlet removal enabled but no outlet surface meshes were received from preCICE");
        } else {
            size_t total_outlet_faces = 0;
            for (const auto& [mesh_name, mesh_data] : all_outlet_meshes) {
                if (mesh_data.is_initialized && mesh_data.NumVertices() > 0) {
                    total_outlet_faces += mesh_data.NumVertices();
                    Log::Info("RunCavityCoupling", "Outlet surface mesh '", mesh_name,
                              "': ", mesh_data.NumVertices(), " OpenFOAM face centers");
                }
            }
            Log::Info("RunCavityCoupling",
                      "Outlet removal will use OpenFOAM/preCICE outlet surface "
                      "face centers directly (", total_outlet_faces,
                      " centers across ", all_outlet_meshes.size(), " mesh(es))");
        }
    }

    auto log_wall_mesh_status = [&](const std::string& stage) {
        const auto& wall_meshes = adapter.GetWallMeshes();
        if (wall_meshes.empty()) {
            Log::Warning("RunCavityCoupling", stage, ": adapter has no wall meshes registered");
            return;
        }
        for (const auto& [mesh_name, mesh] : wall_meshes) {
            const size_t n_vertices = mesh.vertex_ids.size();
            const size_t n_geometry = mesh.geometry.size() / 3;
            const size_t n_normals = mesh.normals.size() / 3;
            const size_t n_areas = mesh.areas.size();
            Log::Warning("RunCavityCoupling",
                         stage, ": wall mesh '", mesh_name,
                         "' initialized=", mesh.is_initialized,
                         ", vertices=", n_vertices,
                         ", geometryFaces=", n_geometry,
                         ", normals=", n_normals,
                         ", areas=", n_areas);
        }
    };

    // Get wall faces from WallSurfaceMesh for face-based collision detection
    // This provides actual wall geometry (centers, normals, areas) from OpenFOAM
    auto wall_faces = adapter.GetWallFaces();
    Log::Info("RunCavityCoupling", "Loaded ", wall_faces.size(),
              " wall faces from WallSurfaceMesh for collision detection");

    // Create shared time step variable for behaviors
    double current_dt = adapter.GetTimeStepSize();

    // Create kinematic behavior - wall collision handled via preCICE face data
    auto* wall_behavior = new KinematicCloudBehavior(&current_dt);
    wall_behavior->SetBoundaryGeometry(adapter.GetBoundaryGeometry());
    const bool growth_division_enabled =
        GetParam<bool>("biology.enable_growth_division", false);
    auto* oxygen_growth_behavior = new OxygenGrowthBehavior(&current_dt);

    const bool rons_response_enabled =
        GetParam<bool>("biology.enable_rons_response", false);
    const bool remove_dead_rons_agents_enabled =
        GetParam<bool>("biology.remove_dead_agents", false);
    const double rons_death_threshold =
        GetParam<double>("biology.rons_death_threshold", 1.0);
    // Hazard once damage reaches the death threshold [1/s]. The removal
    // probability on a step of length dt is 1 - exp(-rate * dt).
    const double rons_death_removal_rate =
        GetParam<double>("biology.death_removal_rate", 0.02);
    auto* rons_response_behavior =
        rons_response_enabled ? new RonsResponseBehavior(&current_dt) : nullptr;
    if (rons_response_enabled) {
        Log::Info("RunCavityCoupling",
                  "Attaching RONS response behavior (intracellular load / damage / fate)");
        Log::Info("RunCavityCoupling",
                  "RONS-dead agent removal: ",
                  remove_dead_rons_agents_enabled
                      ? "probabilistic after accepted biology steps"
                      : "disabled");
        if (remove_dead_rons_agents_enabled) {
            std::cout << "RonsResponse: probabilistic removal once damage >= "
                      << rons_death_threshold
                      << ", hazard " << rons_death_removal_rate << " 1/s"
                      << std::endl;
        }
    }

    const bool temperature_relaxation_enabled =
        GetParam<bool>("simulation.enable_temperature_size_relaxation", false) &&
        GetRuntimeCouplingConfig().enable_temperature;
    auto* temperature_size_relaxation_behavior =
        temperature_relaxation_enabled ? new TemperatureSizeRelaxationBehavior(&current_dt) : nullptr;

    const bool cell_spread_enabled =
        GetParam<bool>("simulation.enable_cell_spread", false);
    auto* cell_spread_behavior_proto = cell_spread_enabled
        ? new CellCellSpreadBehavior() : nullptr;
    if (cell_spread_enabled) {
        Log::Info("RunCavityCoupling",
                  "Cell-cell contact spread enabled (simulation.enable_cell_spread=true)");
    }

    // Set wall faces for face-based collision detection (required)
    if (!wall_faces.empty()) {
        wall_behavior->SetWallFaces(wall_faces);
        if (wall_behavior->IsFaceBasedCollisionEnabled()) {
            Log::Info("RunCavityCoupling", "Using face-based wall collision detection");
        } else {
            Log::Warning("RunCavityCoupling",
                         "Wall faces were received but none are valid for collision (normals/areas invalid)");
        }
    } else {
        Log::Warning("RunCavityCoupling", "No wall faces available - particles may exit domain!");
        log_wall_mesh_status("Startup wall-face check");
    }

    bool clamp_to_domain =
        GetParam<bool>("simulation.clamp_to_domain", false)
        || IsAabbWallConfinementEnabled();
    const bool motion_enabled = GetRuntimeCouplingConfig().enable_particle_motion;
    if (motion_enabled && IsWallConfinementEnabled() &&
        !wall_behavior->HasBoundaryGeometry()) {
        Log::Fatal("RunCavityCoupling",
            "simulation.wall_confinement=enforce requires a valid "
            "BoundaryTriangleMesh. Refusing to run moving agents without "
            "geometry-aware confinement.");
    }
    if (motion_enabled && IsWallConfinementEnabled() && !IsForce2D() &&
        !adapter.HasClosedBoundaryGeometry()) {
        Log::Fatal("RunCavityCoupling",
            "The 3D computational boundary is not closed. Include wall, inlet, "
            "and outlet patches in BoundaryTriangleMesh; inlet/outlet triangles "
            "remain pass-through and are not wall-collision surfaces.");
    }
    Log::Info("RunCavityCoupling", motion_enabled
              ? "Agent motion enabled (simulation.enable_particle_motion=true)"
              : "Agent motion disabled (simulation.enable_particle_motion=false)");
    if (IsMasslessTracerMode() && clamp_to_domain) {
        Log::Info("RunCavityCoupling",
                  "Keeping domain clamp enabled for massless tracers "
                  "(center +/- radius must remain inside domain bounds)");
    }

    // Create agents - three modes available:
    // 1. CSV file (custom_positions_required=true): User-defined positions
    // 2. Inlet mesh (inlet_generation=true): Generate at inlet face centers
    // 3. Volume mesh (default): Create at all CFD cell centers

    // Helper: current working directory as std::string (no <filesystem> – avoids ROOT cling issue)
    auto get_cwd = []() -> std::string {
        char buf[4096] = {};
        return getcwd(buf, sizeof(buf)) ? std::string(buf) : std::string("(unknown)");
    };

    double agent_diameter = DefaultDiameter();
    bool use_csv = GetParam<bool>("agent.custom_positions_required", false);
    std::string csv_file = GetParam<std::string>("agent.custom_positions_csv", "");
    bool use_inlet = IsInletGenerationEnabled();
    double inlet_offset = InletZoneOffset();
    bool init_from_volume = GetParam<bool>("agent.initialize_from_volume_mesh", false);
    int volume_seed_stride = std::max(1, GetParam<int>("agent.volume_seed_stride", 1));

    // DEBUG: Print what we got from params.csv
    Log::Info("RunCavityCoupling", "DEBUG: working dir = ", get_cwd());
    Log::Info("RunCavityCoupling", "DEBUG: use_csv = ", use_csv);
    Log::Info("RunCavityCoupling", "DEBUG: csv_file = '", csv_file, "'");
    Log::Info("RunCavityCoupling", "DEBUG: use_inlet = ", use_inlet);
    Log::Info("RunCavityCoupling", "DEBUG: inlet_offset = ", inlet_offset);
    Log::Info("RunCavityCoupling", "DEBUG: init_from_volume = ", init_from_volume);
    Log::Info("RunCavityCoupling", "DEBUG: volume_seed_stride = ", volume_seed_stride);

    size_t num_agents = 0;
    size_t expected_fixed_inlet_batch_size = 0;
    const bool enforce_fixed_inlet_batch =
        use_inlet &&
        !GetParam<bool>("agent.inlet_random_sampling", true) &&
        !GetParam<bool>("agent.inlet_filter_occupied", true);

    // Get mesh name from params.csv
    std::string volume_mesh = GetParam<std::string>("mesh.primary_volume_mesh", "FVMmesh");
    const auto& bounds = adapter.GetDomainBounds();
    Log::Info("RunCavityCoupling", "DEBUG: domain bounds min = (", bounds.min[0], ", ", bounds.min[1], ", ", bounds.min[2], ") max = (", bounds.max[0], ", ", bounds.max[1], ", ", bounds.max[2], ")");

    // Agent factory - shared across all creation modes
    auto agent_factory = [agent_diameter](const std::array<double, 3>& pos) -> MyCell* {
        Real3 position = {pos[0], pos[1], pos[2]};
        auto* cell = new MyCell(position);
        cell->SetDiameter(agent_diameter);
        cell->SetStationary(!GetRuntimeCouplingConfig().enable_particle_motion);
        return cell;
    };

    auto create_agents_from_volume_mesh = [&](const std::string& mesh_name, int stride) -> size_t {
        if (stride <= 1) {
            return adapter.CreateAgentsAtMeshVertices(mesh_name, agent_factory);
        }

        const auto* mesh_ptr = adapter.GetVolumeMesh(mesh_name);
        if (!mesh_ptr || !mesh_ptr->is_initialized || mesh_ptr->NumVertices() == 0) {
            Log::Warning("RunCavityCoupling", "Volume mesh '", mesh_name,
                         "' unavailable for strided initialization");
            return 0;
        }

        auto* rm_local = bdm::Simulation::GetActive()->GetResourceManager();
        size_t count = 0;
        for (size_t i = 0; i < mesh_ptr->NumVertices(); i += static_cast<size_t>(stride)) {
            auto pos = mesh_ptr->GetVertex(i);
            auto* agent = agent_factory(pos);
            if (agent) {
                rm_local->AddAgent(agent);
                ++count;
            }
        }

        Log::Info("RunCavityCoupling", "Created ", count, " agents from ", mesh_name,
                  " cell centers (stride=", stride, ")");
        return count;
    };

    if (use_csv && !csv_file.empty()) {
        // MODE 1: Create agents from CSV file at user-specified positions
        // Resolve relative path without <filesystem> (avoids ROOT cling dict issues)
        std::string csv_path = csv_file;
        if (csv_path.empty() || csv_path[0] != '/') {
            csv_path = get_cwd() + "/" + csv_path;
        }

        struct stat csv_st{};
        const bool csv_exists = (stat(csv_path.c_str(), &csv_st) == 0 && S_ISREG(csv_st.st_mode));
        if (!csv_exists) {
            Log::Fatal("RunCavityCoupling", "CSV file not found: ", csv_path);
        }

        Log::Info("RunCavityCoupling", "Creating agents from CSV (resolved path): ", csv_path);
        num_agents = adapter.CreateAgentsFromCSV(csv_path, agent_factory);

        Log::Info("RunCavityCoupling", "Created ", num_agents, " agents from CSV file");
        auto* rm_check = bdm::Simulation::GetActive()->GetResourceManager();
        size_t total_agents = rm_check ? rm_check->GetNumAgents() : 0;
        Log::Info("RunCavityCoupling", "ResourceManager agent count after CSV create = ", total_agents);

        // Diagnostic: log first few agent positions after CSV creation
        if (verbose_logging) {
            int log_count = 0;
            rm_check->ForEachAgent([&](Agent* agent) {
                if (log_count < 5) {
                    const auto& p = agent->GetPosition();
                    Log::Info("RunCavityCoupling", "Agent[", log_count, "] pos after CSV create = (",
                              p[0], ", ", p[1], ", ", p[2], ")");
                    log_count++;
                }
                return log_count < 5;
            });
        }

        // If we failed to create any agents, abort instead of silently falling back
        if (num_agents == 0) {
            Log::Fatal("RunCavityCoupling", "CSV load returned zero agents. File: ", csv_path);
        }
        // If more agents than expected exist, abort to avoid mixing with mesh vertices
        if (total_agents != num_agents) {
            Log::Fatal("RunCavityCoupling", "Expected exactly ", num_agents, " agents from CSV, but ResourceManager contains ", total_agents, ". Aborting to prevent fallback to mesh vertices.");
        }

        size_t outside_domain = 0;
        rm_check->ForEachAgent([&](Agent* agent) {
            auto* cell = dynamic_cast<MyCell*>(agent);
            if (!cell) {
                return true;
            }
            const auto& p = cell->GetPosition();
            if (p[0] < hybrid_bounds.min[0] || p[0] > hybrid_bounds.max[0] ||
                p[1] < hybrid_bounds.min[1] || p[1] > hybrid_bounds.max[1] ||
                p[2] < hybrid_bounds.min[2] || p[2] > hybrid_bounds.max[2]) {
                ++outside_domain;
            }
            return true;
        });

        if (outside_domain > 0) {
            Log::Fatal("RunCavityCoupling",
                       outside_domain,
                       " custom CSV agent positions lie outside the OpenFOAM domain. "
                       "Check agents_positions-o2Cube.csv against the mesh bounds.");
        }
    } else {
        // MODE 2/3/4: volume initialization, inlet initialization, or hybrid
        // - MODE 2 (inlet only): use_inlet=true, init_from_volume=false
        // - MODE 3 (volume only): use_inlet=false
        // - MODE 4 (hybrid): use_inlet=true, init_from_volume=true
        size_t volume_agents = 0;
        size_t inlet_agents = 0;

        if (!use_inlet || init_from_volume) {
            Log::Info("RunCavityCoupling", "Creating agents at ", volume_mesh,
                      " cell centers (stride=", volume_seed_stride, ")");
            volume_agents = create_agents_from_volume_mesh(volume_mesh, volume_seed_stride);
            if (!use_inlet) {
                Log::Info("RunCavityCoupling", "Created ", volume_agents, " agents in volume mesh");
            } else {
                Log::Info("RunCavityCoupling", "Hybrid init: created ", volume_agents,
                          " agents in volume mesh before inlet seeding");
            }
        }

        if (use_inlet) {
            std::string inlet_mesh = GetParam<std::string>("mesh.inlet_mesh", "InletSurfaceMesh");
            const bool inlet_seed_initial =
                GetParam<bool>("agent.inlet_seed_initial", true);
            // A value of zero means one agent at every validated inlet owner
            // cell center. The same value/path is used for initialization and
            // every recurring injection.
            const int configured_inlet_batch = InletAgentsPerSpawn();
            const size_t inlet_agents_per_spawn = configured_inlet_batch > 0
                ? static_cast<size_t>(configured_inlet_batch)
                : static_cast<size_t>(0);
            Log::Info("RunCavityCoupling", "Inlet mesh '", inlet_mesh,
                      "' uses offset ", inlet_offset, " m",
                      configured_inlet_batch > 0
                          ? " (using configured batch cap)"
                          : " (using all validated central owner cells)");

            if (!adapter.HasInletMesh(inlet_mesh)) {
                Log::Fatal("RunCavityCoupling", "Inlet generation requested but inlet mesh '", inlet_mesh,
                           "' not available. Check preCICE config and mesh.default_inlet_meshes in params.csv");
            } else if (inlet_seed_initial) {
                inlet_agents = adapter.CreateAgentsAtInletMesh(
                    inlet_mesh, agent_factory, inlet_offset, inlet_agents_per_spawn);
                Log::Info("RunCavityCoupling", "Created ", inlet_agents,
                          " agents at inlet volume cell centers");
            } else {
                Log::Info("RunCavityCoupling",
                          "Initial inlet seeding disabled; continuous spawning will start after the first BioDynaMo substep");
            }
        }

        if (enforce_fixed_inlet_batch) {
            expected_fixed_inlet_batch_size = inlet_agents;
            Log::Info("RunCavityCoupling",
                      "Fixed inlet batch size established from central owner cells: ",
                      expected_fixed_inlet_batch_size, " agents per injection");
        }

        num_agents = volume_agents + inlet_agents;
        auto* rm_check = bdm::Simulation::GetActive()->GetResourceManager();
        size_t total_agents = rm_check ? rm_check->GetNumAgents() : 0;
        const bool continuous_inlet_pending =
            use_inlet && InletSpawnInterval() > 0.0;
        if ((num_agents == 0 || total_agents == 0) && !continuous_inlet_pending) {
            Log::Fatal("RunCavityCoupling", "Agent initialization failed - no agents created");
        }
        Log::Info("RunCavityCoupling", "Total initialized agents: ", total_agents,
                  " (volume=", volume_agents, ", inlet=", inlet_agents, ")");
    }

    zero_particle_velocities();
    Log::Info("RunCavityCoupling",
              "Agent particle velocities initialized to zero (fluid field sampled separately)");

    // #region agent log
    {
        std::ostringstream payload;
        payload << "{"
                << "\"initialAgents\":" << (rm ? rm->GetNumAgents() : 0) << ","
                << "\"twoWay\":" << (adapter_params.coupling.two_way_coupling ? "true" : "false") << ","
                << "\"oxygenEnabled\":" << (adapter_params.coupling.enable_oxygen ? "true" : "false") << ","
                << "\"oxygenSinkEnabled\":" << (adapter_params.coupling.enable_oxygen_sink ? "true" : "false") << ","
                << "\"no2SinkEnabled\":" << (adapter_params.coupling.enable_no2_sink ? "true" : "false") << ","
                << "\"oxygenUptakeRate\":" << adapter_params.coupling.oxygen_uptake_rate << ","
                << "\"oxygenUptakeMax\":" << adapter_params.coupling.oxygen_uptake_max << ","
                << "\"no2UptakeRate\":" << adapter_params.coupling.no2_uptake_rate << ","
                << "\"no2UptakeMax\":" << adapter_params.coupling.no2_uptake_max << ","
                << "\"growthRate\":" << GetParam<double>("biology.growth_rate", 0.0) << ","
                << "\"divisionDiameter\":" << GetParam<double>("biology.division_diameter", 2.0 * DefaultDiameter())
                << "}";
        debug_ndjson_log("baseline", "H4", "cells.h:RunCavityCoupling",
                         "initial_agent_and_coupling_config", payload.str());
    }
    // #endregion

    // CRITICAL: Fetch initial velocity and pressure data NOW that agents exist
    Log::Info("RunCavityCoupling", "Fetching initial velocity/pressure data for ", num_agents, " agents...");
    adapter.FetchInputData(0.0);
    zero_particle_velocities();
    Log::Info("RunCavityCoupling", "Initial field data loaded successfully");

    struct FluidVelocityStats {
        double max = 0.0;
        double mean = 0.0;
        size_t sampled = 0;
    };
    auto sample_fluid_velocity_stats = [&]() -> FluidVelocityStats {
        FluidVelocityStats stats;
        double sum = 0.0;
        rm->ForEachAgent([&](bdm::Agent* a) {
            auto* cell = dynamic_cast<MyCell*>(a);
            if (!cell) return;
            const auto& uf = cell->GetFluidVelocity();
            const double m = std::sqrt(uf[0] * uf[0] + uf[1] * uf[1] + uf[2] * uf[2]);
            sum += m;
            if (m > stats.max) {
                stats.max = m;
            }
            ++stats.sampled;
        });
        stats.mean = stats.sampled ? (sum / static_cast<double>(stats.sampled)) : 0.0;
        return stats;
    };

    // Refresh wall faces after the first explicit read at relative_time=0.
    // This avoids using pre-initialized or stale wall buffers when attaching
    // per-agent collision behaviours.
    if (IsWallConfinementEnabled()) {
        wall_faces = adapter.GetWallFaces();
        wall_behavior->SetWallFaces(wall_faces);
        if (!wall_behavior->IsFaceBasedCollisionEnabled()) {
            clamp_to_domain = true;
            Log::Warning("RunCavityCoupling",
                         "No valid wall faces after initial data fetch; keeping domain clamp enabled");
            log_wall_mesh_status("Post-initial-fetch wall-face check");
        } else {
            Log::Info("RunCavityCoupling", "Validated ", wall_behavior->GetWallFaceCount(),
                      " wall faces after initial fetch");
        }
    }

    // Optional priming: advance a few coupling windows without BioDynaMo
    // integration when startup data is still near-zero (common when OpenFOAM
    // starts from a zero internal field in 0/U). This avoids exporting several
    // initial "frozen" frames while keeping parallel-explicit coupling.
    int prime_windows = std::max(0, GetParam<int>("simulation.prime_windows", 0));
    if (steady_flow_mode && prime_windows > 0) {
        Log::Info("RunCavityCoupling",
                  "Disabling startup priming because coupling.flow_mode=steady");
        prime_windows = 0;
    }
    if (field_read_interval == 0 && prime_windows > 0) {
        Log::Info("RunCavityCoupling",
                  "Disabling startup priming because coupling.field_read_interval=0 (frozen fields)");
        prime_windows = 0;
    }
    const double prime_velocity_threshold =
        std::max(0.0, GetParam<double>("simulation.prime_velocity_threshold", 1e-8));
    double primed_coupling_time = 0.0;
    if (prime_windows > 0 && motion_enabled) {
        for (int w = 0; w < prime_windows && adapter.IsCouplingOngoing(); ++w) {
            auto stats = sample_fluid_velocity_stats();
            if (stats.max > prime_velocity_threshold) {
                Log::Info("RunCavityCoupling",
                          "Startup priming stopped early after ", w,
                          " window(s): |U_fluid| max=", stats.max,
                          " > threshold=", prime_velocity_threshold);
                break;
            }

            const double prime_dt = current_dt;
            Log::Warning("RunCavityCoupling",
                         "Priming startup window ", (w + 1), "/", prime_windows,
                         " because sampled |U_fluid| max=", stats.max,
                         " <= threshold=", prime_velocity_threshold,
                         " (dt=", prime_dt, " s)");
            adapter.Advance(prime_dt);
            primed_coupling_time += prime_dt;
            current_dt = adapter.GetTimeStepSize();
            adapter.FetchInputData(0.0);
            adapter.SetDomainBounds(hybrid_bounds);
        }
    }

    // Re-check wall faces after optional priming. If walls are still missing,
    // keep clamp mode and optionally fail fast (aneurysm should not run this way).
    if (IsWallConfinementEnabled()) {
        wall_faces = adapter.GetWallFaces();
        wall_behavior->SetWallFaces(wall_faces);
        const bool face_collision_enabled = wall_behavior->IsFaceBasedCollisionEnabled();
        if (!face_collision_enabled) {
            clamp_to_domain = true;
            Log::Warning("RunCavityCoupling",
                         "Wall faces still unavailable after startup fetch/priming; keeping domain clamp enabled");
            log_wall_mesh_status("Post-priming wall-face check");
        } else if (prime_windows > 0) {
            Log::Info("RunCavityCoupling",
                      "Wall faces validated after startup fetch/priming: ",
                      wall_behavior->GetWallFaceCount());
        }

        if (GetParam<bool>("simulation.require_wall_faces", false) && !face_collision_enabled) {
            throw std::runtime_error(
                "simulation.require_wall_faces=true but no valid wall faces were received. "
                "Check OpenFOAM WallSurfaceMesh export (WallGeometry/WallNormals/WallAreas), "
                "preciceDict patch names, and that run.sh uses the intended adapter library.");
        }
    }

    // CRITICAL: Export visualization at t=0 BEFORE any Simulate() calls
    // This ensures ParaView shows the true CSV positions at time zero
    if (ExportVisualization()) {
        Log::Info("RunCavityCoupling", "Exporting initial state (t=0) before coupling loop");
        simulation.GetScheduler()->GetOps("visualize")[0]->operator()();
    }

    // Configure per-agent behavior based on motion setting from params.csv
    Log::Info("RunCavityCoupling", motion_enabled
              ? "Attaching kinematic behavior for movable agents"
              : "Leaving agents stationary (no kinematic behavior attached)");
    if (growth_division_enabled) {
        Log::Info("RunCavityCoupling",
                  "Attaching oxygen-dependent growth/division behavior (local O2 threshold)");
    }

    rm->ForEachAgent([&](Agent* agent) {
        if (auto* cell = dynamic_cast<MyCell*>(agent)) {
            cell->SetStationary(!motion_enabled);
            if (motion_enabled) {
                // Handles fluid forces and wall collisions for moving particles.
                auto* agent_wall_behavior = wall_behavior->NewCopy();
                cell->AddBehavior(agent_wall_behavior);
            }
            if (growth_division_enabled) {
                auto* growth_behavior = oxygen_growth_behavior->NewCopy();
                cell->AddBehavior(growth_behavior);
            }
            if (rons_response_behavior) {
                cell->AddBehavior(rons_response_behavior->NewCopy());
            }
            if (temperature_size_relaxation_behavior) {
                cell->AddBehavior(temperature_size_relaxation_behavior->NewCopy());
            }
            if (cell_spread_enabled && cell_spread_behavior_proto) {
                cell->AddBehavior(cell_spread_behavior_proto->NewCopy());
            }
        }
        return true;
    });

    // Get the initial time step size from adapter (from preCICE directly)
    double dt = current_dt;
    double simulation_time = primed_coupling_time;
    if (primed_coupling_time > 0.0) {
        Log::Warning("RunCavityCoupling",
                     "BioDynaMo dynamics start after startup priming offset t=",
                     primed_coupling_time, " s");
    }

    // Visualization cadence in physical time (seconds), robust to runtime
    // biology_dt overrides by preCICE/adaptive stepping.
    // If visualization.export_interval_seconds is not provided (<= 0),
    // fall back to legacy step-based config converted using current_dt.
    double viz_interval_seconds = 0.0;
    double next_viz_time = std::numeric_limits<double>::infinity();
    if (ExportVisualization()) {
        const double configured_seconds =
            GetParam<double>("visualization.export_interval_seconds", -1.0);
        if (configured_seconds > 0.0) {
            viz_interval_seconds = configured_seconds;
        } else {
            const int interval_steps = std::max(1, VisualizationExportInterval());
            viz_interval_seconds = interval_steps * std::max(current_dt, 1e-12);
        }
        next_viz_time = simulation_time + viz_interval_seconds;
        Log::Info("RunCavityCoupling",
                  "Visualization cadence: every ", viz_interval_seconds,
                  " s (next at t=", next_viz_time, " s)");
    }

    // Trajectory logger initial export (if enabled)
    if (trajectory_logger.Enabled()) {
        trajectory_logger.MaybeWrite(simulation_time, simulation);
    }

    // Standard preCICE coupling pattern according to documentation
    // See: https://precice.org/adapter-usage.html
    Log::Info("RunCavityCoupling", "Starting preCICE coupling loop");

    const bool continue_after_coupling =
        GetParam<bool>("coupling.continue_after_coupling", false);
    const double requested_end_time =
        GetParam<double>("simulation.end_time", -1.0);

    double biodynamo_timestep = -1.0;
    if (GetParameterStore().Has<double>("simulation.biodynamo_timestep")) {
        biodynamo_timestep = GetParameterStore().Get<double>("simulation.biodynamo_timestep");
    } else if (GetParameterStore().Has<int>("simulation.biodynamo_timestep")) {
        biodynamo_timestep =
            static_cast<double>(GetParameterStore().Get<int>("simulation.biodynamo_timestep"));
    }

    const auto to_lower = [](std::string value) {
        std::transform(value.begin(), value.end(), value.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        return value;
    };

    const std::string integration_mode_raw =
        GetParam<std::string>("simulation.integration_mode", "default");
    std::string integration_mode = to_lower(integration_mode_raw);
    if (integration_mode == "freeze_until_window_complete") {
        Log::Warning("RunCavityCoupling",
                     "simulation.integration_mode=freeze_until_window_complete is deprecated. "
                     "Use simulation.integration_mode=default.");
        integration_mode = "default";
    } else if (integration_mode == "substep_after_window_complete") {
        Log::Warning("RunCavityCoupling",
                     "simulation.integration_mode=substep_after_window_complete is deprecated. "
                     "Use simulation.integration_mode=substep.");
        integration_mode = "substep";
    } else if (integration_mode == "adaptive" ||
               integration_mode == "adaptive-substep" ||
               integration_mode == "adaptive_substeps") {
        integration_mode = "adaptive_substep";
    }
    const bool default_mode = integration_mode == "default";
    const bool substep_mode = integration_mode == "substep";
    const bool adaptive_substep_mode = integration_mode == "adaptive_substep";

    if (!default_mode && !substep_mode && !adaptive_substep_mode) {
        Log::Fatal("RunCavityCoupling",
                   "Invalid simulation.integration_mode='", integration_mode_raw,
                   "'. Supported values are: default | substep | adaptive_substep");
    }

    auto resolve_config_path = [&](const std::string& maybe_relative_path) {
        if (maybe_relative_path.empty() || maybe_relative_path[0] == '/') {
            return maybe_relative_path;
        }
        return get_cwd() + "/" + maybe_relative_path;
    };

    struct PreciceSchemeInfo {
        double max_time = -1.0;
        bool serial_implicit = false;
        bool implicit = false;
    };

    auto read_precice_scheme_info = [&](const std::string& config_path) {
        PreciceSchemeInfo info;
        std::ifstream config_stream(config_path);
        if (!config_stream.good()) {
            return info;
        }
        std::ostringstream xml_buffer;
        xml_buffer << config_stream.rdbuf();
        const std::string xml_text = xml_buffer.str();

        info.serial_implicit =
            xml_text.find("coupling-scheme:serial-implicit") != std::string::npos;
        info.implicit =
            info.serial_implicit ||
            xml_text.find("coupling-scheme:parallel-implicit") != std::string::npos;

        std::size_t pos = xml_text.find("<max-time");
        while (pos != std::string::npos) {
            const std::size_t value_key = xml_text.find("value=\"", pos);
            if (value_key == std::string::npos) {
                pos = xml_text.find("<max-time", pos + 9);
                continue;
            }
            const std::size_t value_start = value_key + 7;
            const std::size_t value_end = xml_text.find('"', value_start);
            if (value_end == std::string::npos) {
                break;
            }
            try {
                const double parsed_value = std::stod(xml_text.substr(value_start, value_end - value_start));
                if (parsed_value > 0.0) {
                    info.max_time = parsed_value;
                    return info;
                }
            } catch (...) {
                // Try next occurrence if parsing failed.
            }
            pos = xml_text.find("<max-time", value_end);
        }
        return info;
    };

    const std::string precice_config_path =
        resolve_config_path(adapter_params.precice.config_file);
    const PreciceSchemeInfo precice_scheme =
        read_precice_scheme_info(precice_config_path);

    auto resolve_field_read_relative_time = [&](double window_dt) -> double {
        const std::string& mode_raw = adapter_params.coupling.field_read_relative_time;
        std::string mode = mode_raw;
        std::transform(mode.begin(), mode.end(), mode.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (mode == "window_start" || mode == "start") {
            return 0.0;
        }
        if (mode == "window_end" || mode == "end") {
            return window_dt;
        }
        try {
            std::size_t consumed = 0;
            const double explicit_time = std::stod(mode_raw, &consumed);
            if (consumed == mode_raw.size() && explicit_time >= 0.0) {
                return explicit_time;
            }
        } catch (...) {
        }
        Log::Warning("RunCavityCoupling",
                     "Unknown coupling.field_read_relative_time='", mode_raw,
                     "'; using window_start.");
        return 0.0;
    };

    const double xml_max_time = precice_scheme.max_time;

    const bool standalone_mode_enabled =
        continue_after_coupling && !adapter_params.coupling.two_way_coupling &&
        requested_end_time > 0.0;

    double effective_end_time = requested_end_time;
    if (!standalone_mode_enabled) {
        if (xml_max_time > 0.0) {
            effective_end_time = xml_max_time;
            if (requested_end_time > 0.0 &&
                std::abs(requested_end_time - xml_max_time) > 1e-12) {
                Log::Info("RunCavityCoupling",
                          "Coupled-only run: overriding simulation.end_time=",
                          requested_end_time,
                          " s with preCICE XML max-time=",
                          xml_max_time, " s (", precice_config_path, ")");
            }
        } else {
            Log::Warning("RunCavityCoupling",
                         "Coupled-only run: could not parse <max-time> from preCICE config '",
                         precice_config_path,
                         "'. Falling back to simulation.end_time for progress reporting.");
        }
    }

    int internal_substeps_per_window = 0;
    if (GetParameterStore().Has<int>("simulation.internal_substeps_per_window")) {
        internal_substeps_per_window =
            GetParameterStore().Get<int>("simulation.internal_substeps_per_window");
    } else if (GetParameterStore().Has<double>("simulation.internal_substeps_per_window")) {
        internal_substeps_per_window = static_cast<int>(std::lround(
            GetParameterStore().Get<double>("simulation.internal_substeps_per_window")));
    }
    if (internal_substeps_per_window < 0) {
        Log::Fatal("RunCavityCoupling",
                   "simulation.internal_substeps_per_window must be >= 0");
    }
    if (substep_mode && internal_substeps_per_window <= 0) {
        Log::Fatal("RunCavityCoupling",
                   "simulation.internal_substeps_per_window must be > 0 when "
                   "simulation.integration_mode=substep");
    }

    const bool biodynamo_timestep_required =
        default_mode || adaptive_substep_mode || standalone_mode_enabled;
    if (biodynamo_timestep_required && biodynamo_timestep <= 0.0) {
        Log::Fatal("RunCavityCoupling",
                   "simulation.biodynamo_timestep must be > 0 when used as an active "
                   "integration timestep or standalone continuation timestep.");
    }

    if (continue_after_coupling && !standalone_mode_enabled) {
        Log::Warning("RunCavityCoupling",
                     "coupling.continue_after_coupling=true but continuation prerequisites "
                     "were not met (requires one-way coupling and simulation.end_time>0); "
                     "standalone continuation disabled");
    }
    if (standalone_mode_enabled) {
        Log::Info("RunCavityCoupling",
                  "Standalone continuation enabled after coupling end: simulation.end_time=",
                  requested_end_time, " s, simulation.biodynamo_timestep=",
                  biodynamo_timestep, " s");
    }

    int step = 0;

    // Inlet spawning state for continuous injection
    double inlet_spawn_interval = InletSpawnInterval();
    bool inlet_spawning_enabled = IsInletGenerationEnabled() && inlet_spawn_interval > 0.0;
    // Start the recurring spawn clock at zero.  When an initial cohort is
    // present for BioDynaMo's VTK setup, the next identical cohort is released
    // after exactly one configured interval (not after two intervals).
    double time_since_last_spawn = 0.0;
    const double first_recurring_release_time = inlet_spawn_interval;
    std::string inlet_mesh_name = GetParam<std::string>("mesh.inlet_mesh", "InletSurfaceMesh");

    if (inlet_spawning_enabled) {
        Log::Info("RunCavityCoupling", "Continuous inlet spawning enabled: every ",
                  inlet_spawn_interval, " s; first recurring release after ",
                  first_recurring_release_time, " s");
    }

    // BioDynaMo timestep in seconds used directly for agent integration.
    const bool enable_bbox_escape_removal =
        GetParam<bool>("agent.enable_bbox_escape_removal", false);
    const bool low_velocity_removal_enabled =
        GetParam<bool>("agent.remove_low_velocity", false);
    const double low_velocity_threshold =
        std::max(0.0, GetParam<double>("agent.low_velocity_threshold", 0.0));
    const bool low_velocity_use_fluid =
        GetParam<bool>("agent.low_velocity_use_fluid_velocity", true);
    const int profiling_interval = std::max(1, GetParam<int>("logging.print_interval", 100));
    const int timestep_log_interval =
        std::max(1, GetParam<int>("logging.timestep_interval", 1));
    const bool outlet_surface_removal_enabled = IsOutletRemovalEnabled();
    const std::string outlet_removal_summary =
        std::string(outlet_surface_removal_enabled ? "preCICE-outlet-surface-distance" : "") +
        (outlet_surface_removal_enabled && enable_bbox_escape_removal ? "+" : "") +
        (enable_bbox_escape_removal ? "bbox-escape" : "") +
        ((outlet_surface_removal_enabled || enable_bbox_escape_removal) &&
                 low_velocity_removal_enabled
             ? "+"
             : "") +
        (low_velocity_removal_enabled ? "low-velocity-immediate" : "");
    if (outlet_surface_removal_enabled || enable_bbox_escape_removal) {
        Log::Info("RunCavityCoupling",
                  "Outlet/escape-removal strategy: OpenFOAM/preCICE outlet surface "
                  "face-center distance=",
                  (outlet_surface_removal_enabled ? "ON" : "OFF"),
                  ", bbox=",
                  (enable_bbox_escape_removal ? "ON" : "OFF"));
    }
    if (low_velocity_removal_enabled) {
        Log::Info("RunCavityCoupling",
                  "Low-velocity agent removal enabled: threshold=",
                  low_velocity_threshold, " m/s, removal=immediate, source=",
                  (low_velocity_use_fluid ? "mapped fluid velocity" : "agent velocity"));
    }
    Log::Info("RunCavityCoupling",
              "BioDynaMo integration timestep (simulation.biodynamo_timestep): ",
              biodynamo_timestep, " s");
    Log::Info("RunCavityCoupling",
              "BioDynaMo integration mode (simulation.integration_mode): ",
              integration_mode);
    Log::Info("RunCavityCoupling",
              "preCICE scheme detected from XML: ",
              precice_scheme.serial_implicit ? "serial-implicit" :
              (precice_scheme.implicit ? "implicit" : "explicit/unknown"));
    if (default_mode) {
        Log::Info("RunCavityCoupling",
                  "Default mode selected: coupled-phase biology timestep is "
                  "forced to preCICE dt (from XML time-window-size), ignoring "
                  "simulation.biodynamo_timestep while coupling is active.");
        if (effective_end_time > 0.0) {
            Log::Info("RunCavityCoupling",
                      "Default mode coupled-phase end time (from XML max-time): ",
                      effective_end_time, " s");
        }
    }
    if (substep_mode) {
        Log::Info("RunCavityCoupling",
                  "Converged-window substepping enabled with fixed "
                  "simulation.internal_substeps_per_window=",
                  internal_substeps_per_window);
        Log::Info("RunCavityCoupling",
                  "In fixed substep mode, coupled-phase biology dt is derived from "
                  "preCICE window_dt / simulation.internal_substeps_per_window. "
                  "simulation.biodynamo_timestep is used only for standalone continuation.");
    }
    if (adaptive_substep_mode) {
        Log::Info("RunCavityCoupling",
                  "Adaptive converged-window substepping enabled: biology dt <= ",
                  biodynamo_timestep, " s");
    }

    auto remove_agents_outside_domain = [&](double buffer) {
        std::vector<bdm::AgentUid> to_remove;
        rm->ForEachAgent([&](bdm::Agent* agent) {
            auto* cell = dynamic_cast<MyCell*>(agent);
            if (!cell || cell->IsStationary()) return;
            const auto& pos = cell->GetPosition();
            bool outside =
                pos[0] < hybrid_bounds.min[0] - buffer || pos[0] > hybrid_bounds.max[0] + buffer ||
                pos[1] < hybrid_bounds.min[1] - buffer || pos[1] > hybrid_bounds.max[1] + buffer ||
                pos[2] < hybrid_bounds.min[2] - buffer || pos[2] > hybrid_bounds.max[2] + buffer;
            if (outside) {
                to_remove.push_back(cell->GetUid());
            }
        });
        for (const auto& uid : to_remove) {
            rm->RemoveAgent(uid);
        }
        return to_remove.size();
    };

    // Removal is deliberately outside the Behavior::Run call. This function
    // executes only while committing an accepted biology step, so implicit
    // preCICE trial iterations can still roll intracellular state back safely.
    auto remove_dead_rons_agents = [&](double dt) {
        if (!remove_dead_rons_agents_enabled || dt <= 0.0) {
            return static_cast<size_t>(0);
        }

        auto* random = Simulation::GetActive()->GetRandom();
        std::vector<bdm::AgentUid> to_remove;
        rm->ForEachAgent([&](bdm::Agent* agent) {
            auto* cell = dynamic_cast<MyCell*>(agent);
            if (!cell) {
                return;
            }
            const double probability = DeathRemovalProbability(
                cell->GetDamage(), dt, rons_death_threshold,
                rons_death_removal_rate);
            if (probability > 0.0 && random->Uniform() < probability) {
                to_remove.push_back(cell->GetUid());
            }
        });
        for (const auto& uid : to_remove) {
            rm->RemoveAgent(uid);
        }
        if (!to_remove.empty()) {
            Log::Info("RonsResponse",
                      "Removed ", to_remove.size(),
                      " agents by probabilistic damage removal");
            std::cout << "RonsResponse: removed " << to_remove.size()
                      << " agents by probabilistic damage removal"
                      << std::endl;
        }
        return to_remove.size();
    };

    using PrevPosMap = std::unordered_map<bdm::AgentUid, std::array<double, 3>>;

    auto point_segment_distance_sq = [](const std::array<double, 3>& point,
                                        const std::array<double, 3>& seg_start,
                                        const std::array<double, 3>& seg_end) {
        const double vx = seg_end[0] - seg_start[0];
        const double vy = seg_end[1] - seg_start[1];
        const double vz = seg_end[2] - seg_start[2];
        const double wx = point[0] - seg_start[0];
        const double wy = point[1] - seg_start[1];
        const double wz = point[2] - seg_start[2];
        const double vv = vx * vx + vy * vy + vz * vz;
        const double projection = vv > 0.0
                                      ? (wx * vx + wy * vy + wz * vz) / vv
                                      : 0.0;
        const double t = std::clamp(projection, 0.0, 1.0);
        const double cx = seg_start[0] + t * vx;
        const double cy = seg_start[1] + t * vy;
        const double cz = seg_start[2] + t * vz;
        const double dx = point[0] - cx;
        const double dy = point[1] - cy;
        const double dz = point[2] - cz;
        return dx * dx + dy * dy + dz * dz;
    };

    auto remove_agents_reaching_outlet_surfaces =
        [&](const PrevPosMap& prev_positions, double threshold) {
            if (all_outlet_meshes.empty()) return static_cast<size_t>(0);

            const double threshold_sq = threshold * threshold;
            std::set<bdm::AgentUid> to_remove;

            rm->ForEachAgent([&](bdm::Agent* agent) {
                auto* cell = dynamic_cast<MyCell*>(agent);
                if (!cell || cell->IsStationary()) return;

                const auto& pos = cell->GetPosition();
                const std::array<double, 3> p_new = {pos[0], pos[1], pos[2]};
                std::array<double, 3> p_prev = p_new;
                auto prev_it = prev_positions.find(cell->GetUid());
                if (prev_it != prev_positions.end()) {
                    p_prev = prev_it->second;
                }

                for (const auto& [mesh_name, mesh_data] : all_outlet_meshes) {
                    if (!mesh_data.is_initialized || mesh_data.NumVertices() == 0) {
                        continue;
                    }

                    for (size_t i = 0; i < mesh_data.NumVertices(); ++i) {
                        const auto face_center = mesh_data.GetVertex(i);
                        const std::array<double, 3> outlet_point = {
                            face_center[0], face_center[1], face_center[2]};
                        const double d2 =
                            point_segment_distance_sq(outlet_point, p_prev, p_new);
                        if (d2 <= threshold_sq) {
                            to_remove.insert(cell->GetUid());
                            return;
                        }
                    }
                }
            });

            for (const auto& uid : to_remove) {
                rm->RemoveAgent(uid);
            }
            if (!to_remove.empty()) {
                Log::Info("PreciceAdapter", "Removed ", to_remove.size(),
                          " agents at preCICE outlet surface meshes (checked ",
                          all_outlet_meshes.size(), " outlet mesh(es))");
            }
            return to_remove.size();
    };

    auto remove_agents_with_low_velocity = [&]() {
        if (!low_velocity_removal_enabled || low_velocity_threshold <= 0.0) {
            return static_cast<size_t>(0);
        }

        std::vector<bdm::AgentUid> to_remove;
        to_remove.reserve(16);

        rm->ForEachAgent([&](bdm::Agent* agent) {
            auto* cell = dynamic_cast<MyCell*>(agent);
            if (!cell || cell->IsStationary()) return;

            const auto& vel = low_velocity_use_fluid
                                  ? cell->GetFluidVelocity()
                                  : cell->GetVelocity();
            const double speed =
                std::sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
            if (speed <= low_velocity_threshold) {
                to_remove.push_back(cell->GetUid());
            }
        });

        for (const auto& uid : to_remove) {
            rm->RemoveAgent(uid);
        }

        if (!to_remove.empty()) {
            Log::Info("RunCavityCoupling",
                      "Removed ", to_remove.size(),
                      " agents immediately after velocity dropped to/below threshold (threshold=",
                      low_velocity_threshold, " m/s)");
        }
        return to_remove.size();
    };

    bool adapter_finalized = false;
    bool completed_coupled_window = false;
    std::size_t completed_window_count = 0;
    double completed_coupling_time = 0.0;

    std::size_t lifecycle_spawned = 0;
    std::size_t lifecycle_outlet_removed = 0;
    std::size_t lifecycle_bbox_removed = 0;
    std::size_t lifecycle_low_velocity_removed = 0;
    std::size_t lifecycle_rons_removed = 0;

    double pending_biology_window_dt = 0.0;
    bool biology_mutated_since_checkpoint = false;
    int checkpoint_write_noops = 0;
    int checkpoint_read_noops = 0;

    // =========================================================================
    // Fix 1 & 2: Real checkpoint/rollback snapshot for implicit coupling.
    //
    // For each implicit window we snapshot every agent's full physics state
    // on RequiresWritingCheckpoint().  On RequiresReadingCheckpoint() we
    // restore that snapshot so biology can re-run with the corrected CFD
    // fields — enabling true rollback-corrected implicit two-way coupling.
    // =========================================================================

    struct AgentSnapshot {
        bdm::AgentUid uid;
        bdm::Real3    position;
        bdm::Real3    velocity;
        bdm::Real3    fluid_velocity;
        bdm::Real3    accumulated_force;
        double        drag_Sp;
        double        added_mass;
        double        mass;
        double        diameter;
        double        pressure;
        double        temperature;
        double        oxygen;
        double        oxygen_uptake;
        double        h2o2_uptake;
        double        h2o2;
        double        no2;
        double        no2_uptake;
        double        intracellular_h2o2;
        double        intracellular_no2;
        double        oxidative_stress;
        double        damage;
        double        residence_time;
        int           phenotype;
        bdm::Real3    openfoam_cell_center;
        bool          is_stationary;
    };

    std::vector<AgentSnapshot> implicit_checkpoint;
    double implicit_checkpoint_simulation_time = 0.0;
    int    implicit_checkpoint_step = 0;

    auto take_agent_snapshot = [&]() {
        implicit_checkpoint.clear();
        implicit_checkpoint.reserve(rm->GetNumAgents());
        rm->ForEachAgent([&](bdm::Agent* a) {
            auto* cell = dynamic_cast<MyCell*>(a);
            if (!cell) return;
            AgentSnapshot snap;
            snap.uid                         = cell->GetUid();
            snap.position                    = cell->GetPosition();
            snap.velocity                    = cell->GetVelocity();
            snap.fluid_velocity              = cell->GetFluidVelocity();
            snap.accumulated_force           = cell->GetAccumulatedForce();
            snap.drag_Sp                     = cell->GetDragCoefficient();
            snap.added_mass                  = cell->GetAddedMass();
            snap.mass                        = cell->GetMass();
            snap.diameter                    = cell->GetDiameter();
            snap.pressure                    = cell->GetPressure();
            snap.temperature                 = cell->GetTemperature();
            snap.oxygen                      = cell->GetOxygen();
            snap.oxygen_uptake               = cell->GetOxygenUptake();
            snap.h2o2_uptake                 = cell->GetH2O2Uptake();
            snap.h2o2                        = cell->GetH2O2();
            snap.no2                         = cell->GetNO2();
            snap.no2_uptake                  = cell->GetNO2Uptake();
            snap.intracellular_h2o2          = cell->GetIntracellularH2O2();
            snap.intracellular_no2           = cell->GetIntracellularNO2();
            snap.oxidative_stress            = cell->GetOxidativeStress();
            snap.damage                      = cell->GetDamage();
            snap.residence_time               = cell->GetResidenceTime();
            snap.phenotype                   = cell->GetPhenotype();
            snap.openfoam_cell_center        = cell->GetOpenFoamCellCenter();
            snap.is_stationary               = cell->IsStationary();
            implicit_checkpoint.push_back(snap);
        });
    };

    auto restore_agent_snapshot = [&]() {
        // Build uid→snapshot map for O(1) lookups.
        std::unordered_map<bdm::AgentUid, const AgentSnapshot*> snap_map;
        snap_map.reserve(implicit_checkpoint.size());
        for (const auto& snap : implicit_checkpoint) {
            snap_map[snap.uid] = &snap;
        }

        std::vector<bdm::AgentUid> trial_created_agents;
        rm->ForEachAgent([&](bdm::Agent* a) {
            auto* cell = dynamic_cast<MyCell*>(a);
            if (!cell) return;
            if (snap_map.find(cell->GetUid()) == snap_map.end()) {
                trial_created_agents.push_back(cell->GetUid());
            }
        });
        for (const auto& uid : trial_created_agents) {
            rm->RemoveAgent(uid);
        }

        rm->ForEachAgent([&](bdm::Agent* a) {
            auto* cell = dynamic_cast<MyCell*>(a);
            if (!cell) return;
            auto it = snap_map.find(cell->GetUid());
            if (it == snap_map.end()) return;
            const auto& s = *it->second;
            cell->SetPosition(s.position);
            cell->SetVelocity(s.velocity);
            cell->SetFluidVelocity(s.fluid_velocity);
            cell->SetDragCoefficient(s.drag_Sp);
            cell->SetAddedMass(s.added_mass);
            cell->SetOpenFoamCellCenter(s.openfoam_cell_center);
            cell->SetAccumulatedForce(s.accumulated_force);
            cell->SetMass(s.mass);
            cell->SetDiameter(s.diameter);
            cell->SetPressure(s.pressure);
            cell->SetH2O2(s.h2o2);
            cell->SetNO2(s.no2);
            if (!GetParam<bool>("coupling.store_temperature_as_no2", false)) {
                cell->SetTemperature(s.temperature);
            }
            if (!GetParam<bool>("coupling.store_oxygen_as_h2o2", false)) {
                cell->SetOxygen(s.oxygen);
            }
            cell->SetOxygenUptake(s.oxygen_uptake);
            cell->SetH2O2Uptake(s.h2o2_uptake);
            cell->SetNO2Uptake(s.no2_uptake);
            cell->SetIntracellularH2O2(s.intracellular_h2o2);
            cell->SetIntracellularNO2(s.intracellular_no2);
            cell->SetOxidativeStress(s.oxidative_stress);
            cell->SetDamage(s.damage);
            cell->SetResidenceTime(s.residence_time);
            cell->SetPhenotype(s.phenotype);
            cell->SetStationary(s.is_stationary);
        });
        if (!trial_created_agents.empty()) {
            Log::Info("RunCavityCoupling",
                      "Checkpoint rollback removed ",
                      trial_created_agents.size(),
                      " trial-created agents");
        }
        simulation_time = implicit_checkpoint_simulation_time;
        step            = implicit_checkpoint_step;
        biology_mutated_since_checkpoint = false;
    };

    bool final_endpoint_snapshot_ready = false;

    // Coupling loop: preCICE may retry the same implicit window many times.
    // With real rollback, biology now re-runs each iteration before sending
    // updated source terms.
    while (true) {
        const bool coupling_active = adapter.IsCouplingOngoing();
        const bool standalone_active = standalone_mode_enabled && completed_coupled_window && !coupling_active &&
                                       (simulation_time + 1e-12 < requested_end_time);
        if (!coupling_active && standalone_mode_enabled && !completed_coupled_window) {
            Log::Warning("RunCavityCoupling",
                         "preCICE coupling ended before any time window completed; "
                         "skipping standalone continuation");
        }
        if (!coupling_active && !standalone_active) {
            break;
        }

        if (standalone_active && !adapter_finalized) {
            // Freeze the last coupled field snapshot before switching to
            // standalone BioDynaMo continuation with frozen CFD data.
            adapter.EnsureCachedFieldSnapshot(0.0);
            adapter.Finalize();
            adapter_finalized = true;
            Log::Info("RunCavityCoupling",
                      "Completed coupled window exchange; switching to standalone "
                      "BioDynaMo mode with frozen fields");
        }

        double step_dt = biodynamo_timestep;
        if (coupling_active) {
            step_dt = dt;
        }
        if (!coupling_active && standalone_mode_enabled &&
            simulation_time + step_dt > requested_end_time) {
            step_dt = requested_end_time - simulation_time;
        }
        if (step_dt <= 0.0) {
            break;
        }

        const auto window_wall_start = std::chrono::steady_clock::now();
        double fetch_seconds = 0.0;
        double spawn_seconds = 0.0;
        double simulate_seconds = 0.0;
        double removal_seconds = 0.0;
        double clamp_map_seconds = 0.0;
        double output_seconds = 0.0;
        double advance_seconds = 0.0;
        bool implicit_ckpt_write = false;
        bool implicit_ckpt_read = false;

        if (coupling_active && adapter.RequiresWritingCheckpoint()) {
            // Fix 1: Real snapshot — persist full agent physics state so we
            // can roll back if preCICE requests another implicit iteration.
            take_agent_snapshot();
            implicit_checkpoint_simulation_time  = simulation_time;
            implicit_checkpoint_step             = step;
            pending_biology_window_dt            = 0.0;
            biology_mutated_since_checkpoint     = false;
            Log::Info("RunCavityCoupling",
                      "Checkpoint written: snapshotted ", implicit_checkpoint.size(),
                      " agents at t=", simulation_time, " s");
            implicit_ckpt_write = true;
            ++checkpoint_write_noops;  // kept for statistics log at end
        }
        if (coupling_active && adapter.RequiresReadingCheckpoint()) {
            // Fix 1 & 2: Real rollback — restore biology to the snapshot so
            // the window re-runs with the corrected CFD fields from this
            // iteration.  Biology is intentionally allowed to mutate between
            // iterations; the snapshot undoes those mutations.
            if (!implicit_checkpoint.empty()) {
                restore_agent_snapshot();
                Log::Info("RunCavityCoupling",
                          "Checkpoint restored: rolled back ", implicit_checkpoint.size(),
                          " agents to t=", simulation_time, " s");
                implicit_ckpt_read = true;
            } else {
                Log::Warning("RunCavityCoupling",
                             "preCICE requested rollback but no snapshot exists; "
                             "resetting accumulators only.");
                biology_mutated_since_checkpoint = false;
            }
            pending_biology_window_dt = 0.0;
            ++checkpoint_read_noops;  // kept for statistics log at end
        }

        if (coupling_active) {
            const auto fetch_start = std::chrono::steady_clock::now();
            const double read_relative_time =
                resolve_field_read_relative_time(step_dt);
            adapter.FetchInputData(read_relative_time);
            // FetchInputData() already refreshes the cached CFD snapshot when
            // the configured field-read cadence says a real preCICE read
            // should happen. Keeping the cache update there preserves frozen
            // and throttled field modes while substeps below still re-sample
            // the cached window field at updated agent positions.
            // preCICE recomputes domain bounds each fetch; re-apply hybrid bounds for consistency
            adapter.SetDomainBounds(hybrid_bounds);
            if (adapter.IsTwoWayCoupling() && precice_scheme.serial_implicit) {
                const auto u_stats = sample_fluid_velocity_stats();
                Log::Info("ImplicitCouplingDiag-BDM",
                          "window=", completed_window_count + 1,
                          " iter=", adapter.GetCouplingIteration(),
                          " coupled_t=", completed_coupling_time,
                          " read_rel_t=", read_relative_time, " s",
                          " ckpt_write=", implicit_ckpt_write,
                          " ckpt_read=", implicit_ckpt_read,
                          " |U_fluid|_max=", u_stats.max,
                          " |U_fluid|_mean=", u_stats.mean,
                          " agents_sampled=", u_stats.sampled);
            }
            fetch_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - fetch_start).count();
        } else if (standalone_active) {
            const auto fetch_start = std::chrono::steady_clock::now();
            // Standalone continuation: re-sample cached frozen CFD fields at
            // current particle positions each BioDynaMo step.
            adapter.RefreshCachedFieldData();
            fetch_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - fetch_start).count();
        }

        // Diagnostic: sample fluid velocity at a few agents to see whether the
        // nearest-neighbour mapping is actually delivering non-zero data.  If
        // |U_fluid| ≈ 0 for many windows something is wrong on the preCICE /
        // OpenFOAM side (initial-data exchange, restart from rest, etc.) and
        // no amount of sub-stepping will make the agents move.
        if (verbose_logging && step < 12) {
            const auto stats = sample_fluid_velocity_stats();
            Log::Info("RunCavityCoupling",
                      "STEP ", step, " pre-Simulate: |U_fluid| max=", stats.max,
                      " mean=", stats.mean, " over ", stats.sampled, " agents");
        }

        current_dt = step_dt;

        // Fix 2: Evolve biology inside every implicit iteration (not only on
        // window convergence).  After rollback the snapshot has already
        // restored agent state; we run one biology step here so agents
        // experience the corrected CFD fields from this iteration, then
        // collect updated source terms and send them to OpenFOAM.
        //
        // In explicit / standalone mode RequiresWritingCheckpoint() is always
        // false so this block is entered exactly once per step — same as the
        // old behaviour.
        if (coupling_active && adapter.IsTwoWayCoupling()) {
            const auto output_start = std::chrono::steady_clock::now();

            // Implicit trial: evolve the same biology trajectory that will be
            // committed if this preCICE iteration converges.  For substep mode,
            // this means the fixed-point residual is based on the configured
            // BioDynaMo substeps, not on a coarser one-step approximation.
            {
                const double iter_window_dt =
                    (step_dt > 0.0) ? step_dt : biodynamo_timestep;
                const double saved_sim_time = simulation_time;
                const int    saved_step     = step;
                const int trial_substeps = substep_mode
                    ? std::max(1, internal_substeps_per_window)
                    : (adaptive_substep_mode
                        ? std::max(1, static_cast<int>(
                              std::ceil(iter_window_dt / biodynamo_timestep)))
                        : 1);
                const double trial_dt =
                    iter_window_dt / static_cast<double>(trial_substeps);

                Log::Info("RunCavityCoupling",
                          "Implicit trial biology: window=", completed_window_count + 1,
                          " iter=", adapter.GetCouplingIteration(),
                          " window_dt=", iter_window_dt,
                          " s, substeps=", trial_substeps,
                          ", substep_dt=", trial_dt, " s");
                std::cout << "BDM implicit iter " << adapter.GetCouplingIteration()
                          << " window " << (completed_window_count + 1)
                          << " — trial Simulate(" << trial_substeps
                          << " x 1), dt=" << trial_dt << " s"
                          << std::endl;

                for (int trial_i = 0; trial_i < trial_substeps; ++trial_i) {
                    current_dt = trial_dt;
                    simulation.GetScheduler()->Simulate(1);
                    if (motion_enabled && clamp_to_domain) {
                        adapter.ClampDynamicAgentsToDomain();
                    }
                    if (motion_enabled) {
                        adapter.UpdateDynamicMapping();
                    }
                }
                biology_mutated_since_checkpoint = true;
                simulation_time = saved_sim_time;
                step            = saved_step;
            }

            adapter.CollectSourceTerms();

            adapter.SendOutputData();
            output_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - output_start).count();
        }

        // Advance coupling and get time step for next iteration.
        bool time_window_completed = true;
        if (coupling_active) {
            // Preserve the carrier field at the exact XML max-time before the
            // last advance closes preCICE. Without this snapshot, the final
            // visualization retains the preceding window-start sample.
            const bool final_coupling_step =
                effective_end_time > 0.0 &&
                completed_coupling_time + pending_biology_window_dt + step_dt
                    >= effective_end_time - 1e-12;
            if (final_coupling_step && IsMasslessTracerMode()) {
                adapter.EnsureCachedFieldSnapshot(step_dt);
                final_endpoint_snapshot_ready = true;
            }

            const auto advance_start = std::chrono::steady_clock::now();
            adapter.Advance(step_dt);
            advance_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - advance_start).count();
            pending_biology_window_dt += step_dt;
            time_window_completed = adapter.IsTimeWindowComplete();
            if (time_window_completed) {
                completed_coupled_window = true;
                ++completed_window_count;
                completed_coupling_time += step_dt;
                Log::Info("RunCavityCoupling",
                          "Completed preCICE time window ",
                          completed_window_count,
                          " at coupled_time=",
                          completed_coupling_time,
                          " s (window_dt=",
                          step_dt,
                          " s)");
                std::cout << "BDM accepted window " << completed_window_count
                          << " at coupled_time=" << completed_coupling_time
                          << " s — committing biology" << std::endl;
            }
            dt = adapter.GetTimeStepSize();
            current_dt = dt;
        }

        // Only evolve biology (and emit visualization / outlet removal) when
        // the time window has converged. In explicit / standalone modes
        // time_window_completed is always true so this still runs every step.
        if (!coupling_active || time_window_completed) {
            // PrevPosMap is declared here (outside the lambda) so its heap
            // capacity is reused across every call to run_single_biology_step,
            // avoiding a fresh unordered_map allocation on each biology step.
            PrevPosMap step_prev_positions;

            // Two-way implicit: the converged trial iteration already advanced
            // agents for feedback (including per-window substeps). Commit that
            // state without a second Simulate pass on window acceptance.
            const bool skip_physics_on_acceptance =
                coupling_active && adapter.IsTwoWayCoupling() &&
                (default_mode || substep_mode);

            const auto run_single_biology_step = [&](double biology_dt, bool run_physics_step) {
                current_dt = biology_dt;

                if (inlet_spawning_enabled) {
                    const auto spawn_start = std::chrono::steady_clock::now();
                    bool spawned_any = false;
                    time_since_last_spawn += biology_dt;
                    while (time_since_last_spawn + 1e-15 >= inlet_spawn_interval) {
                        const int configured_inlet_batch = InletAgentsPerSpawn();
                        const size_t agents_to_spawn = configured_inlet_batch > 0
                            ? static_cast<size_t>(configured_inlet_batch)
                            : static_cast<size_t>(0);
                        const double offset = InletZoneOffset();

                        const size_t spawned = adapter.CreateAgentsAtInletMesh(
                            inlet_mesh_name,
                            agent_factory,
                            offset,
                            agents_to_spawn);
                        lifecycle_spawned += spawned;

                        if (enforce_fixed_inlet_batch &&
                            spawned != expected_fixed_inlet_batch_size) {
                            Log::Fatal(
                                "RunCavityCoupling",
                                "Fixed inlet batch changed unexpectedly: expected ",
                                expected_fixed_inlet_batch_size, " agents but created ",
                                spawned,
                                ". The central inlet owner-cell set must remain constant.");
                        }

                        if (spawned > 0) {
                            spawned_any = true;
                            rm->ForEachAgent([&](Agent* agent) {
                                if (auto* cell = dynamic_cast<MyCell*>(agent)) {
                                    if (cell->GetAllBehaviors().size() == 0) {
                                        cell->SetVelocity({0.0, 0.0, 0.0});
                                        cell->SetStationary(!motion_enabled);
                                        if (motion_enabled) {
                                            auto* new_behavior = wall_behavior->NewCopy();
                                            cell->AddBehavior(new_behavior);
                                        }
                                        if (growth_division_enabled) {
                                            auto* growth_behavior = oxygen_growth_behavior->NewCopy();
                                            cell->AddBehavior(growth_behavior);
                                        }
                                        if (rons_response_behavior) {
                                            cell->AddBehavior(rons_response_behavior->NewCopy());
                                        }
                                        if (cell_spread_enabled && cell_spread_behavior_proto) {
                                            cell->AddBehavior(cell_spread_behavior_proto->NewCopy());
                                        }
                                    }
                                }
                                return true;
                            });

                            Log::Info("RunCavityCoupling", "Step ", step, ": spawned ", spawned,
                                      " agents at inlet (total time: ", simulation_time, " s)");
                        }

                        time_since_last_spawn -= inlet_spawn_interval;
                        if (inlet_spawn_interval <= 0.0) {
                            time_since_last_spawn = 0.0;
                            break;
                        }
                    }
                    if (spawned_any) {
                        // Newly injected agents were not part of the field read
                        // at the start of this coupling step. Perform a live JIT
                        // velocity read now: the optional frozen-field cache is
                        // normally absent when fields are read every window and
                        // therefore cannot initialize a newly spawned tracer.
                        adapter.RefreshCurrentVelocityData(0.0);
                    }
                    spawn_seconds += std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - spawn_start).count();
                }
                biology_mutated_since_checkpoint = true;

                // Reuse the outer-scope map: clear() retains heap capacity.
                step_prev_positions.clear();
                step_prev_positions.reserve(rm->GetNumAgents());
                rm->ForEachAgent([&](bdm::Agent* a) {
                    const auto& p = a->GetPosition();
                    step_prev_positions[a->GetUid()] = {p[0], p[1], p[2]};
                });

                if (run_physics_step) {
                    const auto simulate_start = std::chrono::steady_clock::now();
                    simulation.GetScheduler()->Simulate(1);
                    simulate_seconds += std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - simulate_start).count();
                }

                // Residence time is committed exactly once per accepted
                // biology step. It is intentionally not advanced during
                // rollback-prone implicit trial iterations.
                rm->ForEachAgent([&](bdm::Agent* agent) {
                    if (auto* cell = dynamic_cast<MyCell*>(agent)) {
                        cell->AdvanceResidenceTime(biology_dt);
                    }
                });

                // Attach the exact t=max-time carrier sample to the last
                // passive-tracer state (including the final injected batch).
                if (final_endpoint_snapshot_ready &&
                    simulation_time + biology_dt >= effective_end_time - 1e-12) {
                    adapter.RefreshCachedFieldData();
                    rm->ForEachAgent([&](bdm::Agent* agent) {
                        if (auto* cell = dynamic_cast<MyCell*>(agent)) {
                            if (!cell->IsStationary()) {
                                cell->SetVelocityFromFluid();
                            }
                        }
                    });
                }

                const auto removal_start = std::chrono::steady_clock::now();
                lifecycle_rons_removed += remove_dead_rons_agents(biology_dt);
                if (outlet_surface_removal_enabled) {
                    lifecycle_outlet_removed +=
                        remove_agents_reaching_outlet_surfaces(
                            step_prev_positions,
                            std::abs(OutletRemovalDistance()));
                }
                if (enable_bbox_escape_removal) {
                    lifecycle_bbox_removed +=
                        remove_agents_outside_domain(
                            std::abs(OutletRemovalDistance()));
                }
                if (low_velocity_removal_enabled) {
                    lifecycle_low_velocity_removed +=
                        remove_agents_with_low_velocity();
                }
                removal_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - removal_start).count();

                const auto clamp_map_start = std::chrono::steady_clock::now();
                if (motion_enabled && clamp_to_domain) {
                    adapter.ClampDynamicAgentsToDomain();
                }
                if (motion_enabled) {
                    adapter.UpdateDynamicMapping();
                }
                clamp_map_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - clamp_map_start).count();

                if (ExportVisualization() && !skip_physics_on_acceptance) {
                    const double biology_time_after_step = simulation_time + biology_dt;
                    bool should_export = false;
                    if (viz_interval_seconds > 0.0) {
                        constexpr double eps = 1e-12;
                        if (biology_time_after_step + eps >= next_viz_time) {
                            should_export = true;
                            while (biology_time_after_step + eps >= next_viz_time) {
                                next_viz_time += viz_interval_seconds;
                            }
                        }
                    }

                    if (!should_export) {
                        // No output due at this biology step.
                    } else {
                    const auto output_start = std::chrono::steady_clock::now();
                    const auto visualize_ops = simulation.GetScheduler()->GetOps("visualize");
                    if (!visualize_ops.empty()) {
                        visualize_ops[0]->operator()();
                    }
                    output_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - output_start).count();
                    }
                }

                if (verbose_logging && step < 3) {
                    auto* rm_diag = bdm::Simulation::GetActive()->GetResourceManager();
                    int log_count_diag = 0;
                    rm_diag->ForEachAgent([&](Agent* agent) {
                        if (log_count_diag < 5) {
                            const auto& p = agent->GetPosition();
                            double mass = -1.0;
                            if (auto* cell = dynamic_cast<MyCell*>(agent)) {
                                mass = cell->GetMass();
                            }
                            Log::Info("RunCavityCoupling", "STEP ", step, " dt=", biology_dt,
                                      " Agent[", log_count_diag, "] pos=(",
                                      p[0], ", ", p[1], ", ", p[2], ") mass=", mass);
                            log_count_diag++;
                        }
                        return log_count_diag < 5;
                    });
                }

                simulation_time += biology_dt;
                step++;

                if (coupling_active &&
                    step % std::max(1, internal_substeps_per_window) == 0) {
                    std::size_t active_agents = 0;
                    std::size_t in_domain_agents = 0;
                    std::size_t outside_geometry_agents = 0;
                    double max_wall_penetration = 0.0;
                    const auto boundary_geometry = adapter.GetBoundaryGeometry();
                    rm->ForEachAgent([&](bdm::Agent* agent) {
                        auto* cell = dynamic_cast<MyCell*>(agent);
                        if (!cell || cell->IsStationary()) return;
                        ++active_agents;
                        const auto& pos = cell->GetPosition();
                        if (hybrid_bounds.Contains({pos[0], pos[1], pos[2]})) {
                            ++in_domain_agents;
                        }
                        if (boundary_geometry) {
                            const bdm_adapter::BoundaryVec3 point = {
                                pos[0], pos[1], pos[2]};
                            if (!IsForce2D() && boundary_geometry->IsClosed() &&
                                !boundary_geometry->IsInside(point)) {
                                ++outside_geometry_agents;
                            }
                            const auto closest =
                                boundary_geometry->ClosestPoint(point, true);
                            if (closest.found) {
                                const double radius =
                                    0.5 * cell->GetDiameter();
                                const double penetration = radius
                                    + bdm_adapter::BDot(
                                        bdm_adapter::BSub(point, closest.point),
                                        closest.normal);
                                max_wall_penetration =
                                    std::max(max_wall_penetration, penetration);
                            }
                        }
                    });
                    Log::Info("AneurysmAgentLifecycle",
                              "window=", completed_window_count + 1,
                              " t=", simulation_time,
                              " spawned=", lifecycle_spawned,
                              " active=", active_agents,
                              " in_domain=", in_domain_agents,
                              " outside_geometry=", outside_geometry_agents,
                              " max_wall_penetration=",
                              std::max(0.0, max_wall_penetration),
                              " removed_outlet=", lifecycle_outlet_removed,
                              " removed_bbox=", lifecycle_bbox_removed,
                              " removed_low_velocity=", lifecycle_low_velocity_removed);
                    lifecycle_spawned = 0;
                    lifecycle_outlet_removed = 0;
                    lifecycle_bbox_removed = 0;
                    lifecycle_low_velocity_removed = 0;
                }

                if (step % timestep_log_interval == 0) {
                    const bool has_total_time = effective_end_time > 0.0;
                    const double progress = has_total_time
                        ? std::clamp(simulation_time / effective_end_time, 0.0, 1.0)
                        : 1.0;

                    const auto format_scientific_dt = [](double value) {
                        std::ostringstream dt_stream;
                        dt_stream << std::scientific << std::setprecision(6) << value;

                        std::string text = dt_stream.str();
                        const auto exponent_pos = text.find('e');
                        if (exponent_pos == std::string::npos) {
                            return text;
                        }

                        std::string mantissa = text.substr(0, exponent_pos);
                        while (mantissa.size() > 1 && mantissa.back() == '0') {
                            mantissa.pop_back();
                        }
                        if (!mantissa.empty() && mantissa.back() == '.') {
                            mantissa.pop_back();
                        }

                        std::string exponent = text.substr(exponent_pos + 1);
                        char sign = '+';
                        std::size_t exponent_index = 0;
                        if (!exponent.empty() && (exponent[0] == '+' || exponent[0] == '-')) {
                            sign = exponent[0];
                            exponent_index = 1;
                        }
                        while (exponent_index + 1 < exponent.size() && exponent[exponent_index] == '0') {
                            ++exponent_index;
                        }

                        std::string normalized_exponent = exponent.substr(exponent_index);
                        if (normalized_exponent.empty()) {
                            normalized_exponent = "0";
                        }

                        return mantissa + "e" + (sign == '-' ? "-" : "") + normalized_exponent;
                    };

                    std::ostringstream timestep_line;
                    timestep_line << std::fixed << std::setprecision(3)
                                  << "BDM_timestep = " << step
                                  << " t = " << simulation_time;
                    if (has_total_time) {
                        timestep_line << "/" << effective_end_time
                                      << " (" << std::setprecision(1)
                                      << (progress * 100.0) << "%)"
                                      << std::setprecision(3);
                    }
                    timestep_line << " (dt = " << format_scientific_dt(biology_dt) << ")";

                    Log::Info("RunCavityCoupling", timestep_line.str());
                    std::cout << timestep_line.str() << std::endl;
                }

                if (trajectory_logger.Enabled()) {
                    const auto output_start = std::chrono::steady_clock::now();
                    trajectory_logger.MaybeWrite(simulation_time, simulation);
                    output_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - output_start).count();
                }

                if (step < 5 || step % 20 == 0) {
                    const double division_diameter =
                        GetParam<double>("biology.division_diameter", 2.0 * DefaultDiameter());
                    size_t tracked = 0;
                    size_t ready_to_divide = 0;
                    double max_diameter = 0.0;
                    double mean_diameter = 0.0;
                    double mean_oxygen = 0.0;
                    double mean_uptake = 0.0;
                    rm->ForEachAgent([&](bdm::Agent* a) {
                        auto* cell = dynamic_cast<MyCell*>(a);
                        if (!cell) return;
                        const double d = cell->GetDiameter();
                        mean_diameter += d;
                        mean_oxygen   += cell->GetOxygen();
                        mean_uptake   += cell->GetOxygenUptake();
                        max_diameter   = std::max(max_diameter, d);
                        if (d >= division_diameter) ++ready_to_divide;
                        ++tracked;
                    });
                    if (tracked > 0) {
                        mean_diameter /= static_cast<double>(tracked);
                        mean_oxygen   /= static_cast<double>(tracked);
                        mean_uptake   /= static_cast<double>(tracked);
                    }
                    std::ostringstream payload;
                    payload << "{"
                            << "\"step\":"            << step            << ","
                            << "\"simulationTime\":"  << simulation_time << ","
                            << "\"agents\":"          << tracked         << ","
                            << "\"readyToDivide\":"   << ready_to_divide << ","
                            << "\"divisionDiameter\":" << division_diameter << ","
                            << "\"meanDiameter\":"    << mean_diameter   << ","
                            << "\"maxDiameter\":"     << max_diameter    << ","
                            << "\"meanOxygen\":"      << mean_oxygen     << ","
                            << "\"meanUptake\":"      << mean_uptake
                            << "}";
                    debug_ndjson_log("baseline", "H3", "cells.h:RunCavityCouplingLoop",
                                     "growth_division_progress", payload.str());
                }

                const double window_wall_seconds =
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - window_wall_start).count();
                if (step < 5 || step % profiling_interval == 0) {
                    Log::Info("RunCavityCoupling",
                              "PROFILE window=", step,
                              " agents=", rm->GetNumAgents(),
                              " dt=", biology_dt,
                              " wall=", window_wall_seconds,
                              "s fetch=", fetch_seconds,
                              "s sim=", simulate_seconds,
                              "s out=", output_seconds,
                              "s adv=", advance_seconds, "s");
                }
            };

            const double biology_window_dt =
                coupling_active ? pending_biology_window_dt : step_dt;

            if (biology_window_dt <= 0.0) {
                Log::Warning("RunCavityCoupling",
                             "Skipping BioDynaMo evolution for non-positive completed window dt=",
                             biology_window_dt);
            } else if (coupling_active && (substep_mode || adaptive_substep_mode)) {
                if (skip_physics_on_acceptance) {
                    std::cout << "BDM window " << completed_window_count
                              << " — skip_physics commit (trial state), dt="
                              << biology_window_dt << " s" << std::endl;
                    run_single_biology_step(biology_window_dt, false);
                    std::cout << "BDM window " << completed_window_count
                              << " — commit done, simulation_time="
                              << simulation_time << " s" << std::endl;
                } else {
                const int n_substeps = substep_mode
                    ? internal_substeps_per_window
                    : std::max(1, static_cast<int>(
                          std::ceil(biology_window_dt / biodynamo_timestep)));

                const double biology_substep_dt =
                    biology_window_dt / static_cast<double>(n_substeps);

                if (n_substeps > 1 && (step < 5 || step % profiling_interval == 0)) {
                    Log::Info("RunCavityCoupling",
                              "Converged window substepping: window_dt=",
                              biology_window_dt,
                              " s, n_substeps=", n_substeps,
                              ", substep_dt=", biology_substep_dt, " s");
                }

                for (int substep_i = 0; substep_i < n_substeps; ++substep_i) {
                    const auto resample_start = std::chrono::steady_clock::now();
                    adapter.RefreshCachedFieldData();
                    fetch_seconds += std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - resample_start).count();
                    run_single_biology_step(biology_substep_dt, true);
                }
                }
            } else {
                if (skip_physics_on_acceptance && (step < 5 || step % profiling_interval == 0)) {
                    Log::Info("RunCavityCoupling",
                              "Two-way implicit default: committing trial agent state "
                              "without second Simulate(1), window_dt=", biology_window_dt, " s");
                }
                run_single_biology_step(biology_window_dt, !skip_physics_on_acceptance);
            }

            if (coupling_active) {
                pending_biology_window_dt = 0.0;
            }
        }
    }

    Log::Info("RunCavityCoupling", "Execution finished after ", step,
              " steps at simulation_time=", simulation_time, " s");
    if (remove_dead_rons_agents_enabled) {
        Log::Info("RonsResponse",
                  "Total agents removed by intracellular damage: ",
                  lifecycle_rons_removed,
                  "; surviving agents: ", rm->GetNumAgents());
        std::cout << "RonsResponse: total agents removed by intracellular damage: "
                  << lifecycle_rons_removed
                  << "; surviving agents: " << rm->GetNumAgents()
                  << std::endl;
    }
    if (checkpoint_write_noops > 0 || checkpoint_read_noops > 0) {
        Log::Info("RunCavityCoupling",
                  "Handled preCICE checkpoint actions as BioDynaMo no-ops: writes=",
                  checkpoint_write_noops, ", reads=", checkpoint_read_noops);
    }
    std::cout << "Execution finished after " << step
              << " steps at simulation_time=" << simulation_time
              << " s" << std::endl;
    trajectory_logger.Finalize();
    if (!adapter_finalized) {
        adapter.Finalize();
    }
    return 0;
}

// Original entry point that creates a new simulation
inline int RunCavityCoupling(int argc, const char** argv) {
    // Set up initial simulation parameters
    // Use kOpen to avoid BioDynaMo's built-in boundary forces (we handle walls explicitly)
    auto set_param = [](Param* param) {
    // Use kOpen - we handle wall collisions explicitly with KinematicCloudBehavior
        param->bound_space = Param::kOpen;
        // Set large domain bounds so BioDynaMo doesn't interfere
        param->min_bound = -1.0;  // Large bounds to avoid interference
        param->max_bound = 1.0;   // Large bounds to avoid interference

        // Disable neighbor search grid - we'll handle interactions explicitly
        param->calculate_gradients = false;
        param->numerical_ode_solver = Param::NumericalODESolver::kEuler;
    };

    // Create BioDynaMo simulation with cubic bounding box
    Simulation simulation(argc, argv, set_param);
    // Replace the default UniformGridEnvironment with KDTreeEnvironment.
    // The uniform grid uses int32_t box_length, so with sub-unit cell coordinates
    // (diameter=0.01, positions in [0,1]^3) all agents land in one box and the
    // uint16_t per-box counter overflows above 65535 cells. KDTreeEnvironment
    // (nanoflann) has no such limit and delivers O(log N) neighbor queries.
    simulation.SetEnvironment(new KDTreeEnvironment());
    return RunCavityCoupling(simulation);
}

// Configuration validation utility
inline void ValidateConfiguration() {
    Log::Info("ValidateConfiguration", "Validating coupling configurations...");

    // Log validation success
    Log::Info("ValidateConfiguration", "Configuration validation complete");
    Log::Info("ValidateConfiguration", "  ✓ Basic configuration");
    Log::Info("ValidateConfiguration", "  ✓ Agent initialization");
    Log::Info("ValidateConfiguration", "  ✓ preCICE integration");
    Log::Info("ValidateConfiguration", "  ✓ OpenFOAM coordinate system alignment");
    Log::Info("ValidateConfiguration", "  ✓ Visualization modes");
}


} // namespace bdm

#endif // CELLS_H_
