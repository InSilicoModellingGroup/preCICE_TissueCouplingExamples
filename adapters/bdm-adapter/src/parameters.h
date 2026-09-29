#ifndef BDM_ADAPTER_PARAMETERS_H_
#define BDM_ADAPTER_PARAMETERS_H_

#include <limits>
#include <string>
#include <vector>

namespace bdm_adapter {

/**
 * @brief Configuration parameters for the preCICE adapter
 *
 * This structure holds all configuration needed to set up the coupling
 * between BioDynaMo and preCICE. Users can customize these parameters
 * or use the provided defaults.
 */
struct AdapterParameters {
    /**
     * @brief preCICE-specific configuration
     */
    struct PreciceConfig {
        std::string participant_name = "cells";
        std::string config_file = "../precice-config.xml";
        bool enable_profiling = true;
        std::string profiling_variable = "PRECICE_PROFILING";
        std::string profiling_value = "1";
        std::string profiling_mode_variable = "PRECICE_PROFILING_MODE";
        std::string profiling_mode_value = "all";

        // Data field names (must match preCICE config)
        std::string velocity_data_name = "U";
        std::string pressure_data_name = "p";
        std::string pressure_gradient_data_name = "gradP";
        std::string fluid_acceleration_data_name = "DUcDt";
        std::string fluid_vorticity_data_name = "curlUc";
        std::string temperature_data_name = "T";
        std::string position_data_name = "Position";
        std::string wall_geometry_data_name = "WallGeometry";
        std::string wall_normals_data_name = "WallNormals";
        std::string wall_areas_data_name = "WallAreas";
        std::string wall_velocity_data_name = "WallVelocity";  // For moving walls
        std::string boundary_vertex0_data_name = "BoundaryTriangleVertex0";
        std::string boundary_vertex1_data_name = "BoundaryTriangleVertex1";
        std::string boundary_vertex2_data_name = "BoundaryTriangleVertex2";
        std::string boundary_kind_data_name = "BoundaryKind";
        std::string boundary_normals_data_name = "BoundaryNormals";
        std::string boundary_velocity_data_name = "BoundaryVelocity";
        std::string inlet_cell_centers_data_name = "InletCellCenters";
        std::string inlet_areas_data_name = "InletAreas";
        std::string fluid_density_data_name = "rho";
        std::string fluid_dynamic_viscosity_data_name = "mu";
        std::string oxygen_data_name = "O2";          // Oxygen concentration from OpenFOAM
        std::string oxygen_sink_data_name = "O2Sink"; // Oxygen sink feedback to OpenFOAM
        std::string h2o2_data_name = "H2O2";
        std::string no2_data_name = "NO2";
        std::string h2o2_sink_data_name = "H2O2Sink";
        std::string no2_sink_data_name = "NO2Sink";   // Nitrite sink feedback to OpenFOAM
        std::string tke_data_name = "k";              // Turbulent kinetic energy from OpenFOAM
        std::string epsilon_data_name = "epsilon";    // Turbulent dissipation rate from OpenFOAM
        std::string omega_data_name = "omega";        // Specific dissipation rate from OpenFOAM
        std::string nut_data_name = "nut";            // Turbulent viscosity from OpenFOAM
    } precice;

    /**
     * @brief Mesh configuration
     */
    struct MeshConfig {
        std::string particle_mesh = "ABMmesh";
        std::string primary_volume_mesh = "FVMmesh";
        std::vector<std::string> default_volume_meshes = {"FVMmesh"};
        std::vector<std::string> default_wall_meshes = {"WallSurfaceMesh"};
        std::string boundary_triangle_mesh = "BoundaryTriangleMesh";
        std::string foam_case_directory = "../cavity2D";
        std::string block_mesh_dict = "system/blockMeshDict";

        // Inlet/Outlet surface meshes for agent generation and removal
        std::string inlet_mesh = "InletSurfaceMesh";   // Mesh for inlet face centers
        std::string outlet_mesh = "OutletSurfaceMesh"; // Mesh for outlet face centers
        std::vector<std::string> default_inlet_meshes = {};  // Empty by default (optional feature)
        std::vector<std::string> default_outlet_meshes = {}; // Empty by default (optional feature)
    } mesh;

    /**
     * @brief Agent (particle) configuration
     */
    struct AgentConfig {
        double default_diameter = 0.004;
        double corner_diameter = 0.008;
        double bounding_box_margin = 0.0;
        double default_plane_z = 0.005;
        double default_mass = 1.814882032e-22;
        std::string custom_positions_csv = "agents_positions.csv";
        bool custom_positions_required = false;

        // Inlet-based agent generation (alternative to CSV or full mesh)
        bool inlet_generation = false;           // Create agents at inlet mesh vertices
        bool inlet_random_sampling = true;       // true: random capped subset; false: repeatable mapped-cell subset
        double inlet_spawn_interval = 0.0;       // Time between spawns (0 = one-time at init)
        int inlet_agents_per_spawn = 1;          // Agents to create per spawn event
        double inlet_zone_offset = 0.01;         // Offset from inlet face into domain [m]
        double inlet_max_cell_distance = 0.0;    // Max inlet face-to-cell mapping distance [m] (0 = auto)
        double inlet_occupied_tolerance = 0.0;   // Skip inlet cell center if an agent is already this close [m] (0 = auto)
        bool inlet_filter_occupied = true;       // Enforce existing-agent and within-batch spacing (disable for one agent per selected cell)
        double inlet_center_radius_fraction = 1.0; // Keep centers within this fraction of inlet radius (1 = full inlet)
        bool inlet_keep_largest_component = false; // Keep only the largest connected inlet face cluster
        double inlet_component_radius = 0.0;     // Face-center clustering radius [m] (0 = auto)
        std::string inlet_midplane_axis = "none"; // Legacy single inlet slice normal: x, y, z, or none
        std::vector<std::string> inlet_midplane_axes = {}; // Optional union of inlet slice normals (e.g. y;z for a cross)
        double inlet_midplane_tolerance = 0.0;   // Slice half-width [m] (0 = infer from candidate spacing)

        // Outlet-based agent removal
        bool outlet_removal = false;             // Remove agents near outlet
        double outlet_removal_distance = 0.02;   // Distance threshold for removal [m]

        // Wall clearance filtering for CSV-based agent initialization
        // When enabled, positions closer than wall_clearance_distance to any
        // wall-mesh face centre are excluded during CreateAgentsFromCSV.
        bool enable_wall_clearance = false;      // Activate the filter
        double wall_clearance_distance = 0.0;    // Minimum distance from wall [m]
        int wall_clearance_layers = 0;           // Optional layer count to exclude near wall
        double wall_clearance_layer_thickness = 0.0; // Per-layer thickness [m]
    } agent;

    /**
     * @brief Environment variable configuration
     */
    struct EnvironmentConfig {
        std::string volume_meshes_variable = "BDM_VOLUME_MESHES";
        std::string primary_volume_variable = "BDM_PRIMARY_VOLUME_MESH";
        std::string wall_meshes_variable = "BDM_WALL_MESHES";
    } environment;

    /**
     * @brief Coupling mode configuration (one-way vs two-way)
     */
    struct CouplingConfig {
        // Logical flow regime descriptor. This is independent from sampling
        // policy and is intended for high-level behavior toggles in callers.
        // Supported values: "transient" | "steady"
        std::string flow_mode = "transient";

        // Field re-sampling cadence for mapAndReadData/readData calls.
        //   1  -> read every FetchInputData call (default, transient)
        //   0  -> read once then freeze values on agents
        //   N>1 -> read every N FetchInputData calls
        int field_read_interval = 1;

        // Relative read time within the current preCICE coupling window [s],
        // passed to mapAndReadData/readData. Supported modes:
        //   window_start -> 0 (start of window; last committed producer state)
        //   window_end   -> current window size (trial end-of-window field)
        // A non-negative numeric literal is also accepted and used directly.
        std::string field_read_relative_time = "window_start";

        // Wall/geometry mesh policy.
        // If true, wall geometry is read exactly once and then reused.
        bool static_mesh = false;

        // Wall data re-sampling cadence when static_mesh=false.
        //   1  -> read every FetchInputData call
        //   0  -> read once
        //   N>1 -> read every N FetchInputData calls
        int wall_read_interval = 1;

        bool two_way_coupling = false;           // Send positions back to CFD solver
        bool enable_velocity = true;             // Read velocity from CFD
        bool enable_pressure = true;             // Read pressure from CFD
        bool enable_temperature = false;         // Read temperature from CFD
        bool enable_drag_physics = true;         // Apply Schiller-Naumann drag
        bool enable_pressure_gradient = true;    // Apply pressure gradient force
        bool enable_wall_interactions = true;    // Apply DEM wall collisions
        std::string wall_confinement = "enforce";
        std::string wall_velocity_response = "rebound";
        bool enable_particle_particle = false;   // Apply DEM particle-particle
        bool use_massless_tracers = false;       // Particles follow fluid exactly
        bool force_2d = true;                    // Constrain motion to 2D plane
        bool enable_oxygen = false;              // Read oxygen scalar field from CFD
        bool enable_h2o2 = false;                // Read H2O2 directly from CFD
        bool enable_no2 = false;                 // Read NO2- directly from CFD
        bool enable_oxygen_sink = false;         // Write oxygen sink feedback to CFD
        bool enable_h2o2_sink = false;           // Write H2O2Sink feedback to CFD
        bool enable_no2_sink = false;            // Write NO2- sink feedback to CFD
        double oxygen_uptake_rate = 0.0;         // Linear uptake coeff [1/s]: uptake = k * O2 (used when Km = 0)
        double oxygen_uptake_max = 1.0e30;       // Optional cap on per-agent uptake [conc/s]
        double oxygen_uptake_vmax = 0.0;         // MM Vmax [conc/s]: max uptake at saturation (0 = disabled)
        double oxygen_uptake_km = 0.0;           // MM Km [conc]: half-saturation constant (0 = use linear model)
        double h2o2_uptake_rate = 0.0;           // First-order H2O2 uptake coefficient [1/s]
        double h2o2_uptake_max = 1.0e30;         // Optional cap on per-agent H2O2 uptake [conc/s]
        double no2_uptake_rate = 0.0;            // First-order NO2- uptake coefficient [1/s]
        double no2_uptake_max = 1.0e30;          // Optional cap on per-agent NO2- uptake [conc/s]
        bool enable_turbulence = false;          // Read turbulence fields (k, epsilon/omega, nut) from CFD
        bool distribute_feedback = false;        // Spread finite-size agent feedback over nearby FVM cells
        double feedback_kernel_radius_factor = 0.5; // Kernel radius = factor * agent diameter
        int feedback_kernel_min_cells = 1;       // Minimum cells used by feedback kernel fallback
    } coupling;

    /**
     * @brief Fluid properties configuration
     */
    struct FluidConfig {
        double density = std::numeric_limits<double>::quiet_NaN();          // Set by preCICE
        double dynamic_viscosity = std::numeric_limits<double>::quiet_NaN(); // Set by preCICE
        double kinematic_viscosity() const { return dynamic_viscosity / density; }
    } fluid;

    /**
     * @brief Particle/agent physics properties
     */
    struct ParticlePhysicsConfig {
        double youngs_modulus = 1.0e-5;          // [Pa] - Particle stiffness
        double poisson_ratio = 0.30;             // [-] - Poisson ratio
        double density = 1.2;                    // [kg/m^3] - Particle density
        double restitution_coefficient = 0.9;    // [-] - Rebound coefficient
        double friction_coefficient = 0.3;       // [-] - Friction coefficient
        double wall_youngs_modulus = 1.0e-3;     // [Pa] - Wall stiffness
        double wall_poisson_ratio = 0.30;        // [-] - Wall Poisson ratio
        double wall_restitution = 0.8;           // [-] - Wall restitution
    } particle_physics;
};

/**
 * @brief Get default adapter parameters
 * @return Default AdapterParameters instance
 */
inline AdapterParameters GetDefaultAdapterParameters() {
    return AdapterParameters{};
}

/**
 * @brief Build AdapterParameters from a ParameterStore (loaded from CSV)
 *
 * This function reads parameters from the store and populates the
 * AdapterParameters struct. Users only need to modify params.csv.
 *
 * @tparam ParameterStore Type providing Has<T> and Get<T> methods
 * @param store The parameter store to read from
 * @return Configured AdapterParameters
 */
template<typename ParameterStore>
AdapterParameters BuildAdapterParametersFromStore(const ParameterStore& store) {
    AdapterParameters params;

    // Helper lambdas for type-safe parameter extraction
    auto get_string = [&](const std::string& key, std::string& target) {
        if (store.template Has<std::string>(key)) {
            target = store.template Get<std::string>(key);
        }
    };
    auto get_bool = [&](const std::string& key, bool& target) {
        if (store.template Has<bool>(key)) {
            target = store.template Get<bool>(key);
        }
    };
    auto get_double = [&](const std::string& key, double& target) {
        if (store.template Has<double>(key)) {
            target = store.template Get<double>(key);
        }
    };
    auto get_string_list = [&](const std::string& key, std::vector<std::string>& target) {
        if (store.template Has<std::vector<std::string>>(key)) {
            target = store.template Get<std::vector<std::string>>(key);
        }
    };

    // preCICE config
    get_string("precice.participant_name", params.precice.participant_name);
    get_string("precice.config_file", params.precice.config_file);
    get_bool("precice.enable_profiling", params.precice.enable_profiling);
    get_string("precice.profiling_variable", params.precice.profiling_variable);
    get_string("precice.profiling_value", params.precice.profiling_value);
    get_string("precice.profiling_mode_variable", params.precice.profiling_mode_variable);
    get_string("precice.profiling_mode_value", params.precice.profiling_mode_value);
    get_string("precice.velocity_data_name", params.precice.velocity_data_name);
    get_string("precice.pressure_data_name", params.precice.pressure_data_name);
    get_string("precice.pressure_gradient_data_name", params.precice.pressure_gradient_data_name);
    get_string("precice.fluid_acceleration_data_name", params.precice.fluid_acceleration_data_name);
    get_string("precice.fluid_vorticity_data_name", params.precice.fluid_vorticity_data_name);
    get_string("precice.temperature_data_name", params.precice.temperature_data_name);
    get_string("precice.position_data_name", params.precice.position_data_name);
    get_string("precice.wall_geometry_data_name", params.precice.wall_geometry_data_name);
    get_string("precice.wall_normals_data_name", params.precice.wall_normals_data_name);
    get_string("precice.wall_areas_data_name", params.precice.wall_areas_data_name);
    get_string("precice.wall_velocity_data_name", params.precice.wall_velocity_data_name);
    get_string("precice.boundary_vertex0_data_name", params.precice.boundary_vertex0_data_name);
    get_string("precice.boundary_vertex1_data_name", params.precice.boundary_vertex1_data_name);
    get_string("precice.boundary_vertex2_data_name", params.precice.boundary_vertex2_data_name);
    get_string("precice.boundary_kind_data_name", params.precice.boundary_kind_data_name);
    get_string("precice.boundary_normals_data_name", params.precice.boundary_normals_data_name);
    get_string("precice.boundary_velocity_data_name", params.precice.boundary_velocity_data_name);
    get_string("precice.inlet_cell_centers_data_name", params.precice.inlet_cell_centers_data_name);
    get_string("precice.inlet_areas_data_name", params.precice.inlet_areas_data_name);
    get_string("precice.fluid_density_data_name", params.precice.fluid_density_data_name);
    get_string("precice.fluid_dynamic_viscosity_data_name", params.precice.fluid_dynamic_viscosity_data_name);
    get_string("precice.oxygen_data_name", params.precice.oxygen_data_name);
    get_string("precice.oxygen_sink_data_name", params.precice.oxygen_sink_data_name);
    get_string("precice.h2o2_data_name", params.precice.h2o2_data_name);
    get_string("precice.no2_data_name", params.precice.no2_data_name);
    get_string("precice.h2o2_sink_data_name", params.precice.h2o2_sink_data_name);
    get_string("precice.no2_sink_data_name", params.precice.no2_sink_data_name);
    get_string("precice.tke_data_name", params.precice.tke_data_name);
    get_string("precice.epsilon_data_name", params.precice.epsilon_data_name);
    get_string("precice.omega_data_name", params.precice.omega_data_name);
    get_string("precice.nut_data_name", params.precice.nut_data_name);

    // Mesh config
    get_string("mesh.particle_mesh", params.mesh.particle_mesh);
    get_string("mesh.primary_volume_mesh", params.mesh.primary_volume_mesh);
    get_string_list("mesh.default_volume_meshes", params.mesh.default_volume_meshes);
    get_string_list("mesh.default_wall_meshes", params.mesh.default_wall_meshes);
    get_string("mesh.boundary_triangle_mesh", params.mesh.boundary_triangle_mesh);
    get_string("mesh.foam_case_directory", params.mesh.foam_case_directory);
    get_string("mesh.inlet_mesh", params.mesh.inlet_mesh);
    get_string("mesh.outlet_mesh", params.mesh.outlet_mesh);
    get_string_list("mesh.default_inlet_meshes", params.mesh.default_inlet_meshes);
    get_string_list("mesh.default_outlet_meshes", params.mesh.default_outlet_meshes);

    // Agent config
    get_double("agent.default_diameter", params.agent.default_diameter);
    get_double("agent.default_mass", params.agent.default_mass);
    get_double("agent.default_plane_z", params.agent.default_plane_z);
    get_string("agent.custom_positions_csv", params.agent.custom_positions_csv);
    get_bool("agent.custom_positions_required", params.agent.custom_positions_required);
    get_bool("agent.inlet_generation", params.agent.inlet_generation);
    get_bool("agent.inlet_random_sampling", params.agent.inlet_random_sampling);
    get_double("agent.inlet_spawn_interval", params.agent.inlet_spawn_interval);
    auto get_int = [&](const std::string& key, int& target) {
        if (store.template Has<int>(key)) {
            target = store.template Get<int>(key);
        }
    };
    get_int("agent.inlet_agents_per_spawn", params.agent.inlet_agents_per_spawn);
    get_double("agent.inlet_zone_offset", params.agent.inlet_zone_offset);
    get_double("agent.inlet_max_cell_distance", params.agent.inlet_max_cell_distance);
    get_double("agent.inlet_occupied_tolerance", params.agent.inlet_occupied_tolerance);
    get_bool("agent.inlet_filter_occupied", params.agent.inlet_filter_occupied);
    get_double("agent.inlet_center_radius_fraction", params.agent.inlet_center_radius_fraction);
    get_bool("agent.inlet_keep_largest_component", params.agent.inlet_keep_largest_component);
    get_double("agent.inlet_component_radius", params.agent.inlet_component_radius);
    get_string("agent.inlet_midplane_axis", params.agent.inlet_midplane_axis);
    get_string_list("agent.inlet_midplane_axes", params.agent.inlet_midplane_axes);
    get_double("agent.inlet_midplane_tolerance", params.agent.inlet_midplane_tolerance);
    get_bool("agent.outlet_removal", params.agent.outlet_removal);
    get_double("agent.outlet_removal_distance", params.agent.outlet_removal_distance);
    get_bool("agent.enable_wall_clearance", params.agent.enable_wall_clearance);
    get_double("agent.wall_clearance_distance", params.agent.wall_clearance_distance);
    get_int("agent.wall_clearance_layers", params.agent.wall_clearance_layers);
    get_double("agent.wall_clearance_layer_thickness", params.agent.wall_clearance_layer_thickness);

    // Coupling config
    get_string("coupling.flow_mode", params.coupling.flow_mode);
    get_int("coupling.field_read_interval", params.coupling.field_read_interval);
    get_string("coupling.field_read_relative_time", params.coupling.field_read_relative_time);
    get_bool("coupling.static_mesh", params.coupling.static_mesh);
    get_int("coupling.wall_read_interval", params.coupling.wall_read_interval);

    get_bool("coupling.two_way", params.coupling.two_way_coupling);
    get_bool("coupling.enable_velocity", params.coupling.enable_velocity);
    get_bool("coupling.enable_pressure", params.coupling.enable_pressure);
    get_bool("coupling.enable_temperature", params.coupling.enable_temperature);
    get_bool("coupling.enable_drag_physics", params.coupling.enable_drag_physics);
    get_bool("coupling.enable_pressure_gradient", params.coupling.enable_pressure_gradient);
    get_bool("coupling.enable_oxygen", params.coupling.enable_oxygen);
    get_bool("coupling.enable_h2o2", params.coupling.enable_h2o2);
    get_bool("coupling.enable_no2", params.coupling.enable_no2);
    get_bool("coupling.enable_oxygen_sink", params.coupling.enable_oxygen_sink);
    get_bool("coupling.enable_h2o2_sink", params.coupling.enable_h2o2_sink);
    get_bool("coupling.enable_no2_sink", params.coupling.enable_no2_sink);
    get_double("coupling.oxygen_uptake_rate", params.coupling.oxygen_uptake_rate);
    get_double("coupling.oxygen_uptake_max", params.coupling.oxygen_uptake_max);
    get_double("coupling.oxygen_uptake_vmax", params.coupling.oxygen_uptake_vmax);
    get_double("coupling.oxygen_uptake_km", params.coupling.oxygen_uptake_km);
    get_double("coupling.h2o2_uptake_rate", params.coupling.h2o2_uptake_rate);
    get_double("coupling.h2o2_uptake_max", params.coupling.h2o2_uptake_max);
    get_double("coupling.no2_uptake_rate", params.coupling.no2_uptake_rate);
    get_double("coupling.no2_uptake_max", params.coupling.no2_uptake_max);
    get_bool("simulation.enable_wall_interactions", params.coupling.enable_wall_interactions);
    get_string("simulation.wall_confinement", params.coupling.wall_confinement);
    get_string("interaction.wall_velocity_response",
               params.coupling.wall_velocity_response);
    get_bool("simulation.enable_particle_particle", params.coupling.enable_particle_particle);
    get_bool("simulation.use_massless_tracers", params.coupling.use_massless_tracers);
    get_bool("simulation.force_2d", params.coupling.force_2d);
    get_bool("coupling.distribute_feedback", params.coupling.distribute_feedback);
    get_double("coupling.feedback_kernel_radius_factor", params.coupling.feedback_kernel_radius_factor);
    get_int("coupling.feedback_kernel_min_cells", params.coupling.feedback_kernel_min_cells);

    // Particle physics config
    get_double("interaction.youngs_modulus", params.particle_physics.youngs_modulus);
    get_double("interaction.poisson_ratio", params.particle_physics.poisson_ratio);
    get_double("interaction.restitution_coefficient", params.particle_physics.restitution_coefficient);
    get_double("interaction.friction_coefficient", params.particle_physics.friction_coefficient);
    get_double("interaction.wall_youngs_modulus", params.particle_physics.wall_youngs_modulus);
    get_double("interaction.wall_poisson_ratio", params.particle_physics.wall_poisson_ratio);
    get_double("interaction.wall_restitution", params.particle_physics.wall_restitution);

    return params;
}

} // namespace bdm_adapter

#endif // BDM_ADAPTER_PARAMETERS_H_
