#ifndef BDM_ADAPTER_RUNTIME_CONFIG_H_
#define BDM_ADAPTER_RUNTIME_CONFIG_H_

/**
 * @file runtime_config.h
 * @brief Runtime configuration for bdm-adapter coupling
 *
 * This is the SINGLE SOURCE OF TRUTH for all coupling configuration.
 * Users configure via params.csv, and this module provides runtime access.
 *
 * Usage in cells participant:
 *   #include "bdm_adapter/runtime_config.h"
 *   const auto& config = bdm_adapter::GetRuntimeConfig();
 *   if (config.coupling.two_way_coupling) { ... }
 */

#include "parameters.h"
#include <functional>
#include <string>

namespace bdm_adapter {

/**
 * @brief Singleton runtime configuration
 *
 * Initialized once from params.csv, then accessible globally.
 * All coupling components should read from here instead of
 * defining their own configuration structs.
 */
class RuntimeConfig {
public:
    // Configuration sections (direct access - no duplication needed elsewhere)
    AdapterParameters::CouplingConfig coupling;
    AdapterParameters::FluidConfig fluid;
    AdapterParameters::ParticlePhysicsConfig particle_physics;
    AdapterParameters::AgentConfig agent;
    AdapterParameters::PreciceConfig precice;
    AdapterParameters::MeshConfig mesh;

    // Domain bounds - DEPRECATED: DO NOT USE
    // All bounds must come from preCICE mesh data via GetActualDomainBounds()
    // This struct is kept only for backwards compatibility but should not be used
    struct DomainBounds {
        double x_min = std::numeric_limits<double>::quiet_NaN();
        double x_max = std::numeric_limits<double>::quiet_NaN();
        double y_min = std::numeric_limits<double>::quiet_NaN();
        double y_max = std::numeric_limits<double>::quiet_NaN();
        double z_min = std::numeric_limits<double>::quiet_NaN();
        double z_max = std::numeric_limits<double>::quiet_NaN();
        double z_plane = std::numeric_limits<double>::quiet_NaN();
    } domain;

    // Convenience accessors (avoid verbose config.coupling.xxx)
    bool IsTwoWayCoupling() const { return coupling.two_way_coupling; }
    bool IsDragEnabled() const { return coupling.enable_drag_physics; }
    bool IsMasslessMode() const { return coupling.use_massless_tracers; }
    bool IsForce2D() const { return coupling.force_2d; }
    bool IsWallInteractionsEnabled() const { return coupling.enable_wall_interactions; }
    bool IsWallConfinementEnabled() const {
        return coupling.wall_confinement == "enforce";
    }
    const std::string& WallVelocityResponse() const {
        return coupling.wall_velocity_response;
    }
    bool IsParticleParticleEnabled() const { return coupling.enable_particle_particle; }

    // Fluid property accessors
    double FluidDensity() const { return fluid.density; }
    double FluidViscosity() const { return fluid.dynamic_viscosity; }
    double KinematicViscosity() const { return fluid.kinematic_viscosity(); }
};

// Forward declaration of loader (implemented below)
template<typename ParameterStore>
void LoadRuntimeConfig(RuntimeConfig& config, const ParameterStore& store);

/**
 * @brief Get the global runtime configuration (read-only)
 *
 * This is the primary access point for all coupling components.
 */
inline RuntimeConfig& GetMutableRuntimeConfig() {
    static RuntimeConfig instance;
    return instance;
}

inline const RuntimeConfig& GetRuntimeConfig() {
    return GetMutableRuntimeConfig();
}

/**
 * @brief Initialize runtime config from a parameter store
 *
 * Call this once at startup after loading params.csv.
 * Example:
 *   LoadParametersFromCsv("params.csv");
 *   bdm_adapter::InitializeRuntimeConfig(GetParameterStore());
 */
template<typename ParameterStore>
void InitializeRuntimeConfig(const ParameterStore& store) {
    LoadRuntimeConfig(GetMutableRuntimeConfig(), store);
}

/**
 * @brief Load configuration from parameter store
 */
template<typename ParameterStore>
void LoadRuntimeConfig(RuntimeConfig& config, const ParameterStore& store) {
    // Helper lambdas
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
    auto get_string = [&](const std::string& key, std::string& target) {
        if (store.template Has<std::string>(key)) {
            target = store.template Get<std::string>(key);
        }
    };

    // === Coupling Config ===
    get_bool("coupling.two_way", config.coupling.two_way_coupling);
    get_bool("coupling.enable_velocity", config.coupling.enable_velocity);
    get_bool("coupling.enable_pressure", config.coupling.enable_pressure);
    get_bool("coupling.enable_temperature", config.coupling.enable_temperature);
    get_bool("coupling.enable_drag_physics", config.coupling.enable_drag_physics);
    get_bool("coupling.enable_pressure_gradient", config.coupling.enable_pressure_gradient);
    get_bool("coupling.enable_oxygen", config.coupling.enable_oxygen);
    get_bool("coupling.enable_h2o2", config.coupling.enable_h2o2);
    get_bool("coupling.enable_no2", config.coupling.enable_no2);
    get_bool("coupling.enable_oxygen_sink", config.coupling.enable_oxygen_sink);
    get_bool("coupling.enable_h2o2_sink", config.coupling.enable_h2o2_sink);
    get_bool("coupling.enable_no2_sink", config.coupling.enable_no2_sink);
    get_double("coupling.oxygen_uptake_rate", config.coupling.oxygen_uptake_rate);
    get_double("coupling.oxygen_uptake_max", config.coupling.oxygen_uptake_max);
    get_double("coupling.oxygen_uptake_vmax", config.coupling.oxygen_uptake_vmax);
    get_double("coupling.oxygen_uptake_km", config.coupling.oxygen_uptake_km);
    get_double("coupling.h2o2_uptake_rate", config.coupling.h2o2_uptake_rate);
    get_double("coupling.h2o2_uptake_max", config.coupling.h2o2_uptake_max);
    get_double("coupling.no2_uptake_rate", config.coupling.no2_uptake_rate);
    get_double("coupling.no2_uptake_max", config.coupling.no2_uptake_max);
    get_bool("simulation.enable_wall_interactions", config.coupling.enable_wall_interactions);
    get_string("simulation.wall_confinement", config.coupling.wall_confinement);
    get_string("interaction.wall_velocity_response",
               config.coupling.wall_velocity_response);
    get_bool("simulation.enable_particle_particle", config.coupling.enable_particle_particle);
    get_bool("simulation.use_massless_tracers", config.coupling.use_massless_tracers);
    get_bool("simulation.force_2d", config.coupling.force_2d);

    // === Particle Physics Config ===
    get_double("interaction.youngs_modulus", config.particle_physics.youngs_modulus);
    get_double("interaction.poisson_ratio", config.particle_physics.poisson_ratio);
    get_double("interaction.restitution_coefficient", config.particle_physics.restitution_coefficient);
    get_double("interaction.friction_coefficient", config.particle_physics.friction_coefficient);
    get_double("interaction.wall_youngs_modulus", config.particle_physics.wall_youngs_modulus);
    get_double("interaction.wall_poisson_ratio", config.particle_physics.wall_poisson_ratio);
    get_double("interaction.wall_restitution", config.particle_physics.wall_restitution);

    // === Agent Config ===
    get_double("agent.default_diameter", config.agent.default_diameter);
    get_double("agent.default_mass", config.agent.default_mass);
    get_double("agent.default_plane_z", config.agent.default_plane_z);
    get_string("agent.custom_positions_csv", config.agent.custom_positions_csv);
    get_bool("agent.custom_positions_required", config.agent.custom_positions_required);

    // === Domain Config ===
    get_double("domain.x_min", config.domain.x_min);
    get_double("domain.x_max", config.domain.x_max);
    get_double("domain.y_min", config.domain.y_min);
    get_double("domain.y_max", config.domain.y_max);
    get_double("domain.z_min", config.domain.z_min);
    get_double("domain.z_max", config.domain.z_max);
    get_double("domain.z_plane", config.domain.z_plane);

    // === preCICE Config ===
    get_string("precice.participant_name", config.precice.participant_name);
    get_string("precice.config_file", config.precice.config_file);
    get_bool("precice.enable_profiling", config.precice.enable_profiling);
    get_string("precice.velocity_data_name", config.precice.velocity_data_name);
    get_string("precice.pressure_data_name", config.precice.pressure_data_name);
    get_string("precice.pressure_gradient_data_name", config.precice.pressure_gradient_data_name);
    get_string("precice.fluid_acceleration_data_name", config.precice.fluid_acceleration_data_name);
    get_string("precice.fluid_vorticity_data_name", config.precice.fluid_vorticity_data_name);
    get_string("precice.temperature_data_name", config.precice.temperature_data_name);
    get_string("precice.position_data_name", config.precice.position_data_name);
    get_string("precice.wall_geometry_data_name", config.precice.wall_geometry_data_name);
    get_string("precice.wall_normals_data_name", config.precice.wall_normals_data_name);
    get_string("precice.wall_areas_data_name", config.precice.wall_areas_data_name);
    get_string("precice.wall_velocity_data_name", config.precice.wall_velocity_data_name);
    get_string("precice.inlet_cell_centers_data_name", config.precice.inlet_cell_centers_data_name);
    get_string("precice.inlet_areas_data_name", config.precice.inlet_areas_data_name);
    get_string("precice.fluid_density_data_name", config.precice.fluid_density_data_name);
    get_string("precice.fluid_dynamic_viscosity_data_name", config.precice.fluid_dynamic_viscosity_data_name);
    get_string("precice.oxygen_data_name", config.precice.oxygen_data_name);
    get_string("precice.oxygen_sink_data_name", config.precice.oxygen_sink_data_name);
    get_string("precice.h2o2_data_name", config.precice.h2o2_data_name);
    get_string("precice.no2_data_name", config.precice.no2_data_name);
    get_string("precice.h2o2_sink_data_name", config.precice.h2o2_sink_data_name);
    get_string("precice.no2_sink_data_name", config.precice.no2_sink_data_name);
    get_string("precice.tke_data_name", config.precice.tke_data_name);
    get_string("precice.epsilon_data_name", config.precice.epsilon_data_name);
    get_string("precice.omega_data_name", config.precice.omega_data_name);
    get_string("precice.nut_data_name", config.precice.nut_data_name);

    // === Coupling Config ===
    get_bool("coupling.enable_turbulence", config.coupling.enable_turbulence);

    // === Mesh Config ===
    get_string("mesh.particle_mesh", config.mesh.particle_mesh);
    get_string("mesh.primary_volume_mesh", config.mesh.primary_volume_mesh);
    get_string("mesh.foam_case_directory", config.mesh.foam_case_directory);
}

} // namespace bdm_adapter

#endif // BDM_ADAPTER_RUNTIME_CONFIG_H_
