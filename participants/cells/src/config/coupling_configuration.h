#ifndef COUPLING_CONFIG_H_
#define COUPLING_CONFIG_H_

#include <array>
#include <string>
#include <limits>
#include <cmath>

#include "parameter_loader.h"
#include "../fluid_properties.h"

namespace bdm {

/**
 * @brief Template-based coupling configuration with boolean feature flags
 *
 * This struct provides compile-time and runtime configuration for the
 * OpenFOAM-preCICE-BioDynaMo coupling with true/false parameters for
 * different coupling features as requested.
 */
template<bool EnableVelocity = true,
         bool EnablePressure = true,
         bool EnableTemperature = false,
         bool TwoWayCoupling = false,
         bool EnableParticleParticleInteractions = true,
         bool EnableParticleWallInteractions = true,
         bool EnableMotion = true,  // Start with stationary agents
         bool Force2D = true,
         bool EnableDragPhysics = true,  // Enable/disable drag force calculations
         bool UseMasslessTracers = false, // Use massless particles that follow fluid exactly
         bool EnablePressureGradientForce = true, // Enable pressure gradient force F = -V_p * ∇p
         bool EnableTurbulence = false>             // Enable turbulence field coupling (k, epsilon/omega, nut)
struct CouplingConfig {
    // Compile-time feature flags
    static constexpr bool enable_velocity = EnableVelocity;
    static constexpr bool enable_pressure = EnablePressure;
    static constexpr bool enable_temperature = EnableTemperature;
    static constexpr bool two_way_coupling = TwoWayCoupling;
    static constexpr bool enable_particle_particle = EnableParticleParticleInteractions;
    static constexpr bool enable_particle_wall = EnableParticleWallInteractions;
    static constexpr bool enable_motion = EnableMotion;
    static constexpr bool force_2d = Force2D;
    static constexpr bool enable_drag_physics = EnableDragPhysics;
    static constexpr bool use_massless_tracers = UseMasslessTracers;
    static constexpr bool enable_pressure_gradient_force = EnablePressureGradientForce;
    static constexpr bool enable_turbulence = EnableTurbulence;

    // Domain configuration - bounds come from preCICE mesh at runtime
    // Initialized to invalid values to force mesh-derived bounds
    struct DomainConfig {
        // Fluid domain bounds - set from FVMmesh via adapter.GetDomainBoundsFromMesh()
        // Invalid defaults ensure we don't accidentally use hardcoded values
        double x_min = std::numeric_limits<double>::quiet_NaN();
        double y_min = std::numeric_limits<double>::quiet_NaN();
        double z_min = std::numeric_limits<double>::quiet_NaN();
        double x_max = std::numeric_limits<double>::quiet_NaN();
        double y_max = std::numeric_limits<double>::quiet_NaN();
        double z_max = std::numeric_limits<double>::quiet_NaN();
        double z_plane = std::numeric_limits<double>::quiet_NaN();

        // Check if bounds are valid (set from mesh)
        bool IsValid() const {
            return !std::isnan(x_min) && !std::isnan(x_max) &&
                   !std::isnan(y_min) && !std::isnan(y_max) &&
                   !std::isnan(z_min) && !std::isnan(z_max);
        }

        // BioDynaMo bounds equal CFD bounds (no manual margin)
        double BdmXMin() const { return x_min; }
        double BdmYMin() const { return y_min; }
        double BdmZMin() const { return z_min; }
        double BdmXMax() const { return x_max; }
        double BdmYMax() const { return y_max; }
        double BdmZMax() const { return z_max; }
    } domain;

    // preCICE config - defaults are placeholders, actual values come from params.csv
    struct PreciceConfig {
        std::string config_file = "../precice-config.xml";
        std::string participant_name = "cells";
        std::string volume_mesh = "FVMmesh";
        std::string particle_mesh = "ABMmesh";
        std::string velocity_field = "U";
        std::string pressure_field = "p";
        std::string temperature_field = "T";
        std::string position_field = "Position";
    } precice;

    // Agent config - defaults are placeholders, actual values come from params.csv
    // via BuildAgentConfigFromParameters()
    struct AgentConfig {
        double initial_diameter = 0.004;
        double default_temperature = 300.0;
        double default_pressure = 0.0;
        bool place_at_cell_centers = true;
        bool initial_stationary = true;
    } agent;

    // Fluid config - defaults are placeholders, actual values come from params.csv
    // via BuildFluidConfigFromParameters()
    struct FluidConfig {
        double density = 0.1;
        double dynamic_viscosity = 1e-3;
        double kinematic_viscosity() const { return dynamic_viscosity / density; }
    } fluid;

    // Interaction config - defaults are placeholders, actual values come from params.csv
    // via BuildInteractionConfigFromParameters()
    struct InteractionConfig {
        double young_modulus = 1.0e-5;
        double poisson_ratio = 0.30;
        double density = 1.2;
        double friction_coefficient = 0.0;
        double restitution_coefficient = 0.1;
        double damping_factor = 0.02;
        double cohesion_energy = 0.0;
        double wall_young_modulus = 1.0e-5;
        double wall_poisson_ratio = 0.30;
        double interaction_range = 2.0e-3;
        double force_ramp_time = 0.0;
    } interaction;

    // Visualization config - defaults are placeholders, actual values come from params.csv
    // via BuildVisualizationConfigFromParameters()
    struct VisualizationConfig {
        double temperature_scale_min = 295.0;
        double temperature_scale_max = 310.0;
    } visualization;

    // Validation method to ensure configuration consistency
    bool validate() const {
        if constexpr (force_2d) {
            if (domain.z_max - domain.z_min < 1e-6) {
                return false; // Invalid Z dimension for 2D
            }
        }

        if constexpr (enable_temperature && !EnableTemperature) {
            return false; // Compile-time mismatch
        }

        return true;
    }
};

// Predefined configurations for common use cases
// Template parameters: EnableVelocity, EnablePressure, EnableTemperature, TwoWayCoupling,
//                      EnableParticleParticleInteractions, EnableParticleWallInteractions,
//                      EnableMotion, Force2D, EnableDragPhysics, UseMasslessTracers, EnablePressureGradientForce, EnableTurbulence
using OneWayVelocityPressure = CouplingConfig<true, true, false, false, true, true, false, true, true, false, true, false>;
using TwoWayVelocityPressure = CouplingConfig<true, true, false, true, true, true, false, true, true, false, true, false>;
using FullCouplingStationary = CouplingConfig<true, true, true, false, true, true, false, true, true, false, true, false>;
using FullCouplingDynamic = CouplingConfig<true, true, true, true, true, true, true, true, true, false, true, false>;

// Massless tracer configuration: particles move exactly with fluid velocity
// No drag forces, no collision forces, pure advection with the flow
using MasslessTracerConfig = CouplingConfig<true,   // EnableVelocity - receive from OpenFOAM
                                           true,   // EnablePressure - receive from OpenFOAM
                                           true,   // EnableTemperature - receive from OpenFOAM
                                           false,  // TwoWayCoupling - no force feedback to fluid
                                           false,  // EnableParticleParticleInteractions - disabled
                                           false,  // EnableParticleWallInteractions - disabled
                                           true,   // EnableMotion - particles move
                                           true,   // Force2D - keep 2D motion
                                           false,  // EnableDragPhysics - no drag forces
                                           true,   // UseMasslessTracers - massless mode
                                           false,  // EnablePressureGradientForce - disabled for massless
                                           false>; // EnableTurbulence - disabled for massless

// Finite-mass passive tracer configuration: retains drag physics and wall interactions
using PassiveTracerConfig = CouplingConfig<true,   // EnableVelocity - receive from OpenFOAM
                                          true,   // EnablePressure - receive from OpenFOAM
                                          false,  // EnableTemperature - DISABLED (not in cavity preCICE config)
                                          false,  // TwoWayCoupling - no force feedback to fluid
                                          false,  // EnableParticleParticleInteractions - disabled
                                          true,   // EnableParticleWallInteractions - enabled (match OpenFOAM walls)
                                          true,   // EnableMotion - particles move
                                          true,   // Force2D - keep 2D motion
                                          true,   // EnableDragPhysics - apply drag forces
                                          false,  // UseMasslessTracers - finite-mass tracers
                                          true,   // EnablePressureGradientForce - enabled for finite-mass
                                          false>; // EnableTurbulence - disabled by default

// Default configuration for cavity example - finite-mass passive tracers with drag
using CavityConfig = PassiveTracerConfig;

/**
 * @brief Runtime coupling configuration loaded from params.csv
 *
 * This struct provides runtime access to coupling parameters without requiring
 * recompilation. Users can change params.csv to switch between coupling modes.
 */
struct RuntimeCouplingConfig {
    bool two_way_coupling = false;
    bool enable_velocity = true;
    bool enable_pressure = true;
    bool enable_temperature = false;
    bool enable_oxygen = false;
    bool enable_h2o2 = false;
    bool enable_no2 = false;
    bool enable_oxygen_sink = false;
    bool enable_h2o2_sink = false;
    bool enable_no2_sink = false;
    bool enable_drag_physics = true;
    bool enable_pressure_gradient = false;
    bool enable_wall_interactions = true;
    std::string wall_confinement = "enforce";
    std::string wall_velocity_response = "rebound";
    bool enable_particle_particle = false;
    bool use_massless_tracers = false;
    bool force_2d = true;
    bool enable_particle_motion = true;

    // Additional physics forces
    bool enable_gravity = true;
    bool enable_virtual_mass = false;
    bool enable_lift_forces = false;
    bool enable_turbulence_dispersion = false;
    bool enable_turbulence = false;            // Read k, epsilon/omega, nut from CFD

    // Environment parameters
    double gravity_x = 0.0;
    double gravity_y = -9.81;
    double gravity_z = 0.0;

    // Particle properties
    double particle_density = 1000.0;
    double oxygen_uptake_rate = 0.0;
    double oxygen_uptake_max = 1.0e30;

    // Inlet/Outlet agent management (for complex geometries)
    bool inlet_generation = false;        // Create agents at inlet mesh vertices
    double inlet_spawn_interval = 0.0;    // Time between spawns (0 = one-time at init)
    int inlet_agents_per_spawn = 1;       // Agents per spawn event
    double inlet_zone_offset = 0.01;      // Offset from inlet face into domain [m]
    bool outlet_removal = false;          // Remove agents near outlet
    double outlet_removal_distance = 0.02; // Distance threshold for removal [m]
};

inline RuntimeCouplingConfig BuildRuntimeCouplingConfigFromParameters() {
    RuntimeCouplingConfig config{};
    const auto& store = GetParameterStore();

    auto apply_bool = [&](const std::string& key, bool* target) {
        if (store.Has<bool>(key)) {
            *target = store.Get<bool>(key);
        }
    };
    auto apply_string = [&](const std::string& key, std::string* target) {
        if (store.Has<std::string>(key)) {
            *target = store.Get<std::string>(key);
        }
    };

    // Coupling mode
    apply_bool("coupling.two_way", &config.two_way_coupling);
    apply_bool("coupling.enable_velocity", &config.enable_velocity);
    apply_bool("coupling.enable_pressure", &config.enable_pressure);
    apply_bool("coupling.enable_temperature", &config.enable_temperature);
    apply_bool("coupling.enable_oxygen", &config.enable_oxygen);
    apply_bool("coupling.enable_h2o2", &config.enable_h2o2);
    apply_bool("coupling.enable_no2", &config.enable_no2);
    apply_bool("coupling.enable_oxygen_sink", &config.enable_oxygen_sink);
    apply_bool("coupling.enable_h2o2_sink", &config.enable_h2o2_sink);
    apply_bool("coupling.enable_no2_sink", &config.enable_no2_sink);
    apply_bool("coupling.enable_drag_physics", &config.enable_drag_physics);
    apply_bool("coupling.enable_pressure_gradient", &config.enable_pressure_gradient);

    // Additional physics forces
    apply_bool("coupling.enable_gravity", &config.enable_gravity);
    apply_bool("coupling.enable_virtual_mass", &config.enable_virtual_mass);
    apply_bool("coupling.enable_lift_forces", &config.enable_lift_forces);
    apply_bool("coupling.enable_turbulence_dispersion", &config.enable_turbulence_dispersion);
    apply_bool("coupling.enable_turbulence", &config.enable_turbulence);

    // Simulation flags (can override coupling flags)
    apply_bool("simulation.force_2d", &config.force_2d);
    apply_bool("simulation.use_massless_tracers", &config.use_massless_tracers);
    apply_bool("simulation.enable_particle_motion", &config.enable_particle_motion);
    apply_bool("simulation.enable_wall_interactions", &config.enable_wall_interactions);
    apply_string("simulation.wall_confinement", &config.wall_confinement);
    apply_string("interaction.wall_velocity_response",
                 &config.wall_velocity_response);
    apply_bool("simulation.enable_particle_particle", &config.enable_particle_particle);

    // Environment parameters
    auto apply_double = [&](const std::string& key, double* target) {
        if (store.Has<double>(key)) {
            *target = store.Get<double>(key);
        }
    };
    apply_double("environment.gravity_x", &config.gravity_x);
    apply_double("environment.gravity_y", &config.gravity_y);
    apply_double("environment.gravity_z", &config.gravity_z);
    apply_double("agent.particle_density", &config.particle_density);
    apply_double("coupling.oxygen_uptake_rate", &config.oxygen_uptake_rate);
    apply_double("coupling.oxygen_uptake_max", &config.oxygen_uptake_max);

    // Inlet/Outlet agent management
    apply_bool("agent.inlet_generation", &config.inlet_generation);
    apply_double("agent.inlet_spawn_interval", &config.inlet_spawn_interval);
    auto apply_int = [&](const std::string& key, int* target) {
        if (store.Has<int>(key)) {
            *target = store.Get<int>(key);
        }
    };
    apply_int("agent.inlet_agents_per_spawn", &config.inlet_agents_per_spawn);
    apply_double("agent.inlet_zone_offset", &config.inlet_zone_offset);
    apply_bool("agent.outlet_removal", &config.outlet_removal);
    apply_double("agent.outlet_removal_distance", &config.outlet_removal_distance);

    return config;
}

inline const RuntimeCouplingConfig& GetRuntimeCouplingConfig() {
    static const RuntimeCouplingConfig config = BuildRuntimeCouplingConfigFromParameters();
    return config;
}

inline bool IsTwoWayCoupling() {
    return GetRuntimeCouplingConfig().two_way_coupling;
}

inline bool IsVelocityEnabled() {
    return GetRuntimeCouplingConfig().enable_velocity;
}

inline bool IsPressureEnabled() {
    return GetRuntimeCouplingConfig().enable_pressure;
}

inline bool IsDragPhysicsEnabled() {
    return GetRuntimeCouplingConfig().enable_drag_physics;
}

inline bool IsOxygenEnabled() {
    return GetRuntimeCouplingConfig().enable_oxygen;
}

inline bool IsOxygenSinkEnabled() {
    return GetRuntimeCouplingConfig().enable_oxygen_sink;
}

inline bool IsNO2SinkEnabled() {
    return GetRuntimeCouplingConfig().enable_no2_sink;
}

inline bool IsMasslessTracerMode() {
    return GetRuntimeCouplingConfig().use_massless_tracers;
}

inline bool IsForce2D() {
    return GetRuntimeCouplingConfig().force_2d;
}

inline bool IsWallInteractionsEnabled() {
    return GetRuntimeCouplingConfig().enable_wall_interactions;
}

inline bool IsWallConfinementEnabled() {
    const auto& mode = GetRuntimeCouplingConfig().wall_confinement;
    return mode == "enforce";
}

inline bool IsAabbWallConfinementEnabled() {
    return GetRuntimeCouplingConfig().wall_confinement == "aabb";
}

inline bool IsWallConfinementModeValid() {
    const auto& mode = GetRuntimeCouplingConfig().wall_confinement;
    return mode == "enforce" || mode == "aabb" || mode == "off";
}

inline bool IsWallVelocityResponseValid() {
    const auto& response =
        GetRuntimeCouplingConfig().wall_velocity_response;
    return response == "none" || response == "slip"
        || response == "tangential" || response == "rebound"
        || response == "stick";
}

enum class WallVelocityResponse {
    kNone,
    kSlip,
    kRebound,
    kStick
};

inline WallVelocityResponse GetWallVelocityResponse() {
    const auto& response =
        GetRuntimeCouplingConfig().wall_velocity_response;
    if (response == "none") return WallVelocityResponse::kNone;
    if (response == "slip" || response == "tangential") {
        return WallVelocityResponse::kSlip;
    }
    if (response == "stick") return WallVelocityResponse::kStick;
    return WallVelocityResponse::kRebound;
}

inline bool IsParticleParticleEnabled() {
    return GetRuntimeCouplingConfig().enable_particle_particle;
}

inline bool IsPressureGradientEnabled() {
    return GetRuntimeCouplingConfig().enable_pressure_gradient;
}

inline bool IsGravityEnabled() {
    return GetRuntimeCouplingConfig().enable_gravity;
}

inline bool IsVirtualMassEnabled() {
    return GetRuntimeCouplingConfig().enable_virtual_mass;
}

inline bool IsLiftForcesEnabled() {
    return GetRuntimeCouplingConfig().enable_lift_forces;
}

inline bool IsTurbulenceDispersionEnabled() {
    return GetRuntimeCouplingConfig().enable_turbulence_dispersion;
}

inline bool IsTurbulenceEnabled() {
    return GetRuntimeCouplingConfig().enable_turbulence;
}

inline std::array<double, 3> GravityVector() {
    const auto& config = GetRuntimeCouplingConfig();
    return {config.gravity_x, config.gravity_y, config.gravity_z};
}

inline double ParticleDensity() {
    return GetRuntimeCouplingConfig().particle_density;
}

inline double OxygenUptakeRate() {
    return GetRuntimeCouplingConfig().oxygen_uptake_rate;
}

inline double OxygenUptakeMax() {
    return GetRuntimeCouplingConfig().oxygen_uptake_max;
}

// =========================================================================
// Inlet/Outlet Agent Management Helpers
// =========================================================================

inline bool IsInletGenerationEnabled() {
    return GetRuntimeCouplingConfig().inlet_generation;
}

inline double InletSpawnInterval() {
    return GetRuntimeCouplingConfig().inlet_spawn_interval;
}

inline int InletAgentsPerSpawn() {
    return GetRuntimeCouplingConfig().inlet_agents_per_spawn;
}

inline double InletZoneOffset() {
    return GetRuntimeCouplingConfig().inlet_zone_offset;
}

inline bool IsOutletRemovalEnabled() {
    return GetRuntimeCouplingConfig().outlet_removal;
}

inline double OutletRemovalDistance() {
    return GetRuntimeCouplingConfig().outlet_removal_distance;
}

inline CavityConfig::DomainConfig BuildDomainConfigFromParameters() {
    CavityConfig::DomainConfig domain{};
    // Note: domain bounds remain NaN until set from mesh bounds
    // This ensures we don't accidentally use hardcoded values
    return domain;
}

inline CavityConfig::AgentConfig BuildAgentConfigFromParameters() {
    CavityConfig::AgentConfig agent{};
    const auto& store = GetParameterStore();
    auto apply_double = [&](const std::string& key, double* target) {
        if (store.Has<double>(key)) {
            *target = store.Get<double>(key);
        }
    };
    auto apply_bool = [&](const std::string& key, bool* target) {
        if (store.Has<bool>(key)) {
            *target = store.Get<bool>(key);
        }
    };
    apply_double("agent.default_diameter", &agent.initial_diameter);
    apply_double("agent.default_temperature", &agent.default_temperature);
    apply_double("agent.default_pressure", &agent.default_pressure);
    apply_bool("agent.place_at_cell_centers", &agent.place_at_cell_centers);
    apply_bool("agent.initial_stationary", &agent.initial_stationary);
    return agent;
}

inline const CavityConfig& PreciceDefaultConfig() {
    static const CavityConfig config{};
    return config;
}

inline const CavityConfig::DomainConfig& PreciceDomain() {
    static const CavityConfig::DomainConfig domain = BuildDomainConfigFromParameters();
    return domain;
}

inline const CavityConfig::AgentConfig& PreciceAgentConfig() {
    static const CavityConfig::AgentConfig agent = BuildAgentConfigFromParameters();
    return agent;
}

inline CavityConfig::VisualizationConfig BuildVisualizationConfigFromParameters() {
    CavityConfig::VisualizationConfig viz{};
    const auto& store = GetParameterStore();
    auto apply_double = [&](const std::string& key, double* target) {
        if (store.Has<double>(key)) {
            *target = store.Get<double>(key);
        }
    };
    apply_double("visualization.temperature_scale_min", &viz.temperature_scale_min);
    apply_double("visualization.temperature_scale_max", &viz.temperature_scale_max);
    return viz;
}

inline const CavityConfig::VisualizationConfig& PreciceVisualization() {
    static const CavityConfig::VisualizationConfig viz = BuildVisualizationConfigFromParameters();
    return viz;
}

inline CavityConfig::PreciceConfig BuildPreciceConfigFromParameters() {
    CavityConfig::PreciceConfig precice{};
    const auto& store = GetParameterStore();
    auto apply_string = [&](const std::string& key, std::string* target) {
        if (store.Has<std::string>(key)) {
            *target = store.Get<std::string>(key);
        }
    };
    apply_string("precice.config_file", &precice.config_file);
    apply_string("precice.participant_name", &precice.participant_name);
    apply_string("mesh.primary_volume_mesh", &precice.volume_mesh);
    apply_string("mesh.particle_mesh", &precice.particle_mesh);
    apply_string("precice.velocity_data_name", &precice.velocity_field);
    apply_string("precice.pressure_data_name", &precice.pressure_field);
    apply_string("precice.temperature_data_name", &precice.temperature_field);
    apply_string("precice.position_data_name", &precice.position_field);
    return precice;
}

inline const CavityConfig::PreciceConfig& PreciceExchange() {
    static const CavityConfig::PreciceConfig precice = BuildPreciceConfigFromParameters();
    return precice;
}

inline CavityConfig::InteractionConfig BuildInteractionConfigFromParameters() {
    CavityConfig::InteractionConfig interaction{};
    const auto& store = GetParameterStore();
    auto apply_double = [&](const std::string& key, double* target) {
        if (store.Has<double>(key)) {
            *target = store.Get<double>(key);
        }
    };
    apply_double("interaction.youngs_modulus", &interaction.young_modulus);
    apply_double("interaction.poisson_ratio", &interaction.poisson_ratio);
    apply_double("interaction.density", &interaction.density);
    apply_double("interaction.friction_coefficient", &interaction.friction_coefficient);
    apply_double("interaction.restitution_coefficient", &interaction.restitution_coefficient);
    apply_double("interaction.damping_factor", &interaction.damping_factor);
    apply_double("interaction.cohesion_energy", &interaction.cohesion_energy);
    apply_double("interaction.wall_youngs_modulus", &interaction.wall_young_modulus);
    apply_double("interaction.wall_poisson_ratio", &interaction.wall_poisson_ratio);
    apply_double("interaction.interaction_range", &interaction.interaction_range);
    apply_double("interaction.force_ramp_time", &interaction.force_ramp_time);
    return interaction;
}

inline const CavityConfig::InteractionConfig& PreciceInteractionConfig() {
    static const CavityConfig::InteractionConfig interaction = BuildInteractionConfigFromParameters();
    return interaction;
}

inline CavityConfig::FluidConfig BuildFluidConfigFromParameters() {
    CavityConfig::FluidConfig fluid{};
    fluid.density = FluidDensity();
    fluid.dynamic_viscosity = FluidViscosity();
    return fluid;
}

inline const CavityConfig::FluidConfig& PreciceFluidConfig() {
    static const CavityConfig::FluidConfig fluid = BuildFluidConfigFromParameters();
    return fluid;
}

} // namespace bdm

#endif // COUPLING_CONFIG_H_
