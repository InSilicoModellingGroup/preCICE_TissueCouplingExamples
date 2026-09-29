#ifndef BDM_ADAPTER_COUPLING_CONFIGURATION_H_
#define BDM_ADAPTER_COUPLING_CONFIGURATION_H_

#include <array>
#include <cmath>
#include <limits>
#include <string>

namespace bdm_adapter {

/**
 * @brief Compile-time coupling configuration using template parameters
 *
 * This template allows enabling/disabling coupling features at compile time,
 * resulting in zero-overhead abstractions for unused features.
 *
 * Template parameters:
 * @tparam EnableVelocity      Enable velocity field coupling
 * @tparam EnablePressure      Enable pressure field coupling
 * @tparam EnableTemperature   Enable temperature field coupling
 * @tparam EnableTwoWay        Enable two-way coupling (particle -> fluid feedback)
 * @tparam EnableInteractions  Enable particle-particle interactions
 * @tparam EnableMotion        Enable particle motion integration
 * @tparam Force2D             Force 2D simulation (zero z-velocity)
 * @tparam EnableDragPhysics   Enable drag force physics
 * @tparam UseMasslessTracers  Use massless tracer mode (v_particle = v_fluid)
 * @tparam EnablePressureGradientForce Enable pressure gradient force
 * @tparam EnableTurbulence   Enable turbulence field coupling (k, epsilon/omega, nut)
 */
template<
    bool EnableVelocity = true,
    bool EnablePressure = true,
    bool EnableTemperature = false,
    bool EnableTwoWay = false,
    bool EnableInteractions = false,
    bool EnableMotion = true,
    bool Force2D = true,
    bool EnableDragPhysics = true,
    bool UseMasslessTracers = false,
    bool EnablePressureGradientForce = true,
    bool EnableTurbulence = false
>
struct CouplingConfig {
    // Compile-time feature flags
    static constexpr bool enable_velocity = EnableVelocity;
    static constexpr bool enable_pressure = EnablePressure;
    static constexpr bool enable_temperature = EnableTemperature;
    static constexpr bool enable_two_way = EnableTwoWay;
    static constexpr bool enable_interactions = EnableInteractions;
    static constexpr bool enable_motion = EnableMotion;
    static constexpr bool force_2d = Force2D;
    static constexpr bool enable_drag_physics = EnableDragPhysics;
    static constexpr bool use_massless_tracers = UseMasslessTracers;
    static constexpr bool enable_pressure_gradient_force = EnablePressureGradientForce;
    static constexpr bool enable_turbulence = EnableTurbulence;

    /**
     * @brief Fluid properties for drag and pressure gradient calculations
     */
    struct FluidConfig {
        double density = std::numeric_limits<double>::quiet_NaN();           // Set by preCICE
        double dynamic_viscosity = std::numeric_limits<double>::quiet_NaN(); // Set by preCICE

        double KinematicViscosity() const {
            return dynamic_viscosity / density;
        }
    };

    /**
     * @brief Domain configuration
     */
    struct DomainConfig {
        std::array<double, 3> min = {0.0, 0.0, 0.0};
        std::array<double, 3> max = {0.1, 0.1, 0.01};
        double z_plane = 0.005;  // Z-plane for 2D simulations
    };

    /**
     * @brief Visualization configuration
     */
    struct VisualizationConfig {
        bool enabled = true;
        double export_interval = 0.1;
        std::string output_directory = "output";
    };

    // Runtime configurations
    FluidConfig fluid;
    DomainConfig domain;
    VisualizationConfig visualization;

    /**
     * @brief Validate configuration consistency
     */
    bool Validate() const {
        if constexpr (use_massless_tracers && enable_drag_physics) {
            // Warning: massless tracers ignore drag physics
            return true; // Still valid, just redundant
        }
        if (std::isnan(fluid.density) || std::isnan(fluid.dynamic_viscosity) ||
            fluid.density <= 0.0 || fluid.dynamic_viscosity <= 0.0) {
            return false;
        }
        return true;
    }

    /**
     * @brief Get configuration description string
     */
    std::string Description() const {
        std::string desc = "CouplingConfig: ";
        if constexpr (enable_velocity) desc += "V ";
        if constexpr (enable_pressure) desc += "P ";
        if constexpr (enable_temperature) desc += "T ";
        if constexpr (enable_two_way) desc += "2Way ";
        if constexpr (enable_interactions) desc += "Inter ";
        if constexpr (enable_motion) desc += "Motion ";
        if constexpr (force_2d) desc += "2D ";
        if constexpr (enable_drag_physics) desc += "Drag ";
        if constexpr (use_massless_tracers) desc += "Tracer ";
        if constexpr (enable_pressure_gradient_force) desc += "GradP ";
        if constexpr (enable_turbulence) desc += "Turb ";
        return desc;
    }
};

// ============================================================================
// Predefined configurations for common use cases
// ============================================================================

/** @brief One-way coupling: fluid affects particles, particles don't affect fluid */
using OneWayCouplingConfig = CouplingConfig<
    true,   // velocity
    true,   // pressure
    false,  // temperature
    false,  // two-way
    false,  // interactions
    true,   // motion
    true,   // 2D
    true,   // drag physics
    false,  // massless tracers
    true,   // pressure gradient force
    false   // turbulence
>;

/** @brief Massless tracer mode: particles follow fluid exactly */
using TracerConfig = CouplingConfig<
    true,   // velocity
    true,   // pressure (for diagnostics)
    false,  // temperature
    false,  // two-way
    false,  // interactions
    true,   // motion
    true,   // 2D
    false,  // drag physics (disabled)
    true,   // massless tracers (enabled)
    false,  // pressure gradient force (disabled)
    false   // turbulence
>;

/** @brief Full physics: drag, pressure gradient, interactions */
using FullPhysicsConfig = CouplingConfig<
    true,   // velocity
    true,   // pressure
    false,  // temperature
    false,  // two-way (still one-way for now)
    true,   // interactions
    true,   // motion
    true,   // 2D
    true,   // drag physics
    false,  // massless tracers
    true,   // pressure gradient force
    false   // turbulence
>;

/** @brief Two-way coupling (future) */
using TwoWayCouplingConfig = CouplingConfig<
    true,   // velocity
    true,   // pressure
    false,  // temperature
    true,   // two-way
    true,   // interactions
    true,   // motion
    true,   // 2D
    true,   // drag physics
    false,  // massless tracers
    true,   // pressure gradient force
    false   // turbulence
>;

/** @brief 3D simulation config */
using Config3D = CouplingConfig<
    true,   // velocity
    true,   // pressure
    false,  // temperature
    false,  // two-way
    false,  // interactions
    true,   // motion
    false,  // 3D (not forced 2D)
    true,   // drag physics
    false,  // massless tracers
    true,   // pressure gradient force
    false   // turbulence
>;

/** @brief Turbulence-aware coupling: includes k, epsilon/omega, nut fields */
using TurbulenceCouplingConfig = CouplingConfig<
    true,   // velocity
    true,   // pressure
    false,  // temperature
    false,  // two-way
    false,  // interactions
    true,   // motion
    true,   // 2D
    true,   // drag physics
    false,  // massless tracers
    true,   // pressure gradient force
    true    // turbulence (enabled)
>;

// Default configuration alias
using DefaultCouplingConfig = OneWayCouplingConfig;

} // namespace bdm_adapter

#endif // BDM_ADAPTER_COUPLING_CONFIGURATION_H_
