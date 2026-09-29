#ifndef BDM_FLUID_PROPERTIES_H_
#define BDM_FLUID_PROPERTIES_H_

#include <limits>
#include <cmath>

#include "core/util/log.h"

namespace bdm {

namespace detail {
    struct FluidProperties {
        double density = std::numeric_limits<double>::quiet_NaN();
        double dynamic_viscosity = std::numeric_limits<double>::quiet_NaN();
        bool is_set = false;
    };

    inline FluidProperties& MutableFluidProperties() {
        static FluidProperties props;
        return props;
    }
}

inline void EnsureFluidPropertiesReady() {
    auto& props = detail::MutableFluidProperties();
    if (!props.is_set || std::isnan(props.density) || std::isnan(props.dynamic_viscosity)) {
        Log::Fatal("FluidProperties",
                   "Fluid properties not set. They must be provided by OpenFOAM via preCICE "
                   "(rho/mu). Ensure FetchFluidProperties succeeds before accessing them.");
    }
}

// Runtime fluid properties (populated from preCICE data or params.csv fallback)
inline double FluidDensity() { EnsureFluidPropertiesReady(); return detail::MutableFluidProperties().density; }
inline double FluidViscosity() { EnsureFluidPropertiesReady(); return detail::MutableFluidProperties().dynamic_viscosity; }
inline double KinematicViscosity() { return FluidViscosity() / FluidDensity(); }

inline void SetFluidProperties(double density, double dynamic_viscosity) {
    auto& props = detail::MutableFluidProperties();
    props.density = density;
    props.dynamic_viscosity = dynamic_viscosity;
    props.is_set = true;
    Log::Info("FluidProperties", "Set fluid properties: rho=",
              props.density, " kg/m^3, mu=", props.dynamic_viscosity, " Pa.s");
}

}  // namespace bdm

#endif  // BDM_FLUID_PROPERTIES_H_
