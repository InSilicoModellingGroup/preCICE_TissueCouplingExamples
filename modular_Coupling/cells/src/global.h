// =============================================================================
//
//   Coupled Cells - BioDynaMo-OpenFOAM-preCICE coupling
//   Following ABM4bio methodology
//
//   This file provides:
//   1. Common includes
//   2. Math utilities (ABM4bio style)
//   3. Progress bar utilities
//   4. Re-exports coupling configuration functions
//
// =============================================================================

#ifndef CELLS_GLOBAL_H_
#define CELLS_GLOBAL_H_

// =============================================================================
// Standard Library Includes
// =============================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// =============================================================================
// BioDynaMo Includes
// =============================================================================
#include "biodynamo.h"
#include "core/util/log.h"
#include "core/resource_manager.h"
#include "fluid_properties.h"

// =============================================================================
// Project Configuration - SINGLE SOURCE OF TRUTH
// coupling_configuration.h defines: IsTwoWayCoupling(), IsDragPhysicsEnabled(), etc.
// =============================================================================
#include "config/coupling_configuration.h"

// =============================================================================
// Namespace
// =============================================================================
namespace bdm {

// =============================================================================
// Type Aliases (ABM4bio style)
// =============================================================================
using Double3x3 = MathArray<MathArray<double, 3>, 3>;

// =============================================================================
// Error Handling Macros (ABM4bio style)
// =============================================================================
#ifndef ABORT_
#define ABORT_(_msg_) \
    bdm::Log::Fatal(" *** @file " + std::string(__FILE__) + " @line " + std::to_string(__LINE__), \
                    std::string(_msg_));
#endif

#ifndef ASSERT_
#define ASSERT_(_cond_, _msg_) \
    if (!(_cond_)) \
        bdm::Log::Fatal(" *** @file " + std::string(__FILE__) + " @line " + std::to_string(__LINE__), \
                        std::string(_msg_));
#endif

// =============================================================================
// Math Utilities (ABM4bio style)
// =============================================================================
inline double pow2(const double& v) { return v * v; }
inline double pow3(const double& v) { return v * pow2(v); }
inline double pow4(const double& v) { return pow2(pow2(v)); }

inline double L2norm(const Double3& v) {
    return std::sqrt(pow2(v[0]) + pow2(v[1]) + pow2(v[2]));
}

// Alias for compatibility
inline double Magnitude(const Real3& v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

inline Double3 normalize(const Double3& v, double tol = 1.0e-12) {
    double mag = L2norm(v);
    if (mag > tol) {
        return {v[0] / mag, v[1] / mag, v[2] / mag};
    }
    return {0.0, 0.0, 0.0};
}

// Alias for compatibility
inline Real3 Normalize(const Real3& v) {
    double mag = Magnitude(v);
    if (mag < 1e-12) return {0.0, 0.0, 0.0};
    return {v[0]/mag, v[1]/mag, v[2]/mag};
}

inline double Dot(const Real3& a, const Real3& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline Real3 Cross(const Real3& a, const Real3& b) {
    return {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

inline double Distance(const Real3& a, const Real3& b) {
    Real3 diff = {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
    return Magnitude(diff);
}

inline double clamp(double value, double min_val, double max_val) {
    return std::max(min_val, std::min(max_val, value));
}

// =============================================================================
// Parameter Access Helper (for direct params.csv queries)
// =============================================================================
template<typename T>
inline T GetParam(const std::string& key, T default_value) {
    const auto& store = GetParameterStore();
    return store.Has<T>(key) ? store.Get<T>(key) : default_value;
}

// Overload: require explicit configuration in params.csv (no default)
template<typename T>
inline T GetParam(const std::string& key) {
    const auto& store = GetParameterStore();
    if (!store.Has<T>(key)) {
        throw std::runtime_error("Required parameter not found in params.csv: " + key);
    }
    return store.Get<T>(key);
}
// =============================================================================
// Fluid Properties (runtime - see fluid_properties.h)
// =============================================================================

// =============================================================================
// Domain Bounds - Runtime storage for mesh-derived bounds
// These are set at runtime from adapter.GetDomainBoundsFromMesh()
// =============================================================================
namespace detail {
    // Runtime storage for domain bounds (set from preCICE mesh)
    // NO DEFAULT VALUES - must be explicitly set from preCICE mesh data
    // Initialized to NaN to catch uninitialized usage
    struct MeshDerivedBounds {
        double x_min = std::numeric_limits<double>::quiet_NaN();
        double x_max = std::numeric_limits<double>::quiet_NaN();
        double y_min = std::numeric_limits<double>::quiet_NaN();
        double y_max = std::numeric_limits<double>::quiet_NaN();
        double z_min = std::numeric_limits<double>::quiet_NaN();
        double z_max = std::numeric_limits<double>::quiet_NaN();
        bool is_set = false;
    };
    inline MeshDerivedBounds& GetMeshBounds() {
        static MeshDerivedBounds bounds;
        return bounds;
    }
}

// Set domain bounds from mesh (called by cells.h after adapter init)
inline void SetDomainBoundsFromMesh(double x_min, double x_max,
                                    double y_min, double y_max,
                                    double z_min, double z_max) {
    auto& bounds = detail::GetMeshBounds();
    bounds.x_min = x_min; bounds.x_max = x_max;
    bounds.y_min = y_min; bounds.y_max = y_max;
    bounds.z_min = z_min; bounds.z_max = z_max;
    bounds.is_set = true;
}

inline bool AreDomainBoundsSet() { return detail::GetMeshBounds().is_set; }

// Domain bound accessors (return mesh-derived values)
inline double DomainXMin() { return detail::GetMeshBounds().x_min; }
inline double DomainXMax() { return detail::GetMeshBounds().x_max; }
inline double DomainYMin() { return detail::GetMeshBounds().y_min; }
inline double DomainYMax() { return detail::GetMeshBounds().y_max; }
inline double DomainZMin() { return detail::GetMeshBounds().z_min; }
inline double DomainZMax() { return detail::GetMeshBounds().z_max; }
inline double DomainZPlane() {
    return (detail::GetMeshBounds().z_min + detail::GetMeshBounds().z_max) / 2.0;
}

// =============================================================================
// Agent Properties (from params.csv)
// =============================================================================
inline double DefaultDiameter() {
    return GetParam<double>("agent.default_diameter", 0.004);
}

inline double DefaultMass() {
    // If user explicitly sets agent.default_mass, honor it; otherwise derive from density and diameter
    const auto& store = GetParameterStore();
    if (store.Has<double>("agent.default_mass")) {
        return store.Get<double>("agent.default_mass");
    }
    const double d = DefaultDiameter();
    const double rho_p = GetParam<double>("agent.particle_density"); // Agent particle density (required in params.csv)
    const double pi = 3.14159265358979323846;
    const double vol = (pi / 6.0) * d * d * d;
    return rho_p * vol;
}

// =============================================================================
// Interaction Properties (from params.csv)
// =============================================================================
inline double YoungsModulus() {
    return GetParam<double>("interaction.youngs_modulus");
}

inline double PoissonRatio() {
    return GetParam<double>("interaction.poisson_ratio");
}

inline double RestitutionCoefficient() {
    return GetParam<double>("interaction.restitution_coefficient");
}

inline double FrictionCoefficient() {
    return GetParam<double>("interaction.friction_coefficient");
}

inline double WallRestitution() {
    return GetParam<double>("interaction.wall_restitution");
}

// =============================================================================
// Visualization Settings (from params.csv)
// =============================================================================
inline bool ExportVisualization() {
    return GetParam<bool>("visualization.export_visualization");
}

inline int VisualizationExportInterval() {
    return GetParam<int>("visualization.export_interval");
}

inline std::string VisualizationOutputDirectory() {
    return GetParam<std::string>("visualization.output_directory");
}

inline std::string VisualizationSimulationName() {
    return GetParam<std::string>("visualization.simulation_name");
}

// =============================================================================
// Progress Bar (ABM4bio style)
// =============================================================================
inline void StatusBar(std::ostream& os, int x, int n, int bar_width = 50) {
    const double progress = static_cast<double>(x) / static_cast<double>(n);
    const int pos = static_cast<int>(bar_width * progress);
    os << "[";
    for (int i = 0; i < bar_width; i++) {
        if (i < pos) os << "=";
        else if (i == pos) os << ">";
        else os << " ";
    }
    os << "] " << std::fixed << std::setprecision(1) << (progress * 100.0) << "%\r";
    os << std::flush;
    if (x == n) {
        os << std::string(bar_width + 20, ' ') << '\r' << std::flush;
    }
}

inline void TimeStatusBar(std::ostream& os, int x, int n, double t, int bar_width = 50) {
    const double progress = static_cast<double>(x) / static_cast<double>(n);
    const int pos = static_cast<int>(bar_width * progress);
    os << "[";
    for (int i = 0; i < bar_width; i++) {
        if (i < pos) os << "=";
        else if (i == pos) os << ">";
        else os << " ";
    }
    os << "] " << std::fixed << std::setprecision(1) << (progress * 100.0) << "% ";
    os << "(t=" << std::scientific << std::setprecision(3) << t << ")\r";
    os << std::flush;
}

// =============================================================================
// Logging Helpers
// =============================================================================
inline void PrintCouplingConfig() {
    std::cout << "\n=== Coupling Configuration (runtime) ===" << std::endl;
    std::cout << "  coupling.two_way:              " << (IsTwoWayCoupling() ? "true" : "false") << std::endl;
    std::cout << "  coupling.enable_drag_physics:  " << (IsDragPhysicsEnabled() ? "true" : "false") << std::endl;
    std::cout << "  coupling.enable_oxygen:        " << (IsOxygenEnabled() ? "true" : "false") << std::endl;
    std::cout << "  coupling.enable_oxygen_sink:   " << (IsOxygenSinkEnabled() ? "true" : "false") << std::endl;
    std::cout << "  coupling.enable_no2_sink:      " << (IsNO2SinkEnabled() ? "true" : "false") << std::endl;
    std::cout << "  simulation.use_massless_tracers: " << (IsMasslessTracerMode() ? "true" : "false") << std::endl;
    std::cout << "  simulation.force_2d:           " << (IsForce2D() ? "true" : "false") << std::endl;
    std::cout << "  simulation.wall_confinement:  "
              << GetRuntimeCouplingConfig().wall_confinement << std::endl;
    std::cout << "  interaction.wall_velocity_response: "
              << GetRuntimeCouplingConfig().wall_velocity_response << std::endl;
    std::cout << "\n=== Fluid Properties (runtime) ===" << std::endl;
    std::cout << "  density (rho):          " << FluidDensity() << " kg/m^3" << std::endl;
    std::cout << "  dynamic_viscosity (mu): " << FluidViscosity() << " Pa.s" << std::endl;
    std::cout << "\n=== Domain Bounds (from preCICE FVMmesh) ===" << std::endl;
    if (AreDomainBoundsSet()) {
        std::cout << "  X: [" << DomainXMin() << ", " << DomainXMax() << "]" << std::endl;
        std::cout << "  Y: [" << DomainYMin() << ", " << DomainYMax() << "]" << std::endl;
        std::cout << "  Z: [" << DomainZMin() << ", " << DomainZMax() << "]" << std::endl;
        std::cout << "  Z plane: " << DomainZPlane() << std::endl;
    } else {
        std::cout << "  [Not yet set - waiting for preCICE mesh initialization]" << std::endl;
    }
    std::cout << "================================================\n" << std::endl;
}

} // namespace bdm

#endif // CELLS_GLOBAL_H_
