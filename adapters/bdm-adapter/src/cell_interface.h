#ifndef BDM_ADAPTER_CELL_INTERFACE_H_
#define BDM_ADAPTER_CELL_INTERFACE_H_

#include <array>
#include <type_traits>

namespace bdm_adapter {

// Use std::array for interface - compatible with bdm::Real3
using Real3 = std::array<double, 3>;

/**
 * @brief Documentation of required interface for coupled cells
 *
 * Any BioDynaMo cell type that wants to participate in preCICE coupling
 * must implement these methods. The adapter uses duck typing - if your
 * cell has these methods, it will work.
 *
 * Required methods:
 * - GetPosition() const -> Real3
 * - SetFluidVelocity(const Real3&)
 * - GetFluidVelocity() const -> Real3
 * - SetPressure(double)
 * - GetPressure() const -> double
 * - SetPressureGradient(const std::array<double, 3>&)
 * - GetPressureGradient() const -> std::array<double, 3>
 * - IsStationary() const -> bool
 * - GetUid() const -> AgentUid (inherited from bdm::Cell)
 *
 * Optional methods for extended features:
 * - SetTemperature(double) / GetTemperature() const
 * - SetOxygen(double) / GetOxygen() const
 * - SetNO2(double) / GetNO2() const
 * - SetOxygenUptake(double) / GetOxygenUptake() const
 * - SetVelocity(const Real3&) / GetVelocity() const
 * - AddForce(const Real3&) / ClearForces() / GetAccumulatedForce() const
 *
 * Example implementation:
 * @code
 * class MyCell : public bdm::Cell {
 *   void SetFluidVelocity(const Real3& v) { fluid_velocity_ = v; }
 *   const Real3& GetFluidVelocity() const { return fluid_velocity_; }
 *   void SetPressure(double p) { pressure_ = p; }
 *   double GetPressure() const { return pressure_; }
 *   void SetPressureGradient(const std::array<double, 3>& g) { grad_p_ = g; }
 *   std::array<double, 3> GetPressureGradient() const { return grad_p_; }
 *   bool IsStationary() const { return is_stationary_; }
 * private:
 *   Real3 fluid_velocity_{0, 0, 0};
 *   double pressure_ = 0.0;
 *   std::array<double, 3> grad_p_{0, 0, 0};
 *   bool is_stationary_ = false;
 * };
 * @endcode
 */

// Type traits for checking cell capabilities (C++17 compatible)
namespace detail {

// Check if type has SetFluidVelocity method
template<typename T, typename = void>
struct has_set_fluid_velocity : std::false_type {};

template<typename T>
struct has_set_fluid_velocity<T,
    std::void_t<decltype(std::declval<T&>().SetFluidVelocity(std::declval<const Real3&>()))>
> : std::true_type {};

// Check if type has SetPressure method
template<typename T, typename = void>
struct has_set_pressure : std::false_type {};

template<typename T>
struct has_set_pressure<T,
    std::void_t<decltype(std::declval<T&>().SetPressure(std::declval<double>()))>
> : std::true_type {};

// Check if type has SetTemperature method
template<typename T, typename = void>
struct has_set_temperature : std::false_type {};

template<typename T>
struct has_set_temperature<T,
    std::void_t<decltype(std::declval<T&>().SetTemperature(std::declval<double>()))>
> : std::true_type {};

// Check if type has SetOxygen method
template<typename T, typename = void>
struct has_set_oxygen : std::false_type {};

template<typename T>
struct has_set_oxygen<T,
    std::void_t<decltype(std::declval<T&>().SetOxygen(std::declval<double>()))>
> : std::true_type {};

// Check if type has GetOxygen method
template<typename T, typename = void>
struct has_get_oxygen : std::false_type {};

template<typename T>
struct has_get_oxygen<T,
    std::void_t<decltype(std::declval<const T&>().GetOxygen())>
> : std::true_type {};

// Check if type has GetNO2 method
template<typename T, typename = void>
struct has_get_no2 : std::false_type {};

template<typename T>
struct has_get_no2<T,
    std::void_t<decltype(std::declval<const T&>().GetNO2())>
> : std::true_type {};

template<typename T, typename = void>
struct has_set_no2 : std::false_type {};

template<typename T>
struct has_set_no2<T,
    std::void_t<decltype(std::declval<T&>().SetNO2(std::declval<double>()))>
> : std::true_type {};

template<typename T, typename = void>
struct has_get_h2o2 : std::false_type {};

template<typename T>
struct has_get_h2o2<T,
    std::void_t<decltype(std::declval<const T&>().GetH2O2())>
> : std::true_type {};

template<typename T, typename = void>
struct has_set_h2o2 : std::false_type {};

template<typename T>
struct has_set_h2o2<T,
    std::void_t<decltype(std::declval<T&>().SetH2O2(std::declval<double>()))>
> : std::true_type {};

template<typename T, typename = void>
struct has_set_h2o2_uptake : std::false_type {};

template<typename T>
struct has_set_h2o2_uptake<T,
    std::void_t<decltype(std::declval<T&>().SetH2O2Uptake(std::declval<double>()))>
> : std::true_type {};

// Check if type has SetOxygenUptake method
template<typename T, typename = void>
struct has_set_oxygen_uptake : std::false_type {};

template<typename T>
struct has_set_oxygen_uptake<T,
    std::void_t<decltype(std::declval<T&>().SetOxygenUptake(std::declval<double>()))>
> : std::true_type {};

// Check if type has GetOxygenUptake method
template<typename T, typename = void>
struct has_get_oxygen_uptake : std::false_type {};

template<typename T>
struct has_get_oxygen_uptake<T,
    std::void_t<decltype(std::declval<const T&>().GetOxygenUptake())>
> : std::true_type {};

// Check if type has SetPressureGradient method
template<typename T, typename = void>
struct has_set_pressure_gradient : std::false_type {};

template<typename T>
struct has_set_pressure_gradient<T,
    std::void_t<decltype(std::declval<T&>().SetPressureGradient(std::declval<const std::array<double,3>&>()))>
> : std::true_type {};

// Check if type has SetFluidAcceleration method
template<typename T, typename = void>
struct has_set_fluid_acceleration : std::false_type {};

template<typename T>
struct has_set_fluid_acceleration<T,
    std::void_t<decltype(std::declval<T&>().SetFluidAcceleration(std::declval<const Real3&>()))>
> : std::true_type {};

// Check if type has SetFluidVorticity method
template<typename T, typename = void>
struct has_set_fluid_vorticity : std::false_type {};

template<typename T>
struct has_set_fluid_vorticity<T,
    std::void_t<decltype(std::declval<T&>().SetFluidVorticity(std::declval<const Real3&>()))>
> : std::true_type {};

// Check if type has IsStationary method
template<typename T, typename = void>
struct has_is_stationary : std::false_type {};

template<typename T>
struct has_is_stationary<T,
    std::void_t<decltype(std::declval<const T&>().IsStationary())>
> : std::true_type {};

// Check if type has SetTKE method (turbulent kinetic energy)
template<typename T, typename = void>
struct has_set_tke : std::false_type {};

template<typename T>
struct has_set_tke<T,
    std::void_t<decltype(std::declval<T&>().SetTKE(std::declval<double>()))>
> : std::true_type {};

// Check if type has SetEpsilon method (turbulent dissipation rate)
template<typename T, typename = void>
struct has_set_epsilon : std::false_type {};

template<typename T>
struct has_set_epsilon<T,
    std::void_t<decltype(std::declval<T&>().SetEpsilon(std::declval<double>()))>
> : std::true_type {};

// Check if type has SetOmega method (specific dissipation rate)
template<typename T, typename = void>
struct has_set_omega : std::false_type {};

template<typename T>
struct has_set_omega<T,
    std::void_t<decltype(std::declval<T&>().SetOmega(std::declval<double>()))>
> : std::true_type {};

// Check if type has SetNut method (turbulent viscosity)
template<typename T, typename = void>
struct has_set_nut : std::false_type {};

template<typename T>
struct has_set_nut<T,
    std::void_t<decltype(std::declval<T&>().SetNut(std::declval<double>()))>
> : std::true_type {};

} // namespace detail

// Convenience aliases
template<typename T>
inline constexpr bool has_fluid_velocity_v = detail::has_set_fluid_velocity<T>::value;

template<typename T>
inline constexpr bool has_pressure_v = detail::has_set_pressure<T>::value;

template<typename T>
inline constexpr bool has_temperature_v = detail::has_set_temperature<T>::value;

template<typename T>
inline constexpr bool has_oxygen_setter_v = detail::has_set_oxygen<T>::value;

template<typename T>
inline constexpr bool has_oxygen_getter_v = detail::has_get_oxygen<T>::value;

template<typename T>
inline constexpr bool has_no2_getter_v = detail::has_get_no2<T>::value;

template<typename T>
inline constexpr bool has_no2_setter_v = detail::has_set_no2<T>::value;

template<typename T>
inline constexpr bool has_h2o2_getter_v = detail::has_get_h2o2<T>::value;

template<typename T>
inline constexpr bool has_h2o2_setter_v = detail::has_set_h2o2<T>::value;

template<typename T>
inline constexpr bool has_h2o2_uptake_setter_v = detail::has_set_h2o2_uptake<T>::value;

template<typename T>
inline constexpr bool has_oxygen_uptake_setter_v = detail::has_set_oxygen_uptake<T>::value;

template<typename T>
inline constexpr bool has_oxygen_uptake_getter_v = detail::has_get_oxygen_uptake<T>::value;

template<typename T>
inline constexpr bool has_pressure_gradient_v = detail::has_set_pressure_gradient<T>::value;

template<typename T>
inline constexpr bool has_fluid_acceleration_v = detail::has_set_fluid_acceleration<T>::value;

template<typename T>
inline constexpr bool has_fluid_vorticity_v = detail::has_set_fluid_vorticity<T>::value;

template<typename T>
inline constexpr bool has_is_stationary_v = detail::has_is_stationary<T>::value;

template<typename T>
inline constexpr bool has_tke_v = detail::has_set_tke<T>::value;

template<typename T>
inline constexpr bool has_epsilon_v = detail::has_set_epsilon<T>::value;

template<typename T>
inline constexpr bool has_omega_v = detail::has_set_omega<T>::value;

template<typename T>
inline constexpr bool has_nut_v = detail::has_set_nut<T>::value;

} // namespace bdm_adapter

#endif // BDM_ADAPTER_CELL_INTERFACE_H_
