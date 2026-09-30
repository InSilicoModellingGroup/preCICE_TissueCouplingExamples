#ifndef RONS_INTRACELLULAR_MODEL_H_
#define RONS_INTRACELLULAR_MODEL_H_

#include <algorithm>

namespace bdm {

/// Per-cell intracellular state used by the ABM4bio RONS-response ODEs.
struct IntracellularRonsState {
  double h2o2_load = 0.0;  // R_H2O2 [mol/m3]
  double no2_load = 0.0;   // R_NO2  [mol/m3]
  double damage = 0.0;     // D, dimensionless damage index
};

/// Shared coefficients for H2O2 and NO2-.
struct IntracellularRonsCoefficients {
  double uptake_rate = 0.0;              // alpha_ext [1/s]
  double clearance_rate = 0.0;           // beta_R [1/s]
  double damage_accumulation_rate = 0.0; // alpha_D [m3/(mol s)]
  double damage_repair_rate = 0.0;        // beta_D [1/s]
};

struct IntracellularRonsStep {
  IntracellularRonsState state;
  double h2o2_uptake = 0.0;  // alpha_ext C_H2O2 [mol/m3/s]
  double no2_uptake = 0.0;   // alpha_ext C_NO2  [mol/m3/s]
};

/// Advance both species with the same explicit-Euler uptake/clearance law.
inline IntracellularRonsStep AdvanceIntracellularRons(
    const IntracellularRonsState& current,
    double extracellular_h2o2,
    double extracellular_no2,
    double dt,
    const IntracellularRonsCoefficients& coefficients) {
  IntracellularRonsStep result;
  result.state.h2o2_load = std::max(0.0, current.h2o2_load);
  result.state.no2_load = std::max(0.0, current.no2_load);
  result.state.damage = std::max(0.0, current.damage);

  if (dt <= 0.0) {
    return result;
  }

  const double h2o2 = std::max(0.0, extracellular_h2o2);
  const double no2 = std::max(0.0, extracellular_no2);
  const double alpha_ext = std::max(0.0, coefficients.uptake_rate);
  const double beta_r = std::max(0.0, coefficients.clearance_rate);
  const double alpha_d =
      std::max(0.0, coefficients.damage_accumulation_rate);
  const double beta_d = std::max(0.0, coefficients.damage_repair_rate);

  result.h2o2_uptake = alpha_ext * h2o2;
  result.no2_uptake = alpha_ext * no2;

  // Forward Euler: every right-hand side uses the state at time n.
  result.state.h2o2_load = std::max(
      0.0, current.h2o2_load +
               dt * (result.h2o2_uptake - beta_r * current.h2o2_load));
  result.state.no2_load = std::max(
      0.0, current.no2_load +
               dt * (result.no2_uptake - beta_r * current.no2_load));

  const double total_load =
      std::max(0.0, current.h2o2_load) + std::max(0.0, current.no2_load);
  result.state.damage = std::max(
      0.0, current.damage +
               dt * (alpha_d * total_load - beta_d * current.damage));

  return result;
}

/// Probability that an agent is removed during one accepted step of length dt.
/// Damage below the death threshold gives probability 0. At and above that
/// threshold the hazard is constant, so crossing it does not delete the agent
/// immediately: each later step draws again.
inline double DeathRemovalProbability(double damage,
                                      double dt,
                                      double death_threshold,
                                      double removal_rate) {
  if (dt <= 0.0 || removal_rate <= 0.0 || damage < death_threshold) {
    return 0.0;
  }
  return 1.0 - std::exp(-removal_rate * dt);
}

}  // namespace bdm

#endif  // RONS_INTRACELLULAR_MODEL_H_
