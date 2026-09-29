#include <cmath>
#include <cstdlib>
#include <iostream>

#include "rons_intracellular_model.h"

namespace {

bool Near(double a, double b, double tolerance = 1e-12) {
  return std::abs(a - b) <= tolerance;
}

void Require(bool condition, const char* message) {
  if (!condition) {
    std::cerr << "rons_intracellular_model_test failed: " << message << '\n';
    std::exit(1);
  }
}

}  // namespace

int main() {
  const bdm::IntracellularRonsCoefficients coefficients{
      0.02,  // alpha_ext [1/s]
      0.02,  // beta_R [1/s]
      0.15,  // alpha_D [m3/(mol s)]
      0.02   // beta_D [1/s]
  };

  bdm::IntracellularRonsState state;
  auto step =
      bdm::AdvanceIntracellularRons(state, 0.1, 0.1, 0.1, coefficients);

  // Equal extracellular concentrations and coefficients must produce exactly
  // equal uptake and intracellular loads for H2O2 and NO2-.
  Require(Near(step.h2o2_uptake, 0.002), "unexpected H2O2 uptake");
  Require(Near(step.no2_uptake, 0.002), "unexpected NO2 uptake");
  Require(Near(step.state.h2o2_load, 0.0002),
          "unexpected first-step H2O2 load");
  Require(Near(step.state.no2_load, 0.0002),
          "unexpected first-step NO2 load");
  Require(Near(step.state.damage, 0.0), "damage must use the old load");

  state = step.state;
  step = bdm::AdvanceIntracellularRons(state, 0.1, 0.1, 0.1, coefficients);
  Require(Near(step.state.h2o2_load, 0.0003996),
          "unexpected second-step H2O2 load");
  Require(Near(step.state.no2_load, 0.0003996),
          "unexpected second-step NO2 load");
  Require(Near(step.state.damage, 0.000006),
          "unexpected second-step damage");

  // With exposure removed, intracellular load and damage clear
  // proportionally to their current values.
  const auto cleared =
      bdm::AdvanceIntracellularRons(step.state, 0.0, 0.0, 0.1, coefficients);
  Require(cleared.state.h2o2_load < step.state.h2o2_load,
          "H2O2 load did not clear");
  Require(cleared.state.no2_load < step.state.no2_load,
          "NO2 load did not clear");
  Require(cleared.state.damage >= 0.0, "damage became negative");

  // Sustained equal exposure must preserve species symmetry and eventually
  // cross the configured case-4.3 death threshold D=1.
  state = {};
  for (int i = 0; i < 6000; ++i) {
    state =
        bdm::AdvanceIntracellularRons(state, 0.1, 0.1, 0.1, coefficients)
            .state;
  }
  Require(Near(state.h2o2_load, state.no2_load),
          "equal species diverged");
  Require(state.damage >= 1.0,
          "sustained exposure did not reach the death threshold");

  Require(Near(bdm::DeathRemovalProbability(0.99, 0.1, 1.0, 0.02), 0.0),
          "damage below the death threshold must not be removed");
  Require(Near(bdm::DeathRemovalProbability(1.0, 0.0, 1.0, 0.02), 0.0),
          "a zero timestep must not remove an agent");
  const double at_threshold =
      bdm::DeathRemovalProbability(1.0, 0.1, 1.0, 0.02);
  Require(Near(at_threshold, 1.0 - std::exp(-0.002), 1e-12),
          "unexpected removal probability at the death threshold");
  Require(Near(bdm::DeathRemovalProbability(5.0, 0.1, 1.0, 0.02),
               at_threshold),
          "hazard must stay capped above the death threshold");

  std::cout << "rons_intracellular_model_test passed\n";
  return 0;
}
