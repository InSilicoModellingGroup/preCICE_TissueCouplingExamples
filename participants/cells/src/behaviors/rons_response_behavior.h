#ifndef RONS_RESPONSE_BEHAVIOR_H_
#define RONS_RESPONSE_BEHAVIOR_H_

#include "biodynamo.h"
#include "../agents/coupled_cell.h"
#include "../global.h"
#include "rons_intracellular_model.h"

namespace bdm {

/**
 * @brief Tissue RONS response: intracellular loading, damage, and cell fate.
 *
 * Used by tissue_43 / tissue_44. Operates on absolute H2O2 / NO2-
 * concentrations in SI units [mol/m3].
 * Case 4.3 enables the ABM4bio two-ODE model. H2O2 and NO2- then use the
 * same first-order uptake and intracellular-clearance coefficients.
 * Legacy behavior remains available for cases that do not enable that mode.
 *
 * Phenotype: 0=viable, 1=stressed, 2=dying.
 * Phenotype 2 marks an agent whose damage has reached the death threshold.
 * Removal is stochastic and happens later, after the biology step is accepted.
 */
class RonsResponseBehavior : public Behavior {
  BDM_BEHAVIOR_HEADER(RonsResponseBehavior, Behavior, 1);

 public:
  RonsResponseBehavior() = default;
  explicit RonsResponseBehavior(double* dt)
      : dt_ptr_(dt),
        enabled_(GetParam<bool>("biology.enable_rons_response", false)),
        intracellular_ode_enabled_(
            GetParam<bool>("biology.enable_intracellular_rons_ode", false)),
        rons_uptake_rate_(
            GetParam<double>("biology.rons_uptake_rate", 0.02)),
        rons_clearance_rate_(
            GetParam<double>("biology.rons_clearance_rate", 0.02)),
        damage_accumulation_rate_(GetParam<double>(
            "biology.damage_accumulation_rate",
            GetParam<double>("biology.damage_rate", 0.15))),
        h2o2_uptake_vmax_(GetParam<double>("biology.h2o2_uptake_vmax",
            GetParam<double>("coupling.oxygen_uptake_vmax", 0.005))),
        h2o2_uptake_km_(GetParam<double>("biology.h2o2_uptake_km",
            GetParam<double>("coupling.oxygen_uptake_km", 0.03))),
        no2_uptake_rate_(GetParam<double>("biology.no2_uptake_rate", 0.02)),
        stress_h2o2_weight_(GetParam<double>("biology.stress_h2o2_weight", 1.0)),
        stress_no2_weight_(GetParam<double>("biology.stress_no2_weight", 0.35)),
        stress_scale_(GetParam<double>(
            "biology.stress_scale_mol_m3",
            GetParam<double>("biology.stress_scale_uM", 0.05))),
        damage_rate_(GetParam<double>("biology.damage_rate", 0.15)),
        damage_repair_rate_(GetParam<double>("biology.damage_repair_rate", 0.02)),
        stress_threshold_(GetParam<double>("biology.rons_stress_threshold", 0.35)),
        death_threshold_(GetParam<double>("biology.rons_death_threshold", 1.0)) {}

  RonsResponseBehavior(const RonsResponseBehavior& other)
      : Behavior(other),
        dt_ptr_(other.dt_ptr_),
        enabled_(other.enabled_),
        intracellular_ode_enabled_(other.intracellular_ode_enabled_),
        rons_uptake_rate_(other.rons_uptake_rate_),
        rons_clearance_rate_(other.rons_clearance_rate_),
        damage_accumulation_rate_(other.damage_accumulation_rate_),
        h2o2_uptake_vmax_(other.h2o2_uptake_vmax_),
        h2o2_uptake_km_(other.h2o2_uptake_km_),
        no2_uptake_rate_(other.no2_uptake_rate_),
        stress_h2o2_weight_(other.stress_h2o2_weight_),
        stress_no2_weight_(other.stress_no2_weight_),
        stress_scale_(other.stress_scale_),
        damage_rate_(other.damage_rate_),
        damage_repair_rate_(other.damage_repair_rate_),
        stress_threshold_(other.stress_threshold_),
        death_threshold_(other.death_threshold_) {}

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
  }

  void Run(Agent* agent) override {
    if (!enabled_) return;
    auto* cell = dynamic_cast<MyCell*>(agent);
    if (!cell) return;
    if (cell->GetPhenotype() >= 2) return;  // already dead

    const double dt = (dt_ptr_ && *dt_ptr_ > 0.0) ? *dt_ptr_ : 0.0;
    if (dt <= 0.0) return;

    const double h2o2 = std::max(0.0, cell->GetH2O2());
    const double no2 = std::max(0.0, cell->GetNO2());

    if (intracellular_ode_enabled_) {
      const IntracellularRonsState current{
          cell->GetIntracellularH2O2(),
          cell->GetIntracellularNO2(),
          cell->GetDamage()};
      const IntracellularRonsCoefficients coefficients{
          rons_uptake_rate_,
          rons_clearance_rate_,
          damage_accumulation_rate_,
          damage_repair_rate_};
      const auto next =
          AdvanceIntracellularRons(current, h2o2, no2, dt, coefficients);

      // These rates are diagnostics in one-way case 4.3. Case 4.4 uses the
      // same first-order coefficient when assembling its OpenFOAM sinks.
      cell->SetH2O2Uptake(next.h2o2_uptake);
      cell->SetNO2Uptake(next.no2_uptake);
      cell->SetIntracellularH2O2(next.state.h2o2_load);
      cell->SetIntracellularNO2(next.state.no2_load);
      cell->SetOxidativeStress(
          next.state.h2o2_load + next.state.no2_load);
      cell->SetDamage(next.state.damage);

      int phenotype = 0;
      if (next.state.damage >= death_threshold_) {
        phenotype = 2;
      } else if (next.state.damage >= stress_threshold_) {
        phenotype = 1;
      }
      cell->SetPhenotype(phenotype);
      return;
    }

    // Legacy algebraic response retained for cases not yet migrated to the
    // intracellular ODE model.
    double h2o2_uptake = cell->GetOxygenUptake();
    if (h2o2_uptake <= 0.0 && h2o2_uptake_vmax_ > 0.0 && h2o2_uptake_km_ > 0.0) {
      h2o2_uptake = h2o2_uptake_vmax_ * h2o2 / (h2o2_uptake_km_ + h2o2);
      cell->SetOxygenUptake(h2o2_uptake);
    }
    // Case 4.4's adapter evaluates the same first-order kinetics for NO2Sink.
    (void)(no2_uptake_rate_ * no2);

    const double scale = std::max(stress_scale_, 1e-12);
    const double stress =
        stress_h2o2_weight_ * (h2o2 / scale) +
        stress_no2_weight_ * (no2 / scale);
    cell->SetOxidativeStress(stress);

    double damage = cell->GetDamage();
    damage += damage_rate_ * stress * dt;
    damage = std::max(0.0, damage - damage_repair_rate_ * dt);
    cell->SetDamage(damage);

    int phenotype = 0;
    if (damage >= death_threshold_) {
      phenotype = 2;
    } else if (damage >= stress_threshold_ || stress >= stress_threshold_) {
      phenotype = 1;
    }
    cell->SetPhenotype(phenotype);
  }

 private:
  double* dt_ptr_ = nullptr;
  bool enabled_ = false;
  bool intracellular_ode_enabled_ = false;
  double rons_uptake_rate_ = 0.02;
  double rons_clearance_rate_ = 0.02;
  double damage_accumulation_rate_ = 0.15;
  double h2o2_uptake_vmax_ = 0.005;
  double h2o2_uptake_km_ = 0.03;
  double no2_uptake_rate_ = 0.02;
  double stress_h2o2_weight_ = 1.0;
  double stress_no2_weight_ = 0.35;
  double stress_scale_ = 0.05;
  double damage_rate_ = 0.15;
  double damage_repair_rate_ = 0.02;
  double stress_threshold_ = 0.35;
  double death_threshold_ = 1.0;
};

}  // namespace bdm

#endif  // RONS_RESPONSE_BEHAVIOR_H_
