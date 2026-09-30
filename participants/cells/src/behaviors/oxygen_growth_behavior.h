#ifndef OXYGEN_GROWTH_BEHAVIOR_H_
#define OXYGEN_GROWTH_BEHAVIOR_H_

#include <algorithm>
#include <cmath>

#include "../global.h"
#include "../agents/coupled_cell.h"
#include "core/simulation.h"

namespace bdm {

class OxygenGrowthBehavior : public Behavior {
 private:
  double* dt_ptr_;
  bool enabled_;
  double oxygen_threshold_;
  double oxygen_quiescent_threshold_;
  double oxygen_min_;
  double oxygen_max_;
  double oxygen_sensitivity_power_;
  double growth_rate_;
  double growth_randomness_;
  double division_diameter_;
  double division_probability_base_;
  double division_probability_oxygen_power_;
  double division_probability_randomness_;
  double division_volume_ratio_;
  double daughter_offset_fraction_;
  bool keep_mother_position_on_division_;
  double post_division_repulsion_distance_fraction_;
  int post_division_repulsion_steps_;
  double oxygen_necrotic_threshold_;  // below this → necrotic phenotype
  int max_cell_count_;  // hard cap on total agents; 0 = unlimited

 public:
  explicit OxygenGrowthBehavior(double* dt_ptr)
      : dt_ptr_(dt_ptr),
        enabled_(GetParam<bool>("biology.enable_growth_division", false)),
        oxygen_threshold_(GetParam<double>("biology.oxygen_growth_threshold", 0.25)),
        oxygen_quiescent_threshold_(
            GetParam<double>("biology.oxygen_quiescent_threshold",
                             oxygen_threshold_)),
        oxygen_min_(GetParam<double>("biology.oxygen_min", 0.0)),
        oxygen_max_(GetParam<double>("biology.oxygen_max", 1.0)),
        oxygen_sensitivity_power_(
            GetParam<double>("biology.oxygen_sensitivity_power", 1.0)),
        growth_rate_(GetParam<double>("biology.growth_rate", 0.0)),
        growth_randomness_(GetParam<double>("biology.growth_randomness", 0.0)),
        division_diameter_(
            GetParam<double>("biology.division_diameter", 2.0 * DefaultDiameter())),
        division_probability_base_(
            GetParam<double>("biology.division_probability_base", 1.0)),
        division_probability_oxygen_power_(
            GetParam<double>("biology.division_probability_oxygen_power", 1.0)),
        division_probability_randomness_(
            GetParam<double>("biology.division_probability_randomness", 0.0)),
        division_volume_ratio_(GetParam<double>("biology.division_volume_ratio", 1.0)),
        daughter_offset_fraction_(
            GetParam<double>("biology.daughter_offset_fraction", 0.125)),
        keep_mother_position_on_division_(
            GetParam<bool>("biology.keep_mother_position_on_division", true)),
        post_division_repulsion_distance_fraction_(
            GetParam<double>("biology.post_division_repulsion_distance_fraction",
                             0.25)),
        post_division_repulsion_steps_(
            GetParam<int>("biology.post_division_repulsion_steps", 8)),
        oxygen_necrotic_threshold_(
            GetParam<double>("biology.oxygen_necrotic_threshold",
                             GetParam<double>("biology.oxygen_growth_threshold", 0.25))),
        max_cell_count_(GetParam<int>("biology.max_cell_count", 0)) {
    AlwaysCopyToNew();
  }

  OxygenGrowthBehavior(const OxygenGrowthBehavior& other)
      : Behavior(other),
        dt_ptr_(other.dt_ptr_),
        enabled_(other.enabled_),
        oxygen_threshold_(other.oxygen_threshold_),
        oxygen_quiescent_threshold_(other.oxygen_quiescent_threshold_),
        oxygen_min_(other.oxygen_min_),
        oxygen_max_(other.oxygen_max_),
        oxygen_sensitivity_power_(other.oxygen_sensitivity_power_),
        growth_rate_(other.growth_rate_),
        growth_randomness_(other.growth_randomness_),
        division_diameter_(other.division_diameter_),
        division_probability_base_(other.division_probability_base_),
        division_probability_oxygen_power_(other.division_probability_oxygen_power_),
        division_probability_randomness_(other.division_probability_randomness_),
        division_volume_ratio_(other.division_volume_ratio_),
        daughter_offset_fraction_(other.daughter_offset_fraction_),
        keep_mother_position_on_division_(
            other.keep_mother_position_on_division_),
        post_division_repulsion_distance_fraction_(
            other.post_division_repulsion_distance_fraction_),
        post_division_repulsion_steps_(other.post_division_repulsion_steps_),
        oxygen_necrotic_threshold_(other.oxygen_necrotic_threshold_),
        max_cell_count_(other.max_cell_count_) {}

  Behavior* New() const override { return new OxygenGrowthBehavior(dt_ptr_); }

  Behavior* NewCopy() const override { return new OxygenGrowthBehavior(*this); }

  void Run(Agent* agent) override {
    if (!enabled_ || !dt_ptr_ || *dt_ptr_ <= 0.0) {
      return;
    }

    auto* cell = dynamic_cast<MyCell*>(agent);
    if (!cell) {
      return;
    }

    // Continue short-lived mother/daughter separation after division.
    cell->ApplyPostDivisionRepulsion();

    auto* random = Simulation::GetActive()->GetRandom();
    const double oxygen = std::clamp(cell->GetOxygen(), oxygen_min_, oxygen_max_);

    // Phenotype transition: necrotic when below the necrotic threshold.
    if (oxygen < oxygen_necrotic_threshold_) {
      cell->SetPhenotype(1);  // necrotic
    } else {
      cell->SetPhenotype(0);  // normoxic
    }

    const double active_threshold =
        std::max(oxygen_threshold_, oxygen_quiescent_threshold_);
    if (oxygen < active_threshold) {
      return;
    }

    const double oxygen_span =
        std::max(oxygen_max_ - oxygen_quiescent_threshold_, 1e-12);
    const double oxygen_drive_linear =
        std::clamp((oxygen - oxygen_quiescent_threshold_) / oxygen_span, 0.0, 1.0);
    if (oxygen_drive_linear <= 0.0) {
      return;
    }
    const double oxygen_drive =
        std::pow(oxygen_drive_linear, std::max(oxygen_sensitivity_power_, 0.0));

    // Growth is deterministic every active timestep (like BioDynaMo's
    // ChangeVolume).  Only the *magnitude* has noise; the event itself
    // is not gated by a coin-flip so the population grows at a
    // predictable rate.  Division (below) is the stochastic event that
    // creates biological variability and de-synchronises the population.
    if (growth_rate_ > 0.0) {
      const double dt = *dt_ptr_;
      const double current_d = cell->GetDiameter();
      double growth_multiplier = oxygen_drive;
      if (growth_randomness_ > 0.0) {
        const double relative_noise = random->Uniform(-growth_randomness_,
                                                      growth_randomness_);
        growth_multiplier *= std::max(0.0, 1.0 + relative_noise);
      }
      const double new_d = current_d + growth_rate_ * growth_multiplier * dt;
      if (new_d > current_d) {
        cell->SetDiameter(new_d);

        const double rho_p = std::max(ParticleDensity(), 1e-12);
        const double pi = 3.14159265358979323846;
        const double volume = (pi / 6.0) * new_d * new_d * new_d;
        cell->SetMass(rho_p * volume);
      }
    }

    if (cell->GetDiameter() >= division_diameter_) {
      // Carrying-capacity guard: prevent division once global cell count hits the cap.
      if (max_cell_count_ > 0) {
        auto* rm = Simulation::GetActive()->GetResourceManager();
        if (static_cast<int>(rm->GetNumAgents()) >= max_cell_count_) {
          return;
        }
      }

      double division_probability =
          std::clamp(division_probability_base_, 0.0, 1.0) *
          std::pow(oxygen_drive_linear, std::max(division_probability_oxygen_power_, 0.0));
      if (division_probability_randomness_ > 0.0) {
        const double relative_noise = random->Uniform(-division_probability_randomness_,
                                                      division_probability_randomness_);
        division_probability *= std::max(0.0, 1.0 + relative_noise);
      }
      division_probability = std::clamp(division_probability, 0.0, 1.0);
      if (random->Uniform() > division_probability) {
        return;
      }

      const auto mother_position_before_division = cell->GetPosition();
      const auto mother_diameter_before_division = cell->GetDiameter();

      auto* daughter =
          bdm_static_cast<MyCell*>(cell->Divide(division_volume_ratio_));
      if (!daughter) {
        return;
      }

      // BioDynaMo initializes daughters with the agent defaults, so preserve
      // the mother's current coupled and intracellular RONS state.
      daughter->SetOxygen(cell->GetOxygen());
      daughter->SetNO2(cell->GetNO2());
      daughter->SetOxygenUptake(cell->GetOxygenUptake());
      daughter->SetH2O2Uptake(cell->GetH2O2Uptake());
      daughter->SetNO2Uptake(cell->GetNO2Uptake());
      daughter->SetIntracellularH2O2(cell->GetIntracellularH2O2());
      daughter->SetIntracellularNO2(cell->GetIntracellularNO2());
      daughter->SetOxidativeStress(cell->GetOxidativeStress());
      daughter->SetDamage(cell->GetDamage());
      daughter->SetPhenotype(cell->GetPhenotype());
      daughter->SetResidenceTime(0.0);

      // Keep daughter on the same OpenFOAM reference plane as the mother.
      // Without this, 2D z-clamping can place daughters on z=0.
      daughter->SetOpenFoamCellCenter(cell->GetOpenFoamCellCenter());

      if (keep_mother_position_on_division_) {
        cell->SetPosition(mother_position_before_division);
      }

      constexpr double kPi = 3.14159265358979323846;
      const double theta = random->Uniform(0.0, 2.0 * kPi);
      const double z = random->Uniform(-1.0, 1.0);
      const double r_xy = std::sqrt(std::max(0.0, 1.0 - z * z));
      Real3 direction = {std::cos(theta) * r_xy, std::sin(theta) * r_xy, z};
      if (IsForce2D()) {
        direction[2] = 0.0;
      }

      double norm = std::sqrt(direction[0] * direction[0] +
                              direction[1] * direction[1] +
                              direction[2] * direction[2]);
      if (norm < 1e-12) {
        direction = {1.0, 0.0, 0.0};
        norm = 1.0;
      }

      direction[0] /= norm;
      direction[1] /= norm;
      direction[2] /= norm;

      const double offset = std::max(0.0, daughter_offset_fraction_) *
                            std::max(mother_diameter_before_division, 1e-9);

      Real3 daughter_position = mother_position_before_division;
      if (keep_mother_position_on_division_) {
        daughter_position[0] += direction[0] * offset;
        daughter_position[1] += direction[1] * offset;
        daughter_position[2] += direction[2] * offset;
      } else {
        const double half_offset = 0.5 * offset;
        Real3 mother_position = mother_position_before_division;
        mother_position[0] -= direction[0] * half_offset;
        mother_position[1] -= direction[1] * half_offset;
        mother_position[2] -= direction[2] * half_offset;
        daughter_position[0] += direction[0] * half_offset;
        daughter_position[1] += direction[1] * half_offset;
        daughter_position[2] += direction[2] * half_offset;
        if (IsForce2D()) {
          mother_position[2] = mother_position_before_division[2];
        }
        cell->SetPosition(mother_position);
      }
      if (IsForce2D()) {
        daughter_position[2] = mother_position_before_division[2];
      }
      daughter->SetPosition(daughter_position);

      const double repulsion_distance =
          std::max(0.0, post_division_repulsion_distance_fraction_) *
          std::max(mother_diameter_before_division, 1e-9);
      if (post_division_repulsion_steps_ > 0 && repulsion_distance > 0.0) {
        Real3 mother_repulsion_dir = {-direction[0], -direction[1], -direction[2]};
        cell->StartPostDivisionRepulsion(mother_repulsion_dir,
                                         post_division_repulsion_steps_,
                                         repulsion_distance);
        daughter->StartPostDivisionRepulsion(direction, post_division_repulsion_steps_,
                                             repulsion_distance);
        // Apply the first push immediately in the same timestep as division.
        cell->ApplyPostDivisionRepulsion();
        daughter->ApplyPostDivisionRepulsion();
      }
    }
  }
};

}  // namespace bdm

#endif  // OXYGEN_GROWTH_BEHAVIOR_H_
