#ifndef TEMPERATURE_SIZE_RELAXATION_BEHAVIOR_H_
#define TEMPERATURE_SIZE_RELAXATION_BEHAVIOR_H_

#include <algorithm>
#include <cmath>

#include "../global.h"
#include "../agents/coupled_cell.h"
#include "core/simulation.h"

namespace bdm {

class TemperatureSizeRelaxationBehavior : public Behavior {
 private:
  double* dt_ptr_;
  bool enabled_;
  double tau_growth_;
  double tau_shrinkage_;
  bool use_volume_relaxation_;

  static double ClampTargetDiameter(double diameter, double min_diameter,
                                     double max_diameter) {
    const double lower = std::min(min_diameter, max_diameter);
    const double upper = std::max(min_diameter, max_diameter);
    return std::clamp(diameter, lower, upper);
  }

 public:
  explicit TemperatureSizeRelaxationBehavior(double* dt_ptr)
      : dt_ptr_(dt_ptr),
        enabled_(GetParam<bool>("simulation.enable_temperature_size_relaxation", false)),
        tau_growth_(GetParam<double>("simulation.temperature_size_relaxation_tau_growth", 1.0)),
        tau_shrinkage_(GetParam<double>("simulation.temperature_size_relaxation_tau_shrinkage", 1.0)),
        use_volume_relaxation_(GetParam<bool>("simulation.temperature_size_relaxation_use_volume", true)) {
    AlwaysCopyToNew();
  }

  Behavior* New() const override { return new TemperatureSizeRelaxationBehavior(dt_ptr_); }

  Behavior* NewCopy() const override { return new TemperatureSizeRelaxationBehavior(*this); }

  void Run(Agent* agent) override {
    if (!enabled_ || !dt_ptr_ || *dt_ptr_ <= 0.0) {
      return;
    }

    auto* cell = dynamic_cast<MyCell*>(agent);
    if (!cell) {
      return;
    }

    const double temperature = cell->GetTemperature();
    const double temp_min = GetParam<double>(
        "simulation.temperature_size_min",
        GetParam<double>("visualization.temperature_scale_min", 300.0));
    const double temp_max = GetParam<double>(
        "simulation.temperature_size_max",
        GetParam<double>("visualization.temperature_scale_max", 350.0));

    const double base_diameter = std::max(cell->GetDiameter(), DefaultDiameter());
    const double min_diameter = GetParam<double>(
        "simulation.temperature_size_min_diameter", base_diameter);
    const double max_diameter = GetParam<double>(
        "simulation.temperature_size_max_diameter", base_diameter * 2.0);

    if (temp_max <= temp_min) {
      return;
    }

    const double normalized = std::clamp((temperature - temp_min) / (temp_max - temp_min), 0.0, 1.0);
    const double lower_diameter = std::min(min_diameter, max_diameter);
    const double upper_diameter = std::max(min_diameter, max_diameter);
    const double target_diameter = lower_diameter + normalized * (upper_diameter - lower_diameter);

    const double current_diameter = cell->GetDiameter();
    const double tau = target_diameter >= current_diameter ? tau_growth_ : tau_shrinkage_;
    if (tau <= 0.0) {
      cell->SetDiameter(ClampTargetDiameter(target_diameter, min_diameter, max_diameter));
      return;
    }

    const double relaxation_factor = std::exp(-(*dt_ptr_) / tau);
    if (use_volume_relaxation_) {
      constexpr double kPi = 3.14159265358979323846;
      const double current_volume = (kPi / 6.0) * current_diameter * current_diameter * current_diameter;
      const double target_volume = (kPi / 6.0) * target_diameter * target_diameter * target_diameter;
      const double new_volume = target_volume + (current_volume - target_volume) * relaxation_factor;
      const double new_diameter = std::cbrt((6.0 * new_volume) / kPi);
      const double clamped_diameter = ClampTargetDiameter(new_diameter, min_diameter, max_diameter);
      cell->SetDiameter(clamped_diameter);

      const double rho_p = std::max(ParticleDensity(), 1e-12);
      const double volume_for_mass = (kPi / 6.0) * clamped_diameter * clamped_diameter * clamped_diameter;
      cell->SetMass(rho_p * volume_for_mass);
    } else {
      const double new_diameter = target_diameter + (current_diameter - target_diameter) * relaxation_factor;
      cell->SetDiameter(ClampTargetDiameter(new_diameter, min_diameter, max_diameter));
    }
  }
};

}  // namespace bdm

#endif  // TEMPERATURE_SIZE_RELAXATION_BEHAVIOR_H_
