// =============================================================================
// coupled_cell-inline.h - Method implementations for MyCell
// Following ABM4bio pattern: declarations in .h, implementations in -inline.h
// =============================================================================
#ifndef MY_CELL_INLINE_H_
#define MY_CELL_INLINE_H_

#include <algorithm>
#include <atomic>

#include "coupled_cell.h"

namespace bdm {

// =============================================================================
// Constructors
// =============================================================================

inline MyCell::MyCell() : Base() {
    InitializeCell();
}

inline MyCell::MyCell(const Real3& openfoam_center) : Base(openfoam_center) {
    openfoam_cell_center_ = openfoam_center;
    InitializeCell();
}

inline bool MyCell::CanDivide() const {
    if (!GetParam<bool>("biology.enable_growth_division", false)) {
        return false;
    }

    const double threshold = GetParam<double>("biology.oxygen_growth_threshold", 0.25);
    return oxygen_ >= threshold;
}

// =============================================================================
// Velocity Methods
// =============================================================================

inline void MyCell::SetVelocity(const Real3& velocity) {
    velocity_ = velocity;
    data_updated_ = true;
    UpdateVisualization();
}

inline const Real3& MyCell::GetVelocity() const {
    return velocity_;
}

inline void MyCell::SetFluidVelocity(const Real3& velocity) {
    fluid_velocity_ = velocity;
    // 2D simulation: constrain Z component to zero
    if (IsForce2D()) {
        fluid_velocity_[2] = 0.0;
    }
}

inline const Real3& MyCell::GetFluidVelocity() const {
    return fluid_velocity_;
}

inline void MyCell::SetFluidAcceleration(const Real3& acceleration) {
    fluid_acceleration_ = acceleration;
    if (IsForce2D()) {
        fluid_acceleration_[2] = 0.0;
    }
}

inline const Real3& MyCell::GetFluidAcceleration() const {
    return fluid_acceleration_;
}

inline void MyCell::SetFluidVorticity(const Real3& vorticity) {
    fluid_vorticity_ = vorticity;
    if (IsForce2D()) {
        fluid_vorticity_[0] = 0.0;
        fluid_vorticity_[1] = 0.0;
    }
}

inline const Real3& MyCell::GetFluidVorticity() const {
    return fluid_vorticity_;
}

inline void MyCell::SetVelocityFromFluid() {
    velocity_ = fluid_velocity_;
}

inline void MyCell::SetVelocityGradient(const Real3& grad_x, const Real3& grad_y, const Real3& grad_z) {
    velocity_gradient_x_ = grad_x;
    velocity_gradient_y_ = grad_y;
    velocity_gradient_z_ = grad_z;
}

inline std::array<Real3, 3> MyCell::GetVelocityGradient() const {
    return {velocity_gradient_x_, velocity_gradient_y_, velocity_gradient_z_};
}

// =============================================================================
// Pressure Methods
// =============================================================================

inline void MyCell::SetPressure(double pressure) {
    pressure_ = pressure;
    data_updated_ = true;
    UpdateVisualization();
}

inline double MyCell::GetPressure() const {
    return pressure_;
}

inline void MyCell::SetPressureGradient(const std::array<double, 3>& gradient) {
    pressure_gradient_[0] = gradient[0];
    pressure_gradient_[1] = gradient[1];
    pressure_gradient_[2] = gradient[2];
}

inline std::array<double, 3> MyCell::GetPressureGradient() const {
    return {pressure_gradient_[0], pressure_gradient_[1], pressure_gradient_[2]};
}

// =============================================================================
// Temperature Methods
// =============================================================================

inline void MyCell::SetTemperature(double temperature) {
    // tissue_modelling: preCICE T channel carries absolute NO2- [mol/m3]
    if (GetParam<bool>("coupling.store_temperature_as_no2", false)) {
        SetNO2(temperature);
        return;
    }
    temperature_ = temperature;
    data_updated_ = true;
    if (!GetParam<bool>("simulation.enable_temperature_size_relaxation", false)) {
        UpdateSizeFromTemperature();
    }
    UpdateVisualization();
}

inline double MyCell::GetTemperature() const {
    return temperature_;
}

inline void MyCell::SetTemperatureGradient(const std::array<double, 3>& gradient) {
    temperature_gradient_[0] = gradient[0];
    temperature_gradient_[1] = gradient[1];
    temperature_gradient_[2] = gradient[2];
}

inline std::array<double, 3> MyCell::GetTemperatureGradient() const {
    return {temperature_gradient_[0], temperature_gradient_[1], temperature_gradient_[2]};
}

// =============================================================================
// Oxygen Methods
// =============================================================================

inline void MyCell::SetOxygen(double oxygen) {
    // Tissue cases: the preCICE H2O2 data channel is routed through this API.
    if (GetParam<bool>("coupling.store_oxygen_as_h2o2", false)) {
        SetH2O2(oxygen);
        return;
    }
    // Clamp to biology.oxygen_min/max so absolute concentration fields are not forced
    // into the legacy normalized [0,1] o2Cube range when that range is widened.
    const double o_min = GetParam<double>("biology.oxygen_min", 0.0);
    const double o_max = GetParam<double>("biology.oxygen_max", 1.0);
    if (o_max > o_min) {
        oxygen_ = std::clamp(oxygen, o_min, o_max);
    } else {
        oxygen_ = std::max(0.0, oxygen);
    }
}

inline double MyCell::GetOxygen() const {
    // When the oxygen API stores H2O2, sink feedback must read h2o2_.
    if (GetParam<bool>("coupling.store_oxygen_as_h2o2", false)) {
        return h2o2_;
    }
    return oxygen_;
}

inline void MyCell::SetOxygenUptake(double uptake) {
    oxygen_uptake_ = std::max(0.0, uptake);
}

inline double MyCell::GetOxygenUptake() const {
    return oxygen_uptake_;
}

inline void MyCell::SetH2O2Uptake(double uptake) {
    h2o2_uptake_ = std::max(0.0, uptake);
}

inline double MyCell::GetH2O2Uptake() const {
    return h2o2_uptake_;
}

inline void MyCell::SetH2O2(double concentration) {
    const double c_min = GetParam<double>("biology.h2o2_min", 0.0);
    const double c_max = GetParam<double>("biology.h2o2_max", 1000.0);
    if (c_max > c_min) {
        h2o2_ = std::clamp(concentration, c_min, c_max);
    } else {
        h2o2_ = std::max(0.0, concentration);
    }
    data_updated_ = true;
    UpdateVisualization();
}

inline double MyCell::GetH2O2() const {
    return h2o2_;
}

inline void MyCell::SetNO2(double concentration) {
    const double c_min = GetParam<double>("biology.no2_min", 0.0);
    const double c_max = GetParam<double>("biology.no2_max", 1000.0);
    if (c_max > c_min) {
        no2_ = std::clamp(concentration, c_min, c_max);
    } else {
        no2_ = std::max(0.0, concentration);
    }
    data_updated_ = true;
    UpdateVisualization();
}

inline double MyCell::GetNO2() const {
    return no2_;
}

inline void MyCell::SetNO2Uptake(double uptake) {
    no2_uptake_ = std::max(0.0, uptake);
}

inline double MyCell::GetNO2Uptake() const {
    return no2_uptake_;
}

inline void MyCell::SetIntracellularH2O2(double concentration) {
    intracellular_h2o2_ = std::max(0.0, concentration);
}

inline double MyCell::GetIntracellularH2O2() const {
    return intracellular_h2o2_;
}

inline void MyCell::SetIntracellularNO2(double concentration) {
    intracellular_no2_ = std::max(0.0, concentration);
}

inline double MyCell::GetIntracellularNO2() const {
    return intracellular_no2_;
}

inline void MyCell::SetOxidativeStress(double stress) {
    oxidative_stress_ = std::max(0.0, stress);
}

inline double MyCell::GetOxidativeStress() const {
    return oxidative_stress_;
}

inline void MyCell::SetDamage(double damage) {
    damage_ = std::max(0.0, damage);
}

inline double MyCell::GetDamage() const {
    return damage_;
}

// =============================================================================
// Force and Motion Methods
// =============================================================================

inline void MyCell::ResetForces() {
    accumulated_force_ = {0.0, 0.0, 0.0};
    drag_Sp_ = 0.0;
    added_mass_ = 0.0;
}

inline void MyCell::SetDragCoefficient(double Sp) {
    drag_Sp_ = std::max(0.0, Sp);
}

inline double MyCell::GetDragCoefficient() const {
    return drag_Sp_;
}

inline void MyCell::SetAddedMass(double added_mass) {
    added_mass_ = std::max(0.0, added_mass);
}

inline double MyCell::GetAddedMass() const {
    return added_mass_;
}

inline void MyCell::AddForce(const Real3& force) {
    accumulated_force_[0] += force[0];
    accumulated_force_[1] += force[1];
    accumulated_force_[2] += force[2];
}

inline void MyCell::SetAccumulatedForce(const Real3& force) {
    accumulated_force_ = force;
}

inline const Real3& MyCell::GetAccumulatedForce() const {
    return accumulated_force_;
}

inline void MyCell::IntegrateForces(double dt) {
    // Log first few calls to detect dt/position anomalies
    static std::atomic<int> log_calls{0};
    const int call_index = log_calls.load(std::memory_order_relaxed);
    const bool log_this = call_index < 6;
    if (log_this) {
        const auto& p0 = GetPosition();
        const auto& uf = fluid_velocity_;
        Log::Debug("IntegrateForces",
                          "call=", call_index,
                          " dt=", dt,
                          " mass=", GetMass(),
                          " drag_Sp=", drag_Sp_,
                          " force_2d=", IsForce2D(),
                          " fluid_vel=(", uf[0], ",", uf[1], ",", uf[2], ")",
                          " vel=(", velocity_[0], ",", velocity_[1], ",", velocity_[2], ")",
                          " pos=(", p0[0], ",", p0[1], ",", p0[2], ")");
    }

    if (!is_stationary_ && dt > 0.0) {
        if (IsMasslessTracerMode()) {
            // Massless tracers: explicit advection only (mapped U). Time refinement is
            // done in cells.h via simulation.biodynamo_timestep (scheduler steps per
            // preCICE window), not inner substeps here.
            Real3 pos = GetPosition();
            if (IsForce2D()) {
                pos[0] += velocity_[0] * dt;
                pos[1] += velocity_[1] * dt;
                pos[2] = openfoam_cell_center_[2];
                velocity_[2] = 0.0;
            } else {
                pos[0] += velocity_[0] * dt;
                pos[1] += velocity_[1] * dt;
                pos[2] += velocity_[2] * dt;
            }
            SetPosition(pos);
            if (log_this) {
                Log::Debug("IntegrateForces",
                                  "massless new_pos=(", pos[0], ",", pos[1], ",", pos[2], ")");
            }
        } else {
            const double mass = static_cast<double>(GetMass());
            const double effective_mass = std::max(mass + added_mass_, 1e-30);
            const Real3 U_old = velocity_;
            const bool force_2d = IsForce2D();

            const double b = drag_Sp_ * mass / effective_mass;

            if (b > 1e-12) {
                const double expFactor = std::exp(-b * dt);
                const double dtEff = (1.0 - expFactor) / b;
                velocity_[0] = fluid_velocity_[0] +
                               (U_old[0] - fluid_velocity_[0]) * expFactor +
                               accumulated_force_[0] * dtEff / effective_mass;
                velocity_[1] = fluid_velocity_[1] +
                               (U_old[1] - fluid_velocity_[1]) * expFactor +
                               accumulated_force_[1] * dtEff / effective_mass;

                if (force_2d) {
                    velocity_[2] = 0.0;
                } else {
                    velocity_[2] = fluid_velocity_[2] +
                                   (U_old[2] - fluid_velocity_[2]) * expFactor +
                                   accumulated_force_[2] * dtEff / effective_mass;
                }

                Real3 pos = GetPosition();
                pos[0] += 0.5 * (U_old[0] + velocity_[0]) * dt;
                pos[1] += 0.5 * (U_old[1] + velocity_[1]) * dt;
                if (force_2d) {
                    pos[2] = openfoam_cell_center_[2];
                } else {
                    pos[2] += 0.5 * (U_old[2] + velocity_[2]) * dt;
                }
                SetPosition(pos);
                if (log_this) {
                    const auto& uf = fluid_velocity_;
                    Log::Debug("IntegrateForces",
                                      "semi-implicit new_pos=(",
                                      pos[0], ",", pos[1], ",", pos[2],
                                      ") U_old=(", U_old[0], ",", U_old[1], ",", U_old[2],
                                      ") U_new=(", velocity_[0], ",", velocity_[1], ",", velocity_[2],
                                      ") U_fluid=(", uf[0], ",", uf[1], ",", uf[2], ") dt=", dt);
                }
            } else {
                const double scale = dt / effective_mass;
                velocity_[0] += accumulated_force_[0] * scale;
                velocity_[1] += accumulated_force_[1] * scale;
                if (force_2d) {
                    velocity_[2] = 0.0;
                } else {
                    velocity_[2] += accumulated_force_[2] * scale;
                }

                Real3 pos = GetPosition();
                pos[0] += 0.5 * (U_old[0] + velocity_[0]) * dt;
                pos[1] += 0.5 * (U_old[1] + velocity_[1]) * dt;
                if (force_2d) {
                    pos[2] = openfoam_cell_center_[2];
                } else {
                    pos[2] += 0.5 * (U_old[2] + velocity_[2]) * dt;
                }
                SetPosition(pos);
                if (log_this) {
                    const auto& uf = fluid_velocity_;
                    Log::Debug("IntegrateForces",
                                      "explicit new_pos=(",
                                      pos[0], ",", pos[1], ",", pos[2],
                                      ") U_old=(", U_old[0], ",", U_old[1], ",", U_old[2],
                                      ") U_new=(", velocity_[0], ",", velocity_[1], ",", velocity_[2],
                                      ") U_fluid=(", uf[0], ",", uf[1], ",", uf[2], ") dt=", dt);
                }
            }
        }
    }
    ResetForces();
    if (log_this) {
        log_calls.fetch_add(1, std::memory_order_relaxed);
    }
}

inline void MyCell::StartPostDivisionRepulsion(const Real3& direction, int steps,
                                               double total_distance) {
    post_division_repulsion_steps_remaining_ = 0;
    post_division_repulsion_step_size_ = 0.0;
    post_division_repulsion_dir_ = {0.0, 0.0, 0.0};

    if (steps <= 0 || total_distance <= 0.0) {
        return;
    }

    Real3 dir = direction;
    if (IsForce2D()) {
        dir[2] = 0.0;
    }

    const double norm = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
    if (norm <= 1e-12) {
        return;
    }

    post_division_repulsion_dir_[0] = dir[0] / norm;
    post_division_repulsion_dir_[1] = dir[1] / norm;
    post_division_repulsion_dir_[2] = dir[2] / norm;
    post_division_repulsion_steps_remaining_ = steps;
    post_division_repulsion_step_size_ = total_distance / static_cast<double>(steps);
}

inline void MyCell::ApplyPostDivisionRepulsion() {
    if (post_division_repulsion_steps_remaining_ <= 0 ||
        post_division_repulsion_step_size_ <= 0.0) {
        return;
    }

    Real3 pos = GetPosition();
    pos[0] += post_division_repulsion_dir_[0] * post_division_repulsion_step_size_;
    pos[1] += post_division_repulsion_dir_[1] * post_division_repulsion_step_size_;
    pos[2] += post_division_repulsion_dir_[2] * post_division_repulsion_step_size_;
    if (IsForce2D()) {
        pos[2] = openfoam_cell_center_[2];
    }
    SetPosition(pos);

    --post_division_repulsion_steps_remaining_;
    if (post_division_repulsion_steps_remaining_ <= 0) {
        post_division_repulsion_step_size_ = 0.0;
        post_division_repulsion_dir_ = {0.0, 0.0, 0.0};
    }
}

// =============================================================================
// Private Helpers
// =============================================================================

inline void MyCell::InitializeCell() {
    // Use diameter from params.csv via global.h helper
    initial_diameter_ = DefaultDiameter();
    SetDiameter(initial_diameter_);
    SetMass(DefaultMass());
    velocity_ = {0.0, 0.0, 0.0};
    fluid_velocity_ = {0.0, 0.0, 0.0};
    fluid_acceleration_ = {0.0, 0.0, 0.0};
    fluid_vorticity_ = {0.0, 0.0, 0.0};
    oxygen_ = 0.0;
    oxygen_uptake_ = 0.0;
    h2o2_uptake_ = 0.0;
    h2o2_ = 0.0;
    no2_ = 0.0;
    no2_uptake_ = 0.0;
    intracellular_h2o2_ = 0.0;
    intracellular_no2_ = 0.0;
    oxidative_stress_ = 0.0;
    damage_ = 0.0;
    track_id_ = GetUid().GetIndex();
    residence_time_ = 0.0;

    if (is_stationary_) {
        SetPosition(openfoam_cell_center_);
    }

    UpdateSizeFromTemperature();
}

inline void MyCell::UpdateSizeFromTemperature() {
    if (!GetParam<bool>("simulation.enable_temperature_size_scaling", false)) {
        return;
    }

    const double temp_min = GetParam<double>(
        "simulation.temperature_size_min",
        GetParam<double>("visualization.temperature_scale_min", 300.0));
    const double temp_max = GetParam<double>(
        "simulation.temperature_size_max",
        GetParam<double>("visualization.temperature_scale_max", 350.0));

    const double min_diameter = GetParam<double>(
        "simulation.temperature_size_min_diameter", initial_diameter_);
    const double max_diameter = GetParam<double>(
        "simulation.temperature_size_max_diameter", initial_diameter_ * 2.0);

    if (temp_max <= temp_min) {
        return;
    }

    const double lower_diameter = std::min(min_diameter, max_diameter);
    const double upper_diameter = std::max(min_diameter, max_diameter);

    double normalized = (temperature_ - temp_min) / (temp_max - temp_min);
    normalized = std::clamp(normalized, 0.0, 1.0);

    const double target_diameter =
        lower_diameter + normalized * (upper_diameter - lower_diameter);

    constexpr double kDiameterTolerance = 1e-9;
    if (std::abs(target_diameter - GetDiameter()) > kDiameterTolerance) {
        SetDiameter(target_diameter);
    }
}

inline void MyCell::UpdateVisualization() {
    const auto& config = PreciceVisualization();
    (void)config;
    const std::string color_mode =
        GetParam<std::string>("visualization.color_mode", "velocity");
    if (color_mode == "oxygen_phenotype" || color_mode == "phenotype") {
        UpdateColorByOxygenPhenotype();
    } else if (color_mode == "rons_fate" || color_mode == "fate") {
        UpdateColorByRonsFate();
    } else if (color_mode == "h2o2" || color_mode == "H2O2") {
        UpdateColorByH2O2();
    } else if (color_mode == "no2" || color_mode == "NO2" || color_mode == "temperature") {
        // temperature mode reused for NO2 when store_temperature_as_no2 is on
        if (GetParam<bool>("coupling.store_temperature_as_no2", false)) {
            cell_color_ = {no2_, no2_, no2_};
        } else {
            UpdateColorByTemperature(config);
        }
    } else {
        UpdateColorByVelocity(config);
    }
}

inline void MyCell::UpdateColorByVelocity(const VisualizationConfig& config) {
    (void)config;
    // Store raw velocity components as color to keep actual values (no normalization)
    cell_color_ = {velocity_[0], velocity_[1], velocity_[2]};
}

inline void MyCell::UpdateColorByPressure(const VisualizationConfig& config) {
    (void)config;
    cell_color_ = {pressure_, pressure_, pressure_};
}

inline void MyCell::UpdateColorByTemperature(const VisualizationConfig& config) {
    (void)config;
    cell_color_ = {temperature_, temperature_, temperature_};
}

inline void MyCell::SetPhenotype(int phenotype) {
    phenotype_ = phenotype;
    UpdateVisualization();
}

inline void MyCell::UpdateColorByOxygenPhenotype() {
    if (phenotype_ == 1) {
        // Necrotic: dark brown/red so cells are clearly distinct in ParaView
        cell_color_ = {0.55, 0.10, 0.05};
    } else {
        // Normoxic: green gradient — high O2 = bright green, low O2 = pale
        const double o2 = std::clamp(oxygen_, 0.0, 1.0);
        cell_color_ = {0.05 + 0.25 * (1.0 - o2), 0.40 + 0.55 * o2, 0.05};
    }
}

inline void MyCell::UpdateColorByH2O2() {
    const double c_max = std::max(GetParam<double>("biology.h2o2_max", 1000.0), 1e-12);
    const double n = std::clamp(h2o2_ / c_max, 0.0, 1.0);
    // Blue (low) -> cyan/green (high)
    cell_color_ = {0.05, 0.20 + 0.70 * n, 0.55 + 0.40 * n};
}

inline void MyCell::UpdateColorByRonsFate() {
    // 0 viable (green), 1 stressed (orange), 2 dead (dark red)
    if (phenotype_ >= 2) {
        cell_color_ = {0.45, 0.05, 0.05};
    } else if (phenotype_ == 1) {
        cell_color_ = {0.90, 0.45, 0.05};
    } else {
        cell_color_ = {0.10, 0.70, 0.20};
    }
}

} // namespace bdm

#endif // MY_CELL_INLINE_H_
