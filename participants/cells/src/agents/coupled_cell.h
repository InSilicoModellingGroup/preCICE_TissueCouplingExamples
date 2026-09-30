// =============================================================================
// coupled_cell.h - Class declarations for MyCell
// Following ABM4bio pattern: declarations in .h, implementations in -inline.h
// =============================================================================
#ifndef MY_CELL_H_
#define MY_CELL_H_

#include "../global.h"
#include <algorithm>
#include <cstdint>
#include <limits>

namespace bdm {

/**
 * @brief Coupled cell class for OpenFOAM-preCICE-BioDynaMo simulations
 *
 * This class represents a particle/cell that can exchange data with
 * OpenFOAM via preCICE. Supports velocity, pressure, and temperature fields.
 */
class MyCell : public Cell {
    BDM_AGENT_HEADER(MyCell, Cell, 5);

    using VisualizationConfig = CavityConfig::VisualizationConfig;

private:
    // Physical properties from OpenFOAM
    double temperature_ = 300.0;
    double oxygen_ = 0.0;
    double oxygen_uptake_ = 0.0;
    double h2o2_uptake_ = 0.0;  // Agent H2O2 uptake rate [mol/m3/s]
    double h2o2_ = 0.0;   // Absolute H2O2 [mol/m3] for tissue_modelling
    double no2_ = 0.0;    // Absolute NO2- [mol/m3] for tissue_modelling
    double no2_uptake_ = 0.0;  // Agent NO2- uptake rate [mol/m3/s]
    double intracellular_h2o2_ = 0.0;  // Intracellular H2O2 load [mol/m3]
    double intracellular_no2_ = 0.0;   // Intracellular NO2- load [mol/m3]
    double oxidative_stress_ = 0.0;  // Combined intracellular RONS load [mol/m3]
    double damage_ = 0.0;            // Accumulated dimensionless damage
    double pressure_ = 0.0;
    Real3 velocity_ = {0.0, 0.0, 0.0};
    Real3 fluid_velocity_ = {0.0, 0.0, 0.0};
    Real3 fluid_acceleration_ = {0.0, 0.0, 0.0};
    Real3 fluid_vorticity_ = {0.0, 0.0, 0.0};

    // Turbulence fields from OpenFOAM
    double tke_ = 0.0;          // Turbulent kinetic energy (k)
    double epsilon_ = 0.0;      // Turbulent dissipation rate (epsilon)
    double omega_ = 0.0;        // Specific dissipation rate (omega)
    double nut_ = 0.0;          // Turbulent viscosity (nut)

    // Gradient information
    Real3 velocity_gradient_x_ = {0.0, 0.0, 0.0};
    Real3 velocity_gradient_y_ = {0.0, 0.0, 0.0};
    Real3 velocity_gradient_z_ = {0.0, 0.0, 0.0};
    double pressure_gradient_[3] = {0.0, 0.0, 0.0};
    double temperature_gradient_[3] = {0.0, 0.0, 0.0};

    // Coupling data
    Real3 openfoam_cell_center_ = {0.0, 0.0, 0.0};
    int precice_vertex_id_ = -1;

    // Mechanical properties
    Real3 accumulated_force_ = {0.0, 0.0, 0.0};
    double drag_Sp_ = 0.0;
    double added_mass_ = 0.0;
    Real3 post_division_repulsion_dir_ = {0.0, 0.0, 0.0};
    int post_division_repulsion_steps_remaining_ = 0;
    double post_division_repulsion_step_size_ = 0.0;
    Double3 cell_color_ = {0.5, 0.5, 0.5};
    double initial_diameter_ = 0.0;
    uint64_t track_id_ = 0;
    double residence_time_ = 0.0;  // Time elapsed since this agent was created [s]
    // phenotype_: 0=viable, 1=stressed, 2=dead (RONS damage)
    int phenotype_ = 0;  // 0 = viable / normoxic, 1 = stressed, 2 = dead

    // State flags
    bool is_stationary_ = true;
    bool data_updated_ = false;

public:
    // === Constructors ===
    MyCell();
    explicit MyCell(const Real3& openfoam_center);

    bool CanDivide() const;

    // === Velocity Interface ===
    void SetVelocity(const Real3& velocity);
    const Real3& GetVelocity() const;
    void SetFluidVelocity(const Real3& velocity);
    const Real3& GetFluidVelocity() const;
    void SetFluidAcceleration(const Real3& acceleration);
    const Real3& GetFluidAcceleration() const;
    void SetFluidVorticity(const Real3& vorticity);
    const Real3& GetFluidVorticity() const;
    void SetVelocityFromFluid();
    void SetVelocityGradient(const Real3& grad_x, const Real3& grad_y, const Real3& grad_z);
    std::array<Real3, 3> GetVelocityGradient() const;

    // === Pressure Interface ===
    void SetPressure(double pressure);
    double GetPressure() const;
    void SetPressureGradient(const std::array<double, 3>& gradient);
    std::array<double, 3> GetPressureGradient() const;

    // === Temperature Interface ===
    void SetTemperature(double temperature);
    double GetTemperature() const;
    void SetTemperatureGradient(const std::array<double, 3>& gradient);
    std::array<double, 3> GetTemperatureGradient() const;

    // === Oxygen Interface ===
    void SetOxygen(double oxygen);
    double GetOxygen() const;
    void SetOxygenUptake(double uptake);
    double GetOxygenUptake() const;
    void SetH2O2Uptake(double uptake);
    double GetH2O2Uptake() const;

    // === Tissue reactive species (H2O2 / NO2-) ===
    void SetH2O2(double concentration);
    double GetH2O2() const;
    void SetNO2(double concentration);
    double GetNO2() const;
    void SetNO2Uptake(double uptake);
    double GetNO2Uptake() const;
    void SetIntracellularH2O2(double concentration);
    double GetIntracellularH2O2() const;
    void SetIntracellularNO2(double concentration);
    double GetIntracellularNO2() const;
    void SetOxidativeStress(double stress);
    double GetOxidativeStress() const;
    void SetDamage(double damage);
    double GetDamage() const;

    // === Turbulence Interface ===
    void SetTKE(double k) { tke_ = k; }
    double GetTKE() const { return tke_; }
    void SetEpsilon(double eps) { epsilon_ = eps; }
    double GetEpsilon() const { return epsilon_; }
    void SetOmega(double om) { omega_ = om; }
    double GetOmega() const { return omega_; }
    void SetNut(double nu_t) { nut_ = nu_t; }
    double GetNut() const { return nut_; }

    // === Coupling Interface ===
    void SetPreciceVertexId(int id) { precice_vertex_id_ = id; }
    int GetPreciceVertexId() const { return precice_vertex_id_; }
    void SetOpenFoamCellCenter(const Real3& center) { openfoam_cell_center_ = center; }
    const Real3& GetOpenFoamCellCenter() const { return openfoam_cell_center_; }
    bool IsStationary() const { return is_stationary_; }
    void SetStationary(bool stationary) { is_stationary_ = stationary; }

    // === Force and Motion Interface ===
    void ResetForces();
    void SetDragCoefficient(double Sp);
    double GetDragCoefficient() const;
    void SetAddedMass(double added_mass);
    double GetAddedMass() const;
    void AddForce(const Real3& force);
    void SetAccumulatedForce(const Real3& force);
    const Real3& GetAccumulatedForce() const;
    void IntegrateForces(double dt);
    void StartPostDivisionRepulsion(const Real3& direction, int steps,
                                    double total_distance);
    void ApplyPostDivisionRepulsion();

    // === Compatibility Accessors ===
    int GetMeshVertexId() const { return GetPreciceVertexId(); }
    void SetMeshVertexId(int id) { SetPreciceVertexId(id); }

    // === Phenotype Interface ===
    int GetPhenotype() const { return phenotype_; }
    void SetPhenotype(int phenotype);

    // === Visualization ===
    void SetCellColor(const Double3& color) { cell_color_ = color; }
    const Double3& GetCellColor() const { return cell_color_; }
    uint64_t GetTrackId() const { return track_id_; }
    double GetResidenceTime() const { return residence_time_; }
    void SetResidenceTime(double residence_time) {
        residence_time_ = std::max(0.0, residence_time);
    }
    void AdvanceResidenceTime(double dt) {
        if (dt > 0.0) {
            residence_time_ += dt;
        }
    }

private:
    void InitializeCell();
    void UpdateSizeFromTemperature();
    void UpdateVisualization();
    void UpdateColorByVelocity(const VisualizationConfig& config);
    void UpdateColorByPressure(const VisualizationConfig& config);
    void UpdateColorByTemperature(const VisualizationConfig& config);
    void UpdateColorByOxygenPhenotype();
    void UpdateColorByH2O2();
    void UpdateColorByRonsFate();
};

// Type aliases for compatibility
typedef MyCell CavityCell;
typedef MyCell StationaryCell;
typedef MyCell DynamicCell;

} // namespace bdm

// Include inline implementations
#include "coupled_cell-inline.h"

#endif // MY_CELL_H_
