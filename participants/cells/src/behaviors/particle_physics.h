#ifndef OPENFOAM_INTERACTIONS_H_
#define OPENFOAM_INTERACTIONS_H_

#include "biodynamo.h"
#include "../agents/coupled_cell.h"
#include "../config/coupling_configuration.h"
#include <cmath>
#include <algorithm>
#include <array>

namespace bdm {

/**
 * @brief OpenFOAM-compatible particle interaction models
 *
 * This class implements particle-particle and particle-wall interactions
 * using the same models as OpenFOAM 2406 DEM (Discrete Element Method).
 * Cross-validated with OpenFOAM's collision models.
 */
class OpenFoamInteractions {
public:
    using CellType = MyCell;

    // Material properties - loaded from params.csv via global.h helpers
    struct MaterialProperties {
        double youngs_modulus;        // [Pa] - Particle stiffness
        double poisson_ratio;         // [-] - Poisson ratio
        double density;               // [kg/m^3] - Particle density (rho0)
        double friction_coefficient;  // [-] - Friction coefficient
        double restitution_coefficient; // [-] - Rebound coefficient
        double cohesion_energy;       // [J/m^3] - Cohesion energy density

        // Damping parameters
        double normal_damping;        // [-] - Normal damping factor
        double tangential_damping;    // [-] - Tangential damping factor

        // Default constructor loads from params.csv
        MaterialProperties()
            : youngs_modulus(YoungsModulus())
            , poisson_ratio(PoissonRatio())
            , density(GetParam<double>("agent.particle_density"))
            , friction_coefficient(FrictionCoefficient())
            , restitution_coefficient(RestitutionCoefficient())
            , cohesion_energy(GetParam<double>("interaction.cohesion_energy"))
            , normal_damping(GetParam<double>("interaction.normal_damping"))
            , tangential_damping(GetParam<double>("interaction.tangential_damping"))
        {}

        // Derived properties
        double shear_modulus() const {
            return youngs_modulus / (2.0 * (1.0 + poisson_ratio));
        }
    };

    struct WallProperties {
        double youngs_modulus;        // [Pa] - Wall stiffness
        double poisson_ratio;         // [-] - Wall Poisson ratio
        double friction_coefficient;  // [-] - Wall friction coefficient
        double restitution_coefficient; // [-] - Wall restitution coefficient

        // Default constructor loads from params.csv
        WallProperties()
            : youngs_modulus(GetParam<double>("interaction.wall_youngs_modulus"))
            , poisson_ratio(GetParam<double>("interaction.wall_poisson_ratio"))
            , friction_coefficient(GetParam<double>("interaction.wall_friction_coefficient"))
            , restitution_coefficient(WallRestitution())
        {}
    };

private:
    MaterialProperties particle_props_;
    WallProperties wall_props_;

    // Wall positions removed - now using face-based collision via preCICE mesh data

public:
    // Default constructor - all values loaded from params.csv
    OpenFoamInteractions()
        : particle_props_(), wall_props_() {
        // All defaults come from params.csv via constructors
    }

    OpenFoamInteractions(const MaterialProperties& particle_props,
                        const WallProperties& wall_props)
        : particle_props_(particle_props), wall_props_(wall_props) {
        // Wall collision now handled via preCICE face mesh
    }

    // === Particle-Particle Interactions ===

    void ComputeParticleParticleForce(CellType* agent1, CellType* agent2) const {

        const Real3 pos1 = agent1->GetPosition();
        const Real3 pos2 = agent2->GetPosition();
        const double r1 = agent1->GetDiameter() * 0.5;
        const double r2 = agent2->GetDiameter() * 0.5;

        // Relative position vector
        Real3 delta = {pos1[0] - pos2[0], pos1[1] - pos2[1], pos1[2] - pos2[2]};
        double distance = std::sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);

        // Contact detection
        double overlap = (r1 + r2) - distance;
        if (overlap <= 0.0) return;  // No contact

        // Unit normal vector (from agent2 to agent1)
        Real3 normal = {delta[0]/distance, delta[1]/distance, delta[2]/distance};

        // Relative velocity
        const Real3 vel1 = agent1->GetVelocity();
        const Real3 vel2 = agent2->GetVelocity();
        Real3 rel_vel = {vel1[0] - vel2[0], vel1[1] - vel2[1], vel1[2] - vel2[2]};

        // Normal and tangential relative velocities
        double v_normal = rel_vel[0]*normal[0] + rel_vel[1]*normal[1] + rel_vel[2]*normal[2];
        Real3 v_tangent = {
            rel_vel[0] - v_normal * normal[0],
            rel_vel[1] - v_normal * normal[1],
            rel_vel[2] - v_normal * normal[2]
        };

        // Effective properties (OpenFOAM DEM approach)
        double eff_radius = (r1 * r2) / (r1 + r2);
        double eff_mass = GetEffectiveMass(agent1, agent2);
        double eff_modulus = GetEffectiveModulus(particle_props_, particle_props_);

        // Normal force (Hertz-Mindlin model)
        double sqrt_overlap = std::sqrt(overlap * eff_radius);
        double normal_stiffness = (4.0/3.0) * eff_modulus * sqrt_overlap;
        double normal_force_mag = normal_stiffness * overlap;

        // Normal damping force (OpenFOAM style)
        double normal_damping_coeff = particle_props_.normal_damping *
                                     std::sqrt(normal_stiffness * eff_mass);
        double normal_damping_force = -normal_damping_coeff * v_normal;

        // Total normal force
        double total_normal_force = normal_force_mag + normal_damping_force;
        Real3 normal_force = {
            total_normal_force * normal[0],
            total_normal_force * normal[1],
            total_normal_force * normal[2]
        };

        // Tangential force (friction)
        Real3 tangent_force = {0.0, 0.0, 0.0};
        // Always enable particle-particle interactions in cavity case
        {
            tangent_force = ComputeFrictionForce(v_tangent, total_normal_force,
                                               particle_props_.friction_coefficient,
                                               particle_props_.tangential_damping,
                                               normal_stiffness, eff_mass);
        }

        // Total contact force
        Real3 total_force = {
            normal_force[0] + tangent_force[0],
            normal_force[1] + tangent_force[1],
            normal_force[2] + tangent_force[2]
        };

        // Apply forces (Newton's 3rd law)
        agent1->AddForce(total_force);
        agent2->AddForce({-total_force[0], -total_force[1], -total_force[2]});
    }

    // === Particle-Wall Interactions (Axis-Aligned Fallback) ===
    // Disabled for generic geometries. Wall collisions should be handled via
    // face-based collision (see KinematicCloudBehavior::ComputeFaceBasedWallForces)
    // using wall geometry and velocities provided by OpenFOAM/preCICE.
    void ComputeParticleWallForces(CellType* agent) const {
        (void)agent;
        Log::Warning("OpenFoamInteractions", "Axis-aligned wall fallback is disabled. "
                     "Provide wall faces via preCICE and use face-based collision.");
    }
    double GetEffectiveMass(const CellType* agent1, const CellType* agent2) const {
        double m1 = agent1->GetMass();
        double m2 = agent2->GetMass();
        return (m1 * m2) / (m1 + m2);
    }

    double GetEffectiveModulus(const MaterialProperties& props1,
                              const MaterialProperties& props2) const {
        double term1 = (1.0 - props1.poisson_ratio * props1.poisson_ratio) / props1.youngs_modulus;
        double term2 = (1.0 - props2.poisson_ratio * props2.poisson_ratio) / props2.youngs_modulus;
        return 1.0 / (term1 + term2);
    }

    Real3 ComputeFrictionForce(const Real3& tangent_velocity, double normal_force,
                              double friction_coeff, double damping_coeff,
                              double stiffness, double eff_mass) const {

        double v_tangent_mag = std::sqrt(tangent_velocity[0]*tangent_velocity[0] +
                                        tangent_velocity[1]*tangent_velocity[1] +
                                        tangent_velocity[2]*tangent_velocity[2]);

        if (v_tangent_mag < 1e-12) return {0.0, 0.0, 0.0};

        // Unit tangent vector
        Real3 tangent_unit = {
            tangent_velocity[0] / v_tangent_mag,
            tangent_velocity[1] / v_tangent_mag,
            tangent_velocity[2] / v_tangent_mag
        };

        // Maximum friction force (Coulomb criterion)
        double max_friction = friction_coeff * std::abs(normal_force);

        // Tangential stiffness (typically 2/7 of normal stiffness for spheres)
        double tangent_stiffness = (2.0/7.0) * stiffness;

        // Damping force in tangential direction
        double damping_force = damping_coeff * std::sqrt(tangent_stiffness * eff_mass) * v_tangent_mag;

        // Apply friction limit
        double friction_force = std::min(damping_force, max_friction);

        return {
            -friction_force * tangent_unit[0],
            -friction_force * tangent_unit[1],
            -friction_force * tangent_unit[2]
        };
    }

    // Axis-aligned collision helper is unused now that face-based collision is required.
    void CheckWallCollision(CellType*, const Real3&, double, int, double, const Real3&) const {}

    double GetEffectiveWallModulus() const {
        double particle_term = (1.0 - particle_props_.poisson_ratio * particle_props_.poisson_ratio) /
                              particle_props_.youngs_modulus;
        double wall_term = (1.0 - wall_props_.poisson_ratio * wall_props_.poisson_ratio) /
                          wall_props_.youngs_modulus;
        return 1.0 / (particle_term + wall_term);
    }
};

/**
 * @brief BioDynaMo behavior for applying OpenFOAM-style interactions
 */
class OpenFoamInteractionBehavior : public Behavior {
private:
    mutable OpenFoamInteractions interactions_;
    double* dt_ptr_;

public:
    explicit OpenFoamInteractionBehavior(double* dt_ptr) : dt_ptr_(dt_ptr) {}

    // Wall positions constructor removed - use face-based collision via preCICE

    OpenFoamInteractionBehavior(const OpenFoamInteractionBehavior& other)
        : Behavior(other), interactions_(other.interactions_), dt_ptr_(other.dt_ptr_) {}

    Behavior* New() const override {
        return new OpenFoamInteractionBehavior(dt_ptr_);
    }

    Behavior* NewCopy() const override {
        return new OpenFoamInteractionBehavior(*this);
    }

    void Run(Agent* agent) override {
        if (!dt_ptr_ || *dt_ptr_ <= 0.0) {
            return;
        }

        auto* cell = dynamic_cast<MyCell*>(agent);
        if (!cell) {
            return;
        }

        // Apply fluid drag force (OpenFOAM style coupling)
        ApplyFluidDragForce(cell);

        // Particle-particle interactions are handled by InteractionForce
        // Wall collisions must be handled via face-based collision (see KinematicCloudBehavior)
        // Integrate forces (now includes fluid drag and particle contacts)
        if (dt_ptr_ && *dt_ptr_ > 0.0) {
            cell->IntegrateForces(*dt_ptr_);
        }
    }
    // Apply fluid drag force based on OpenFOAM velocity field
    void ApplyFluidDragForce(MyCell* cell) {
        if (!dt_ptr_ || *dt_ptr_ <= 0.0) return;

        // Get fluid velocity at particle location (from OpenFOAM via preCICE)
        const Real3 fluid_velocity = cell->GetFluidVelocity();

        // Current particle velocity (kinematic state)
        const Real3 particle_velocity = cell->GetVelocity();

        // Relative velocity (fluid - particle)
        Real3 rel_velocity = {
            fluid_velocity[0] - particle_velocity[0],
            fluid_velocity[1] - particle_velocity[1],
            fluid_velocity[2] - particle_velocity[2]
        };

        double rel_speed_sq = rel_velocity[0]*rel_velocity[0] +
                              rel_velocity[1]*rel_velocity[1] +
                              rel_velocity[2]*rel_velocity[2];

        if (rel_speed_sq < 1e-16) {
            return; // No slip -> no drag
        }

        // Linear Stokes-like drag coefficient - viscosity provided by OpenFOAM via preCICE
        const double particle_diameter = cell->GetDiameter();
        const double fluid_viscosity = FluidViscosity();
        double drag_coeff = 3.0 * M_PI * fluid_viscosity * particle_diameter;

        Real3 drag_force = {
            drag_coeff * rel_velocity[0],
            drag_coeff * rel_velocity[1],
            drag_coeff * rel_velocity[2]
        };

        cell->AddForce(drag_force);
    }
};

/**
 * @brief Simple advection behavior that moves agents with the fluid velocity field
 */
class FluidAdvectionBehavior : public Behavior {
private:
    double* dt_ptr_;

public:
    explicit FluidAdvectionBehavior(double* dt_ptr) : dt_ptr_(dt_ptr) {}

    FluidAdvectionBehavior(const FluidAdvectionBehavior& other)
        : Behavior(other), dt_ptr_(other.dt_ptr_) {}

    Behavior* New() const override {
        return new FluidAdvectionBehavior(dt_ptr_);
    }

    Behavior* NewCopy() const override {
        return new FluidAdvectionBehavior(*this);
    }

    void Run(Agent* agent) override {
        if (!dt_ptr_ || *dt_ptr_ <= 0.0) {
            return;
        }

        auto* cell = dynamic_cast<MyCell*>(agent);
        if (!cell || cell->IsStationary()) {
            return;
        }

        const double dt = *dt_ptr_;
        const Real3 velocity = cell->GetVelocity();
        Real3 position = cell->GetPosition();

        position[0] += velocity[0] * dt;
        position[1] += velocity[1] * dt;
        position[2] += velocity[2] * dt;

        // Enforce 2D plane by keeping the original Z if necessary
        if (std::abs(position[2] - cell->GetOpenFoamCellCenter()[2]) > 1e-12) {
            position[2] = cell->GetOpenFoamCellCenter()[2];
        }

        cell->SetPosition(position);
    }
};

/**
 * @brief BioDynaMo InteractionForce for particle-particle interactions
 *
 * This class implements the InteractionForce interface to handle
 * particle-particle collisions using OpenFOAM-compatible models.
 */
class ParticleInteractionForce : public InteractionForce {
private:
    OpenFoamInteractions interactions_;

public:
    ParticleInteractionForce() = default;

    ParticleInteractionForce(const OpenFoamInteractions::MaterialProperties& particle_props,
                            const OpenFoamInteractions::WallProperties& wall_props)
        : interactions_(particle_props, wall_props) {}

    InteractionForce* NewCopy() const {
        return new ParticleInteractionForce(*this);
    }

    void ForceBetweenSpheres(const Agent* sphere_1, const Agent* sphere_2) const {
        // Cast to our cell types
        const auto* cell1 = dynamic_cast<const MyCell*>(sphere_1);
        const auto* cell2 = dynamic_cast<const MyCell*>(sphere_2);

        if (!cell1 || !cell2) return;

        // Use const_cast for force computation (BioDynaMo interface limitation)
        auto* nc_cell1 = const_cast<MyCell*>(cell1);
        auto* nc_cell2 = const_cast<MyCell*>(cell2);

        interactions_.ComputeParticleParticleForce(nc_cell1, nc_cell2);
    }
};

} // namespace bdm

#endif // OPENFOAM_INTERACTIONS_H_
