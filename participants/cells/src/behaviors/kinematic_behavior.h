#ifndef BEHAVIORS_H_
#define BEHAVIORS_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <vector>

#include "../global.h"
#include "../agents/coupled_cell.h"
#include "particle_physics.h"
#include "../environment/precice_mesh_mapping.h"
#include "boundary_geometry.h"

namespace bdm {

/**
 * @brief Unified behavior combining OpenFOAM kinematic cloud dynamics with wall interactions
 *
 * This behavior replaces separate advection and drag behaviors with a single,
 * comprehensive implementation that handles:
 * 1. Fluid drag forces from OpenFOAM velocity field
 * 2. Wall collision forces using OpenFOAM DEM models
 * 3. Force integration and particle motion
 */
/**
 * @brief Wall face structure for face-based collision detection
 *
 * Mirrors bdm_adapter::PreciceAdapter::WallFace for use in behaviors.
 * Contains wall face geometry from OpenFOAM via preCICE WallSurfaceMesh.
 *
 * Wall velocity is automatically obtained from OpenFOAM boundary conditions:
 * - fixedValue: Uses the specified velocity (e.g., moving walls)
 * - noSlip: Zero velocity (stationary walls)
 * - Other BCs: Extracted from the computed boundary field
 */
struct WallFaceData {
    std::array<double, 3> center;   // Face center position
    std::array<double, 3> normal;   // OpenFOAM patch.nf() unit normal, outward from fluid domain
    double area;                     // Face area
    std::array<double, 3> velocity; // Face velocity (for moving walls)

    // Check if this is a moving wall
    bool IsMoving() const {
        return (std::abs(velocity[0]) > 1e-12 ||
                std::abs(velocity[1]) > 1e-12 ||
                std::abs(velocity[2]) > 1e-12);
    }
};

class KinematicCloudBehavior : public Behavior {
private:
    mutable OpenFoamInteractions interactions_;
    double* dt_ptr_;
    std::vector<WallFaceData> wall_faces_;  // Wall faces from OpenFOAM
    bool use_face_based_collision_ = false;  // True if wall faces available
    std::shared_ptr<const bdm_adapter::BoundaryGeometry> boundary_geometry_;

public:
    explicit KinematicCloudBehavior(double* dt_ptr) : dt_ptr_(dt_ptr) {}

    KinematicCloudBehavior(const KinematicCloudBehavior& other)
        : Behavior(other), interactions_(other.interactions_), dt_ptr_(other.dt_ptr_),
          wall_faces_(other.wall_faces_), use_face_based_collision_(other.use_face_based_collision_),
          boundary_geometry_(other.boundary_geometry_) {}

    Behavior* New() const override {
        return new KinematicCloudBehavior(dt_ptr_);
    }

    Behavior* NewCopy() const override {
        return new KinematicCloudBehavior(*this);
    }

    void Run(Agent* agent) override {
        if (!dt_ptr_ || *dt_ptr_ <= 0.0) {
            return;
        }

        auto* cell = dynamic_cast<MyCell*>(agent);
        if (!cell) {
            return;
        }

        // Keep the pre-integration position so collision handling can
        // detect swept wall crossings in this sub-step (tunneling guard).
        const Real3 pre_integrate_position = cell->GetPosition();

        // Apply fluid velocity (either as force or direct assignment depending on mode)
        // Uses params.csv: coupling.enable_drag_physics, simulation.use_massless_tracers
        if (IsDragPhysicsEnabled() || IsMasslessTracerMode()) {
            ApplyFluidDragForce(cell);
        }

        // Apply pressure gradient force (important when ρ_p ≈ ρ_f)
        // Uses params.csv: coupling.enable_pressure_gradient
        if (IsPressureGradientEnabled() && !IsMasslessTracerMode()) {
            ApplyPressureGradientForce(cell);
        }

        // Apply gravity force: F_g = m * g
        // Uses params.csv: coupling.enable_gravity, environment.gravity_*
        if (IsGravityEnabled() && !IsMasslessTracerMode()) {
            ApplyGravityForce(cell);
        }

        // Apply virtual/added mass force (important when ρ_p ≈ ρ_f)
        // Uses params.csv: coupling.enable_virtual_mass
        if (IsVirtualMassEnabled() && !IsMasslessTracerMode()) {
            ApplyVirtualMassForce(cell);
        }

        // Apply lift forces (Saffman shear-lift, Magnus rotation-lift)
        // Uses params.csv: coupling.enable_lift_forces
        if (IsLiftForcesEnabled() && !IsMasslessTracerMode()) {
            ApplyLiftForces(cell);
        }

        // Apply turbulent dispersion force (stochastic)
        // Uses params.csv: coupling.enable_turbulence_dispersion
        if (IsTurbulenceDispersionEnabled() && !IsMasslessTracerMode()) {
            ApplyTurbulentDispersionForce(cell);
        }

        // Integrate forces (for massless tracers, this just updates position from velocity)
        if (dt_ptr_ && *dt_ptr_ > 0.0) {
            cell->IntegrateForces(*dt_ptr_);
        }

        if (IsWallConfinementEnabled() && boundary_geometry_ &&
            !boundary_geometry_->Empty()) {
            ResolveBoundaryPath(cell, pre_integrate_position, *dt_ptr_);
        } else if (IsWallConfinementEnabled() &&
                   use_face_based_collision_ && !wall_faces_.empty()) {
            // Compatibility fallback for old configurations. Complex cases
            // fail fast at startup unless exact triangles are available.
            ComputeFaceBasedWallForces(cell, &pre_integrate_position);
        }
    }

    // Wall bounds are now set via SetWallFaces() using preCICE mesh data
    // No hardcoded axis-aligned walls - all wall geometry comes from OpenFOAM

    /**
     * @brief Set wall faces from OpenFOAM via preCICE for face-based collision
     *
     * This enables accurate wall collision detection for any geometry
     * using actual wall face centers, normals, and velocities from
     * OpenFOAM's WallSurfaceMesh exported via preCICE.
     *
     * IMPORTANT: This is required for moving walls. The wall velocity is
     * automatically obtained from OpenFOAM boundary conditions (e.g., fixedValue U).
     * Without this data, the fallback assumes all walls are stationary.
     *
     * @tparam WallFaceType Type with center, normal, area, velocity members
     * @param faces Vector of wall faces from adapter.GetWallFaces()
     */
    template<typename WallFaceType>
    void SetWallFaces(const std::vector<WallFaceType>& faces) {
        wall_faces_.clear();
        wall_faces_.reserve(faces.size());
        size_t moving_count = 0;
        size_t dropped_count = 0;
        for (const auto& face : faces) {
            const double area = face.area;
            const double nx = face.normal[0];
            const double ny = face.normal[1];
            const double nz = face.normal[2];
            const double nmag = std::sqrt(nx * nx + ny * ny + nz * nz);
            if (!std::isfinite(face.center[0]) ||
                !std::isfinite(face.center[1]) ||
                !std::isfinite(face.center[2]) ||
                !std::isfinite(area) ||
                area <= 0.0 ||
                !std::isfinite(nmag) ||
                nmag <= 1e-12) {
                ++dropped_count;
                continue;
            }

            WallFaceData wf;
            wf.center = face.center;
            wf.normal = {nx / nmag, ny / nmag, nz / nmag};
            wf.area = area;
            wf.velocity = face.velocity;
            wall_faces_.push_back(wf);
            if (wf.IsMoving()) ++moving_count;
        }
        use_face_based_collision_ = !wall_faces_.empty();

        if (!use_face_based_collision_) {
            Log::Warning("KinematicCloudBehavior",
                              "No valid wall faces available (received ",
                              faces.size(), ", dropped ", dropped_count,
                              "). Face-based wall collisions are disabled.");
        } else if (moving_count > 0) {
            Log::Info("KinematicCloudBehavior", "Set ", wall_faces_.size(),
                          " wall faces (", moving_count, " moving, dropped ",
                          dropped_count, ")");
        } else {
            Log::Info("KinematicCloudBehavior", "Set ", wall_faces_.size(),
                          " wall faces (all stationary, dropped ",
                          dropped_count, ")");
        }
    }

    /**
     * @brief Check if face-based collision is enabled
     */
    bool IsFaceBasedCollisionEnabled() const {
        return use_face_based_collision_ && !wall_faces_.empty();
    }

    /**
     * @brief Get number of wall faces
     */
    size_t GetWallFaceCount() const {
        return wall_faces_.size();
    }

    void SetBoundaryGeometry(
        std::shared_ptr<const bdm_adapter::BoundaryGeometry> geometry) {
        boundary_geometry_ = std::move(geometry);
    }

    bool HasBoundaryGeometry() const {
        return boundary_geometry_ && !boundary_geometry_->Empty();
    }

private:
    void ResolveBoundaryPath(MyCell* cell, const Real3& start, double dt) {
        if (!boundary_geometry_) return;
        const double radius = std::max(0.0, cell->GetDiameter() * 0.5);
        const double slop = std::max(1e-10, radius * 1e-6);
        bdm_adapter::BoundaryVec3 segmentStart = {
            start[0], start[1], start[2]};
        const Real3 proposed = cell->GetPosition();
        bdm_adapter::BoundaryVec3 segmentEnd = {
            proposed[0], proposed[1], proposed[2]};
        Real3 velocity = cell->GetVelocity();

        const auto initialWall =
            boundary_geometry_->ClosestPoint(segmentStart, true);
        if (initialWall.found) {
            const double penetration = radius
                + bdm_adapter::BDot(
                    bdm_adapter::BSub(segmentStart, initialWall.point),
                    initialWall.normal);
            if (penetration > 0.0) {
                segmentStart = bdm_adapter::BSub(
                    segmentStart,
                    bdm_adapter::BScale(
                        initialWall.normal, penetration + slop));
            }
        }

        const int maxContacts = std::max(
            1, GetParam<int>("interaction.wall_max_contacts", 8));
        for (int iteration = 0; iteration < maxContacts; ++iteration) {
            const auto contact = boundary_geometry_->FindFirstWallContact(
                segmentStart, segmentEnd, radius, slop * 0.1);
            if (!contact.found) break;

            const auto displacement =
                bdm_adapter::BSub(segmentEnd, segmentStart);
            const double outwardPath =
                bdm_adapter::BDot(displacement, contact.normal);
            if (contact.fraction <= 1e-9 && outwardPath <= 0.0) {
                break;
            }

            const auto hitCenter = bdm_adapter::BAdd(
                segmentStart,
                bdm_adapter::BScale(displacement, contact.fraction));
            const auto insideCenter = bdm_adapter::BSub(
                contact.point,
                bdm_adapter::BScale(contact.normal, radius + slop));

            bdm_adapter::BoundaryVec3 remaining =
                bdm_adapter::BSub(segmentEnd, hitCenter);
            bdm_adapter::BoundaryVec3 relativeVelocity = {
                velocity[0] - contact.velocity[0],
                velocity[1] - contact.velocity[1],
                velocity[2] - contact.velocity[2]};
            bdm_adapter::BoundaryVelocityMode responseMode =
                bdm_adapter::BoundaryVelocityMode::kRebound;
            switch (GetWallVelocityResponse()) {
                case WallVelocityResponse::kNone:
                    responseMode = bdm_adapter::BoundaryVelocityMode::kNone;
                    break;
                case WallVelocityResponse::kSlip:
                    responseMode = bdm_adapter::BoundaryVelocityMode::kSlip;
                    break;
                case WallVelocityResponse::kStick:
                    responseMode = bdm_adapter::BoundaryVelocityMode::kStick;
                    break;
                case WallVelocityResponse::kRebound:
                    responseMode = bdm_adapter::BoundaryVelocityMode::kRebound;
                    break;
            }
            const auto response = bdm_adapter::ApplyBoundaryResponse(
                responseMode, relativeVelocity, remaining, contact.normal,
                contact.velocity,
                std::max(0.0, 1.0 - contact.fraction) * dt,
                GetParam<double>("interaction.wall_restitution", 0.8),
                GetParam<double>(
                    "interaction.wall_friction_coefficient", 0.0));
            relativeVelocity = response.relativeVelocity;
            remaining = response.remainingDisplacement;

            velocity = {
                contact.velocity[0] + relativeVelocity[0],
                contact.velocity[1] + relativeVelocity[1],
                contact.velocity[2] + relativeVelocity[2]};
            segmentStart = insideCenter;
            segmentEnd = bdm_adapter::BAdd(insideCenter, remaining);
        }

        Real3 corrected = {segmentEnd[0], segmentEnd[1], segmentEnd[2]};
        if (IsForce2D()) {
            corrected[2] = DomainZPlane();
            velocity[2] = 0.0;
        }
        cell->SetPosition(corrected);
        cell->SetVelocity(velocity);
    }

    // Apply fluid drag force based on OpenFOAM velocity field
    void ApplyFluidDragForce(MyCell* cell) {
        if (!dt_ptr_ || *dt_ptr_ <= 0.0) return;

        // Get fluid velocity at particle location (from OpenFOAM via preCICE)
        const Real3 fluid_velocity = cell->GetFluidVelocity();

        // Massless tracer mode: particles move exactly with fluid (params.csv)
        if (IsMasslessTracerMode()) {
            cell->SetVelocity(fluid_velocity);
            return;
        }

        // Check if drag physics is enabled (params.csv)
        if (!IsDragPhysicsEnabled()) {
            return;
        }

        // Standard drag physics mode (default behavior)
        // Uses semi-implicit/analytical integration (OpenFOAM style)
        // Instead of adding drag force explicitly, we compute the implicit coefficient Sp
        // and let IntegrateForces() handle the analytical integration

        const Real3 particle_velocity = cell->GetVelocity();

        // Fluid properties from params.csv
        const double rho_f = FluidDensity();      // [kg/m^3]
        const double mu_f  = FluidViscosity();    // [Pa.s]

        const double d = cell->GetDiameter();
        const double m = cell->GetMass();

        if (d > 0.0 && m > 0.0) {
            const double pi = 3.14159265358979323846;
            const double vol = (pi / 6.0) * d * d * d;
            const double rho_p = m / vol;

            // Compute relative velocity magnitude for Reynolds number
            const bool force_2d = IsForce2D();
            Real3 u_rel = {
                fluid_velocity[0] - particle_velocity[0],
                fluid_velocity[1] - particle_velocity[1],
                force_2d ? 0.0 : (fluid_velocity[2] - particle_velocity[2])
            };
            const double mag_u_rel = std::sqrt(u_rel[0]*u_rel[0] + u_rel[1]*u_rel[1] + u_rel[2]*u_rel[2]);

            if (rho_p > 0.0) {
                // Compute Reynolds number
                const double Re = (mag_u_rel > 1e-12) ?
                    rho_f * mag_u_rel * d / mu_f : 0.0;

                // Schiller-Naumann drag correlation (same as OpenFOAM SphereDragForce)
                // CdRe = Cd * Re, where Cd is the drag coefficient
                double CdRe;
                if (Re > 1000.0) {
                    // Newton regime: Cd ≈ 0.424
                    CdRe = 0.424 * Re;
                } else {
                    // Stokes to intermediate regime
                    // Cd = 24/Re * (1 + Re^(2/3)/6)
                    CdRe = 24.0 * (1.0 + (1.0/6.0) * std::pow(Re, 2.0/3.0));
                }

                // Implicit drag coefficient Sp [1/s]
                // From OpenFOAM: Sp = (3/4) * (mu/rho_p) * CdRe / d^2
                // This appears in the equation: dU/dt = Sp * (Uc - Up)
                const double Sp = 0.75 * mu_f * CdRe / (rho_p * d * d);

                // Store the implicit coefficient for semi-implicit integration
                // The integration will solve: dU/dt = Sp*(Uc - U) + F_explicit/m
                cell->SetDragCoefficient(Sp);

                // Note: Do NOT add drag as explicit force - it's handled implicitly
                // Only non-drag forces (gravity, wall collision, etc.) should use AddForce()
                return;
            }
        }

        // Fallback: no drag if particle properties invalid
        cell->SetDragCoefficient(0.0);

    }

    /**
     * @brief Apply pressure gradient force (OpenFOAM PressureGradientForce equivalent)
     *
     * The pressure gradient force accounts for the fluid pressure acting on the particle surface.
     * For a particle immersed in a pressure gradient, this force equals the pressure gradient
     * times the particle volume (buoyancy-like effect).
     *
     * F_p = -V_p * ∇p
     *
     * where:
     *   V_p = particle volume = (π/6) * d³
     *   ∇p = pressure gradient vector
     *
     * This force is important when particle density is similar to fluid density (ρ_p ≈ ρ_f).
     * It's typically small for dense particles (ρ_p >> ρ_f) but significant for:
     *   - Neutrally buoyant particles
     *   - Bubbles
     *   - Particles in high-pressure gradient regions
     */
    void ApplyPressureGradientForce(MyCell* cell) {
        // Check params.csv settings
        if (!IsPressureGradientEnabled()) return;
        if (IsMasslessTracerMode()) return;

        const double d = cell->GetDiameter();
        if (d <= 0.0) return;

        // Particle volume
        const double pi = 3.14159265358979323846;
        const double V_p = (pi / 6.0) * d * d * d;
        const double rho_f = FluidDensity();
        const Real3& fluid_acc = cell->GetFluidAcceleration();
        const bool force_2d = IsForce2D();

        Real3 pressure_force = {
            rho_f * V_p * fluid_acc[0],
            rho_f * V_p * fluid_acc[1],
            force_2d ? 0.0 : rho_f * V_p * fluid_acc[2]
        };

        // Only add force if it's significant
        const double mag_force = std::sqrt(pressure_force[0]*pressure_force[0] +
                                           pressure_force[1]*pressure_force[1] +
                                           pressure_force[2]*pressure_force[2]);
        if (mag_force > 1e-20) {
            cell->AddForce(pressure_force);
        }
    }

    /**
     * @brief Apply gravity force
     *
     * Simple body force: F_g = m * g
     *
     * where:
     *   m = particle mass
     *   g = gravity vector from params.csv (environment.gravity_x/y/z)
     *
     * For buoyancy, the net gravity force is:
     *   F_net = (rho_p - rho_f) * V_p * g
     *
     * This is equivalent to gravity + pressure gradient force in hydrostatic case.
     */
    void ApplyGravityForce(MyCell* cell) {
        const double m = cell->GetMass();
        if (m <= 0.0) return;

        const auto g_arr = GravityVector();
        const Real3 g = {g_arr[0], g_arr[1], g_arr[2]};
        const bool force_2d = IsForce2D();

        Real3 gravity_force = {
            m * g[0],
            m * g[1],
            force_2d ? 0.0 : m * g[2]
        };

        cell->AddForce(gravity_force);
    }

    /**
     * @brief Apply virtual (added) mass force
     *
     * The virtual mass force accounts for the inertia of the fluid displaced
     * by the accelerating particle. Important when ρ_p ≈ ρ_f.
     *
     * F_vm = C_vm * V_p * ρ_f * (DU_f/Dt - dU_p/dt)
     *
     * where:
     *   C_vm = virtual mass coefficient (0.5 for sphere)
     *   V_p = particle volume
     *   ρ_f = fluid density
     *   DU_f/Dt = material derivative of fluid velocity (includes convective term)
     *   dU_p/dt = particle acceleration
     *
     * Simplified formulation (ignoring convective derivative):
     *   F_vm ≈ C_vm * V_p * ρ_f * (a_f - a_p)
     *
     * Note: This implementation uses the simplified form. For accurate results
     * in highly unsteady flows, the full material derivative should be used.
     */
    void ApplyVirtualMassForce(MyCell* cell) {
        const double d = cell->GetDiameter();
        const double m = cell->GetMass();
        if (d <= 0.0 || m <= 0.0) return;

        const double pi = 3.14159265358979323846;
        const double V_p = (pi / 6.0) * d * d * d;
        const double rho_f = FluidDensity();
        const double C_vm = 0.5;  // Virtual mass coefficient for sphere
        const Real3& fluid_acc = cell->GetFluidAcceleration();
        const bool force_2d = IsForce2D();

        cell->SetAddedMass(C_vm * rho_f * V_p);

        Real3 vm_force = {
            C_vm * V_p * rho_f * fluid_acc[0],
            C_vm * V_p * rho_f * fluid_acc[1],
            force_2d ? 0.0 : C_vm * V_p * rho_f * fluid_acc[2]
        };

        // Only add if significant
        const double mag = std::sqrt(vm_force[0]*vm_force[0] + vm_force[1]*vm_force[1] + vm_force[2]*vm_force[2]);
        if (mag > 1e-20) {
            cell->AddForce(vm_force);
        }
    }

    /**
     * @brief Apply Saffman-Mei-style lift using OpenFOAM vorticity.
     *
     * Active model:
     *   u_rel = u_f - v_p
     *   Re    = rho_f |u_rel| d / mu_f
     *   Rew   = rho_f |omega_f| d^2 / mu_f
     *   beta  = 0.5 Rew/Re
     *   alpha = 0.3314 sqrt(beta)
     *   f     = (1 - alpha) exp(-0.1 Re) + alpha
     *   Cld   = 6.46 f                                  for Re < 40
     *         = 6.46 * 0.0524 * sqrt(beta Re)            for Re >= 40
     *   Cl    = 3 Cld / (2 pi sqrt(Rew))
     *   F_L   = rho_f V_p Cl (u_rel x omega_f)
     *
     * omega_f is the OpenFOAM fvc::curl(U) field transferred as curlUc.
     */
    void ApplyLiftForces(MyCell* cell) {
        const double d = cell->GetDiameter();
        if (d <= 0.0) return;

        const Real3 fluid_vel = cell->GetFluidVelocity();
        const Real3 particle_vel = cell->GetVelocity();
        const Real3 fluid_vorticity = cell->GetFluidVorticity();
        const bool force_2d = IsForce2D();

        // Relative velocity
        Real3 u_rel = {
            fluid_vel[0] - particle_vel[0],
            fluid_vel[1] - particle_vel[1],
            force_2d ? 0.0 : fluid_vel[2] - particle_vel[2]
        };

        const double mag_u_rel = std::sqrt(u_rel[0]*u_rel[0] + u_rel[1]*u_rel[1] + u_rel[2]*u_rel[2]);
        if (mag_u_rel < 1e-12) return;

        const double rho_f = FluidDensity();
        const double mu_f = FluidViscosity();
        const double pi = 3.14159265358979323846;
        const double V_p = (pi / 6.0) * d * d * d;
        const double vort_mag = std::sqrt(fluid_vorticity[0]*fluid_vorticity[0] +
                                          fluid_vorticity[1]*fluid_vorticity[1] +
                                          fluid_vorticity[2]*fluid_vorticity[2]);
        if (rho_f <= 0.0 || mu_f <= 0.0 || vort_mag < 1e-12) return;

        const double Re = rho_f * mag_u_rel * d / std::max(mu_f, 1e-30);
        const double Rew = rho_f * vort_mag * d * d / std::max(mu_f, 1e-30);
        const double beta = 0.5 * (Rew / std::max(Re, 1e-30));
        const double alpha = 0.3314 * std::sqrt(std::max(beta, 0.0));
        const double f = (1.0 - alpha) * std::exp(-0.1 * Re) + alpha;
        const double Cld = (Re < 40.0)
            ? 6.46 * f
            : 6.46 * 0.0524 * std::sqrt(std::max(beta * Re, 0.0));
        const double Cl = 3.0 / (2.0 * pi * std::sqrt(std::max(Rew, 1e-30))) * Cld;

        Real3 cross = {
            u_rel[1] * fluid_vorticity[2] - u_rel[2] * fluid_vorticity[1],
            u_rel[2] * fluid_vorticity[0] - u_rel[0] * fluid_vorticity[2],
            u_rel[0] * fluid_vorticity[1] - u_rel[1] * fluid_vorticity[0]
        };

        Real3 lift_force = {
            rho_f * V_p * Cl * cross[0],
            rho_f * V_p * Cl * cross[1],
            force_2d ? 0.0 : rho_f * V_p * Cl * cross[2]
        };

        cell->AddForce(lift_force);
    }

    /**
     * @brief Apply turbulent dispersion force
     *
     * Models the effect of turbulent velocity fluctuations on particle motion.
     * Uses a stochastic approach where the turbulent force is modeled as:
     *
     * F_turb = β * m * (u' / τ_t)
     *
     * where:
     *   β = random coefficient ~ N(0,1)
     *   u' = turbulent velocity fluctuation ~ sqrt(2k/3) for isotropic turbulence
     *   k = turbulent kinetic energy
     *   τ_t = turbulent time scale ~ k/ε or L/u'
     *
     * Since we don't have k-ε from OpenFOAM (icoFoam is laminar), this uses a
     * simplified model based on local flow conditions:
     *   u' ~ 0.1 * |U| (turbulence intensity)
     *   τ_t ~ d / |u_rel| (particle crossing time)
     */
    void ApplyTurbulentDispersionForce(MyCell* cell) {
        const double d = cell->GetDiameter();
        const double m = cell->GetMass();
        if (d <= 0.0 || m <= 0.0) return;

        const Real3 fluid_vel = cell->GetFluidVelocity();
        const Real3 particle_vel = cell->GetVelocity();
        const bool force_2d = IsForce2D();

        const double U_mag = std::sqrt(fluid_vel[0]*fluid_vel[0] + fluid_vel[1]*fluid_vel[1] + fluid_vel[2]*fluid_vel[2]);
        if (U_mag < 1e-12) return;

        // Relative velocity magnitude
        Real3 u_rel = {
            fluid_vel[0] - particle_vel[0],
            fluid_vel[1] - particle_vel[1],
            force_2d ? 0.0 : fluid_vel[2] - particle_vel[2]
        };
        const double mag_u_rel = std::sqrt(u_rel[0]*u_rel[0] + u_rel[1]*u_rel[1] + u_rel[2]*u_rel[2]);

        // Turbulence parameters (simplified model)
        const double turb_intensity = 0.1;  // 10% turbulence intensity
        const double u_prime = turb_intensity * U_mag;

        // Turbulent time scale (particle crossing time)
        double tau_t = 0.01;  // Default time scale
        if (mag_u_rel > 1e-12) {
            tau_t = d / mag_u_rel;
        }
        tau_t = std::max(tau_t, 1e-6);  // Prevent division issues

        // Generate random perturbation using linear congruential generator
        // Simple but fast random number generation for stochastic forces
        static thread_local unsigned int turb_seed = 12345 + static_cast<unsigned int>(reinterpret_cast<uintptr_t>(cell) & 0xFFFF);

        auto rand_uniform = []() -> double {
            turb_seed = turb_seed * 1103515245 + 12345;
            return static_cast<double>(turb_seed & 0x7FFFFFFF) / 0x7FFFFFFF;
        };

        // Box-Muller transform for normal distribution N(0,1)
        double u1 = std::max(rand_uniform(), 1e-12);
        double u2 = rand_uniform();
        double n1 = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * 3.14159265358979 * u2);
        double n2 = std::sqrt(-2.0 * std::log(u1)) * std::sin(2.0 * 3.14159265358979 * u2);
        u1 = std::max(rand_uniform(), 1e-12);
        u2 = rand_uniform();
        double n3 = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * 3.14159265358979 * u2);

        // Stochastic force components using standard normal random numbers
        Real3 turb_force = {
            m * u_prime / tau_t * n1,
            m * u_prime / tau_t * n2,
            force_2d ? 0.0 : m * u_prime / tau_t * n3
        };

        cell->AddForce(turb_force);
    }

    /**
     * @brief Face-based wall collision detection using actual wall geometry from OpenFOAM
     *
     * Uses wall face centers, normals, and velocities received from OpenFOAM via preCICE
     * WallSurfaceMesh. Works for any geometry (complex shapes, not just axis-aligned boxes).
     *
     * WALL VELOCITY HANDLING:
     * - Wall velocities are automatically obtained from OpenFOAM boundary conditions
     * - Moving walls (fixedValue BC): Uses the BC velocity for relative velocity calculation
     * - Stationary walls (noSlip BC): Zero velocity
     * - No user configuration needed - data comes directly from OpenFOAM via preCICE
     *
     * For each wall face:
     * 1. Compute signed distance from particle to wall plane
     * 2. If particle overlaps wall (distance < radius), apply collision response
     * 3. Use relative velocity (v_particle - v_wall) for correct collision physics
     * 4. Use wall normal for accurate force direction
     */
    void ComputeFaceBasedWallForces(MyCell* cell, const Real3* previous_pos = nullptr) {
        const double radius = cell->GetDiameter() * 0.5;
        if (radius <= 0.0) {
            return;
        }

        const double wall_restitution = GetParam<double>("interaction.wall_restitution", 0.8);
        const double wall_friction = GetParam<double>("interaction.wall_friction_coefficient", 0.0);
        const double position_slop = std::max(1e-9, radius * 1e-6);

        for (int iter = 0; iter < 4; ++iter) {
            const Real3 pos = cell->GetPosition();
            const Real3 vel = cell->GetVelocity();
            const WallFaceData* best_face = nullptr;
            Real3 best_normal = {0.0, 0.0, 0.0};
            double best_penetration = 0.0;

            for (const auto& face : wall_faces_) {
                Real3 nhat = {face.normal[0], face.normal[1], face.normal[2]};
                const double nmag = std::sqrt(nhat[0]*nhat[0] + nhat[1]*nhat[1] + nhat[2]*nhat[2]);
                if (nmag <= 1e-12) {
                    continue;
                }
                nhat[0] /= nmag;
                nhat[1] /= nmag;
                nhat[2] /= nmag;

                const double dx = pos[0] - face.center[0];
                const double dy = pos[1] - face.center[1];
                const double dz = pos[2] - face.center[2];
                const double signed_dist = dx * nhat[0] + dy * nhat[1] + dz * nhat[2];
                const double penetration = radius + signed_dist;
                const bool local_overlap = (penetration > 0.0 && signed_dist <= radius);

                // Continuous collision detection: catch crossings where the
                // particle swept through the wall plane during this sub-step.
                bool crossed = false;
                double lateral_dx = dx;
                double lateral_dy = dy;
                double lateral_dz = dz;
                if (previous_pos) {
                    const double pdx = (*previous_pos)[0] - face.center[0];
                    const double pdy = (*previous_pos)[1] - face.center[1];
                    const double pdz = (*previous_pos)[2] - face.center[2];
                    const double prev_signed_dist =
                        pdx * nhat[0] + pdy * nhat[1] + pdz * nhat[2];
                    crossed = (prev_signed_dist < -radius && signed_dist >= -radius);

                    if (crossed) {
                        const double denom = signed_dist - prev_signed_dist;
                        double t_hit = (std::abs(denom) > 1e-14)
                                           ? ((-radius - prev_signed_dist) / denom)
                                           : 1.0;
                        t_hit = std::clamp(t_hit, 0.0, 1.0);
                        lateral_dx = pdx + t_hit * (dx - pdx);
                        lateral_dy = pdy + t_hit * (dy - pdy);
                        lateral_dz = pdz + t_hit * (dz - pdz);
                    }
                }

                if (!local_overlap && !crossed) {
                    continue;
                }

                double lateral_signed_dist =
                    lateral_dx * nhat[0] + lateral_dy * nhat[1] + lateral_dz * nhat[2];
                double lateral_dist_sq =
                    lateral_dx * lateral_dx +
                    lateral_dy * lateral_dy +
                    lateral_dz * lateral_dz -
                    lateral_signed_dist * lateral_signed_dist;
                lateral_dist_sq = std::max(0.0, lateral_dist_sq);
                const double face_size = std::sqrt(std::max(face.area, 0.0));
                double lateral_limit = radius + 2.5 * face_size;
                lateral_limit = std::max(lateral_limit, radius * 1.25);
                if (lateral_dist_sq > lateral_limit * lateral_limit) {
                    continue;
                }

                if (!best_face || penetration > best_penetration) {
                    best_face = &face;
                    best_normal = nhat;
                    best_penetration = penetration;
                }
            }

            if (!best_face) {
                break;
            }

            const auto& wall_vel = best_face->velocity;
            Real3 v_rel = {
                vel[0] - wall_vel[0],
                vel[1] - wall_vel[1],
                vel[2] - wall_vel[2]
            };
            const double v_rel_normal =
                v_rel[0] * best_normal[0] +
                v_rel[1] * best_normal[1] +
                v_rel[2] * best_normal[2];

            if (v_rel_normal > 0.0) {
                Real3 v_rel_t = {
                    v_rel[0] - v_rel_normal * best_normal[0],
                    v_rel[1] - v_rel_normal * best_normal[1],
                    v_rel[2] - v_rel_normal * best_normal[2]
                };
                Real3 new_vel = {
                    wall_vel[0] + v_rel[0] - (1.0 + wall_restitution) * v_rel_normal * best_normal[0] - wall_friction * v_rel_t[0],
                    wall_vel[1] + v_rel[1] - (1.0 + wall_restitution) * v_rel_normal * best_normal[1] - wall_friction * v_rel_t[1],
                    wall_vel[2] + v_rel[2] - (1.0 + wall_restitution) * v_rel_normal * best_normal[2] - wall_friction * v_rel_t[2]
                };
                if (IsForce2D()) {
                    new_vel[2] = 0.0;
                }
                cell->SetVelocity(new_vel);
            }

            Real3 new_pos = {
                pos[0] - (best_penetration + position_slop) * best_normal[0],
                pos[1] - (best_penetration + position_slop) * best_normal[1],
                pos[2] - (best_penetration + position_slop) * best_normal[2]
            };
            if (IsForce2D()) {
                new_pos[2] = DomainZPlane();
            }
            cell->SetPosition(new_pos);
        }
    }

};

/**
 * @brief Legacy behaviors for backward compatibility
 */

// Automatic Fluid-Particle Coupling Behavior
// Uses fluid velocity transferred from OpenFOAM to advect particles
class AutomaticFluidParticleCoupling : public Behavior {
 public:
  AutomaticFluidParticleCoupling(double* dt, const bool* enable_motion,
                                 const bool* is_2d = nullptr,
                                 const double* z_plane = nullptr)
      : dt_(dt),
        enable_motion_(enable_motion),
        is_2d_(is_2d),
        z_plane_(z_plane) {}
  AutomaticFluidParticleCoupling(const AutomaticFluidParticleCoupling& other)
      : Behavior(other),
        dt_(other.dt_),
        enable_motion_(other.enable_motion_),
        is_2d_(other.is_2d_),
        z_plane_(other.z_plane_) {}

  Behavior* New() const override {
    return new AutomaticFluidParticleCoupling(dt_, enable_motion_, is_2d_,
                                              z_plane_);
  }
  Behavior* NewCopy() const override {
    return new AutomaticFluidParticleCoupling(*this);
  }

  void Run(Agent* a) override {
    auto* cell = dynamic_cast<MyCell*>(a);
    if (!cell || dt_ == nullptr) {
      return;
    }
    if (enable_motion_ != nullptr && !(*enable_motion_)) {
      return;
    }
    const auto& velocity = cell->GetVelocity();
    auto new_pos = cell->GetPosition();
    double dt = *dt_;
    new_pos[0] += velocity[0] * dt;
    new_pos[1] += velocity[1] * dt;
    if (is_2d_ != nullptr && *is_2d_) {
      if (z_plane_ != nullptr) {
        new_pos[2] = *z_plane_;
      }
    } else {
      new_pos[2] += velocity[2] * dt;
    }
    cell->SetPosition(new_pos);
  }

 private:
  double* dt_;
  const bool* enable_motion_;
  const bool* is_2d_;
  const double* z_plane_;
};

// Temperature behavior - adjust size based on temperature
class TemperatureBehavior : public Behavior {
 public:
  explicit TemperatureBehavior(double factor = 0.0) : factor_(factor) {}
  TemperatureBehavior(const TemperatureBehavior& other)
      : Behavior(other), factor_(other.factor_) {}

  Behavior* New() const override {
    return new TemperatureBehavior(factor_);
  }
  Behavior* NewCopy() const override {
    return new TemperatureBehavior(*this);
  }

  void Run(Agent* a) override {
    auto* cell = dynamic_cast<MyCell*>(a);
    if (!cell) return;
    const double T = cell->GetTemperature();
    (void)T;       // Temperature behaviour no longer alters size
    (void)factor_; // Growth disabled; retain member to keep interface stable
  }

 private:
  double factor_;  // growth factor per temperature unit
};

// Placeholder pressure behavior
class PressureBehavior : public Behavior {
 public:
  PressureBehavior() = default;
  PressureBehavior(const PressureBehavior& other) : Behavior(other) {}

  Behavior* New() const override { return new PressureBehavior(); }
  Behavior* NewCopy() const override {
    return new PressureBehavior(*this);
  }

  void Run(Agent* a) override {
    (void)a; // no-op placeholder
  }
};

}  // namespace bdm

#endif  // BEHAVIORS_H_
