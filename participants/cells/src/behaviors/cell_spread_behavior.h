#ifndef CELL_SPREAD_BEHAVIOR_H_
#define CELL_SPREAD_BEHAVIOR_H_

#include <algorithm>
#include <cmath>

#include "../global.h"
#include "../agents/coupled_cell.h"
#include "core/simulation.h"

namespace bdm {

/**
 * @brief Cell-cell contact overlap/adhesion behavior.
 *
 * Runs every ABM step.  For each cell it iterates over all other cells,
 * detects overlaps (centre-to-centre distance < sum of radii), applies a
 * position-level correction that pushes the two centres apart, and can apply
 * a weak short-range adhesion correction for close non-overlapping cells.
 *
 * Because the correction is applied directly to SetPosition() (like
 * ApplyPostDivisionRepulsion), it works even when
 * simulation.enable_particle_motion = false (is_stationary_ = true).
 *
 * Parameters (params.csv):
 *   biology.cell_spread_overlap_strength  [0-1]  fraction of overlap resolved
 *                                                 per step (default 0.5).
 *   biology.cell_spread_search_radius_factor     neighbourhood cutoff as a
 *                                                 multiple of this cell's
 *                                                 diameter (default 2.0).
 *   biology.cell_spread_adhesion_strength        weak attraction strength for
 *                                                 close non-overlapping cells.
 *   biology.cell_spread_adhesion_range_factor    adhesion cutoff as a multiple
 *                                                 of the sum of radii.
 *   biology.cell_spread_max_displacement_fraction
 *                                                 cap per-step correction as a
 *                                                 fraction of this cell radius.
 *
 * Thread-safety note: each Run() only writes to its own agent via
 * SetPosition(); reads from other agents are read-only.  This is safe with
 * BioDynaMo's parallel behaviour scheduler.
 */
class CellCellSpreadBehavior : public Behavior {
 private:
 double spread_strength_;
  double adhesion_strength_;
  double adhesion_range_factor_;
  double max_displacement_fraction_;
  double search_radius_factor_;
  // Run every N BDM steps (1 = every step, 5 = every 5th step, etc.).
  // All cells see the same global step counter so the behaviour either
  // fires for all cells or skips for all cells in a given step — keeps
  // the O(N²) cost predictable and tunable via params.csv.
  int run_interval_;

 public:
  CellCellSpreadBehavior()
      : spread_strength_(
            GetParam<double>("biology.cell_spread_overlap_strength", 0.5)),
        adhesion_strength_(
            GetParam<double>("biology.cell_spread_adhesion_strength", 0.0)),
        adhesion_range_factor_(
            GetParam<double>("biology.cell_spread_adhesion_range_factor", 1.2)),
        max_displacement_fraction_(
            GetParam<double>("biology.cell_spread_max_displacement_fraction", 0.25)),
        search_radius_factor_(
            GetParam<double>("biology.cell_spread_search_radius_factor", 2.0)),
        run_interval_(
            GetParam<int>("biology.cell_spread_run_interval", 1)) {
    AlwaysCopyToNew();  // daughters inherit this behavior automatically
  }

  CellCellSpreadBehavior(const CellCellSpreadBehavior& other)
      : Behavior(other),
        spread_strength_(other.spread_strength_),
        adhesion_strength_(other.adhesion_strength_),
        adhesion_range_factor_(other.adhesion_range_factor_),
        max_displacement_fraction_(other.max_displacement_fraction_),
        search_radius_factor_(other.search_radius_factor_),
        run_interval_(other.run_interval_) {}

  Behavior* New() const override { return new CellCellSpreadBehavior(); }
  Behavior* NewCopy() const override { return new CellCellSpreadBehavior(*this); }

  void Run(Agent* agent) override {
    // Skip if not on the scheduled interval — reduces O(N²) cost while
    // still resolving overlap frequently enough for numerical stability.
    if (run_interval_ > 1) {
      const auto* sched = bdm::Simulation::GetActive()->GetScheduler();
      if (sched && (sched->GetSimulatedSteps() % static_cast<uint64_t>(run_interval_) != 0)) return;
    }

    auto* cell = dynamic_cast<MyCell*>(agent);
    if (!cell) return;

    const Real3 my_pos = cell->GetPosition();
    const double my_r = cell->GetDiameter() * 0.5;
    // Neighbourhood cutoff — skip distant cells early to save N² time
    const double search_r = my_r * search_radius_factor_;

    Real3 displacement = {0.0, 0.0, 0.0};

    auto* rm = Simulation::GetActive()->GetResourceManager();
    rm->ForEachAgent([&](Agent* other_agent) {
      if (other_agent->GetUid() == agent->GetUid()) return;
      auto* other = dynamic_cast<MyCell*>(other_agent);
      if (!other) return;

      const Real3 other_pos = other->GetPosition();
      const double other_r = other->GetDiameter() * 0.5;

      // Quick bounding-box pre-rejection before the sqrt
      const double dx = my_pos[0] - other_pos[0];
      const double dy = my_pos[1] - other_pos[1];
      const double dz = my_pos[2] - other_pos[2];
      const double sum_r = my_r + other_r;
      if (std::abs(dx) > search_r || std::abs(dy) > search_r ||
          std::abs(dz) > search_r) {
        return;
      }

      const double dist_sq = dx * dx + dy * dy + dz * dz;
      const double dist = std::sqrt(dist_sq);
      const double overlap = sum_r - dist;
      const double adhesion_range =
          std::max(sum_r, sum_r * std::max(1.0, adhesion_range_factor_));

      if (dist < 1e-12) {
        // Coincident centres: push in a deterministic direction to break tie
        displacement[0] += spread_strength_ * std::max(sum_r, 1e-12) * 0.5;
        return;
      }

      if (overlap > 0.0) {
        // Each cell resolves half the overlap so the pair together closes it.
        const double push = spread_strength_ * overlap * 0.5;
        displacement[0] += (dx / dist) * push;
        displacement[1] += (dy / dist) * push;
        displacement[2] += (dz / dist) * push;
      } else if (adhesion_strength_ > 0.0 && dist < adhesion_range) {
        const double adhesion_band = std::max(adhesion_range - sum_r, 1e-12);
        const double closeness = (adhesion_range - dist) / adhesion_band;
        const double pull =
            adhesion_strength_ * std::clamp(closeness, 0.0, 1.0) *
            adhesion_band * 0.5;
        displacement[0] -= (dx / dist) * pull;
        displacement[1] -= (dy / dist) * pull;
        displacement[2] -= (dz / dist) * pull;
      }
    });

    const bool moved = (std::abs(displacement[0]) > 1e-15 ||
                        std::abs(displacement[1]) > 1e-15 ||
                        std::abs(displacement[2]) > 1e-15);
    if (!moved) return;

    const double max_disp =
        std::max(0.0, max_displacement_fraction_) * std::max(my_r, 1e-12);
    const double disp_norm = std::sqrt(displacement[0] * displacement[0] +
                                       displacement[1] * displacement[1] +
                                       displacement[2] * displacement[2]);
    if (max_disp > 0.0 && disp_norm > max_disp) {
      const double scale = max_disp / disp_norm;
      displacement[0] *= scale;
      displacement[1] *= scale;
      displacement[2] *= scale;
    }

    Real3 new_pos = my_pos;
    new_pos[0] += displacement[0];
    new_pos[1] += displacement[1];
    if (IsForce2D()) {
      // Keep 2D z-plane intact
      new_pos[2] = cell->GetOpenFoamCellCenter()[2];
    } else {
      new_pos[2] += displacement[2];
    }
    cell->SetPosition(new_pos);
  }
};

}  // namespace bdm

#endif  // CELL_SPREAD_BEHAVIOR_H_
