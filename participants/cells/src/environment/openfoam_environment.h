#ifndef FINITE_VOLUME_ENVIRONMENT_H_
#define FINITE_VOLUME_ENVIRONMENT_H_

#include "biodynamo.h"
#include <array>
#include <cmath>
#include <initializer_list>

namespace bdm {

// Lightweight load-balance stub that keeps rootcling from attempting to
// generate a dictionary for a private nested type.
class FiniteVolumeLoadBalanceInfo : public LoadBalanceInfo {
 public:
  void CallHandleIteratorConsumer(
      uint64_t start, uint64_t end,
      Functor<void, Iterator<AgentHandle>*>& consumer) const override {
    (void)start;
    (void)end;
    (void)consumer;
  }
};

// Environment stub that relies on OpenFOAM's finite-volume mesh instead of
// BioDynaMo's spatial data structures. It intentionally disables all neighbor
// queries and only keeps track of the overall simulation bounds so that other
// BioDynaMo subsystems can operate without accessing the traditional
// environments.
class FiniteVolumeEnvironment : public Environment {
 public:
  FiniteVolumeEnvironment(const Real3& domain_min, const Real3& domain_max) {
    AssignDomainBounds(domain_min, domain_max, /*mark_out_of_sync=*/false);
  }

  void SetDomainBounds(const Real3& domain_min, const Real3& domain_max) {
    AssignDomainBounds(domain_min, domain_max, /*mark_out_of_sync=*/true);
  }

  void ForEachNeighbor(Functor<void, Agent*, real_t>& lambda, const Agent& query,
                       real_t squared_radius) override {
    IgnoreNeighborQuery("Agent-Agent radius", lambda, query, squared_radius);
  }

  void ForEachNeighbor(Functor<void, Agent*>& lambda, const Agent& query,
                       void* criteria) override {
    IgnoreNeighborQuery("Agent-Agent criteria", lambda, query, criteria);
  }

  void ForEachNeighbor(Functor<void, Agent*, real_t>& lambda,
                       const Real3& query_position, real_t squared_radius,
                       const Agent* query_agent = nullptr) override {
    IgnoreNeighborQuery("Position-Agent radius", lambda, query_position,
                        squared_radius, query_agent);
  }

  void Clear() override {}

  std::array<int32_t, 6> GetDimensions() const override { return dimensions_; }

  std::array<int32_t, 2> GetDimensionThresholds() const override {
    return dimension_thresholds_;
  }

  LoadBalanceInfo* GetLoadBalanceInfo() override { return &load_balance_info_; }

  class DummyNeighborMutexBuilder : public NeighborMutexBuilder {
   public:
    class DummyNeighborMutex : public NeighborMutex {
     public:
      void lock() override {}
      void unlock() override {}
    };

    NeighborMutex* GetMutex(uint64_t box_idx) override {
      (void)box_idx;
      return &mutex_;
    }

   private:
    DummyNeighborMutex mutex_;
  };

  NeighborMutexBuilder* GetNeighborMutexBuilder() override {
    return &neighbor_mutex_builder_;
  }

 protected:
  void UpdateImplementation() override {
    // Keep bounding box consistent with current agent positions while still
    // respecting the initial mesh bounds defined by OpenFOAM.
    std::array<real_t, 6> bounds = {domain_min_[0], domain_max_[0],
                                    domain_min_[1], domain_max_[1],
                                    domain_min_[2], domain_max_[2]};
    largest_object_size_ = 0.0;
    largest_object_size_squared_ = 0.0;
    this->CalcSimDimensionsAndLargestAgent(&bounds);

    domain_min_[0] = bounds[0];
    domain_max_[0] = bounds[1];
    domain_min_[1] = bounds[2];
    domain_max_[1] = bounds[3];
    domain_min_[2] = bounds[4];
    domain_max_[2] = bounds[5];

    UpdateDimensionCaches(domain_min_, domain_max_);
    has_grown_ = false;
  }

 private:
  void UpdateDimensionCaches(const Real3& domain_min, const Real3& domain_max) {
    for (int i = 0; i < 3; ++i) {
      dimensions_[2 * i] =
          static_cast<int32_t>(std::floor(domain_min[i]));
      dimensions_[2 * i + 1] =
          static_cast<int32_t>(std::ceil(domain_max[i]));
    }
    dimension_thresholds_[0] = dimensions_[0];
    dimension_thresholds_[1] = dimensions_[1];
  }

  template <typename... Args>
  void IgnoreNeighborQuery(const char* tag, Args&&... ignored_args) const {
    (void)std::initializer_list<int>{((void)ignored_args, 0)...};
    Log::Debug("FiniteVolumeEnvironment",
               "Neighbor query (", tag,
               ") ignored because OpenFOAM handles spatial relationships.");
  }

  void AssignDomainBounds(const Real3& domain_min, const Real3& domain_max,
                          bool mark_out_of_sync) {
    domain_min_ = domain_min;
    domain_max_ = domain_max;
    UpdateDimensionCaches(domain_min_, domain_max_);
    if (mark_out_of_sync) {
      this->MarkAsOutOfSync();
    }
  }

  Real3 domain_min_;
  Real3 domain_max_;
  std::array<int32_t, 6> dimensions_{};
  std::array<int32_t, 2> dimension_thresholds_{};
  FiniteVolumeLoadBalanceInfo load_balance_info_;
  DummyNeighborMutexBuilder neighbor_mutex_builder_;
};

}  // namespace bdm

#endif  // FINITE_VOLUME_ENVIRONMENT_H_
