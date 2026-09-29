#ifndef AGENT_MAPPING_H_
#define AGENT_MAPPING_H_

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <functional>
#include <unordered_map>
#include <utility>
#include <vector>
#include <string>

#include "biodynamo.h"
#ifndef __CLING__
#include "precice/Participant.hpp"
#include "../agents/coupled_cell.h"

namespace bdm {

/**
 * @brief Agent mapping using BioDynaMo's spatial structures
 *
 * Leverages BioDynaMo's built-in spatial queries and ResourceManager
 * for efficient agent-mesh mapping operations.
 */
class AgentMapping {
 public:
  struct MeshData {
    MeshData() = default;
    explicit MeshData(std::string mesh_name, bool volume_mesh = true)
        : name(std::move(mesh_name)), is_volume_mesh(volume_mesh) {}

    std::vector<precice::VertexID> vertex_ids;
    std::vector<double> coordinates;  // flat buffer, stride = coordinate_stride
    size_t vertex_count = 0;
    int dimensions = 3;
    int coordinate_stride = 3;
    bool use_z_override = false;
    double z_override = 0.0;
    std::string name;
    bool is_initialized = false;
    bool is_volume_mesh = false;

    double GetCoordinate(size_t idx, int component) const {
      if (component < 0 || component >= 3) {
        if (component == 2 && use_z_override) {
          return z_override;
        }
        return 0.0;
      }
      if (idx >= vertex_count) {
        return 0.0;
      }
      const size_t stride = static_cast<size_t>(std::max(1, coordinate_stride));
      const size_t offset = idx * stride + static_cast<size_t>(component);
      return offset < coordinates.size() ? coordinates[offset] : 0.0;
    }

    Real3 GetPosition(size_t idx) const {
      Real3 pos{0.0, 0.0, 0.0};
      pos[0] = GetCoordinate(idx, 0);
      pos[1] = GetCoordinate(idx, 1);
      pos[2] = GetCoordinate(idx, 2);
      if (use_z_override) {
        pos[2] = z_override;
      }
      return pos;
    }

    Real3 AlignToPlane(const Real3& input) const {
      Real3 adjusted = input;
      if (use_z_override) {
        adjusted[2] = z_override;
      }
      return adjusted;
    }

    void ConfigurePlanarOverride(bool enable, double plane) {
      use_z_override = enable;
      if (enable) {
        z_override = plane;
      }
    }
  };

  // Use BioDynaMo's Environment for spatial queries
  static MeshData ReadMeshFromPrecice(precice::Participant& interface,
                                      const std::string& mesh_name) {
  MeshData mesh(mesh_name);
  mesh.is_volume_mesh = true;

    Log::Info("AgentMapping", "Reading mesh: ", mesh_name);

    // Determine mesh dimensionality and allocate storage
    mesh.vertex_count = interface.getMeshVertexSize(mesh_name);
    if (mesh.vertex_count == 0) {
      Log::Warning("AgentMapping", "No vertices in mesh ", mesh_name);
      return mesh;
    }

    int dims = 3;
    try {
      dims = interface.getMeshDimensions(mesh_name);
    } catch (const std::exception& e) {
      Log::Warning("AgentMapping",
                   "Falling back to 3D coordinates for mesh '", mesh_name,
                   "' after failing to query dimensions: ", e.what());
      dims = 3;
    }
    mesh.dimensions = std::max(1, dims);

    // Always allocate with stride 3 to avoid buffer overrun when preCICE writes
    // 3 coordinates per vertex even for 2D meshes (as noted in the comment below)
    const int coordinate_stride_for_allocation = 3;
    mesh.vertex_ids.resize(mesh.vertex_count);
    std::vector<double> raw_coordinates(
        mesh.vertex_count * static_cast<size_t>(coordinate_stride_for_allocation), 0.0);

    // Retrieve vertex IDs and coordinates using modern preCICE API. Some
    // preCICE builds always write three coordinates per vertex, even for meshes
    // that report a lower dimensionality (e.g. 2D). We therefore read into a
    // buffer with a stride of three and then keep that stride internally to
    // avoid buffer overruns and the resulting memory corruption.
    interface.getMeshVertexIDsAndCoordinates(
        mesh_name,
        precice::span<precice::VertexID>(mesh.vertex_ids.data(),
                                         mesh.vertex_ids.size()),
        precice::span<double>(raw_coordinates.data(),
                              raw_coordinates.size()));

    mesh.coordinates = std::move(raw_coordinates);
    mesh.is_initialized = mesh.vertex_count > 0;
    if (mesh.vertex_count > 0) {
      const size_t stride = mesh.coordinates.size() / mesh.vertex_count;
      mesh.coordinate_stride = stride > 0 ? static_cast<int>(stride)
                                          : std::max(1, mesh.dimensions);
    }

    Log::Info("AgentMapping", "Read ", mesh.vertex_count,
              " vertices for mesh ", mesh_name, " (dims=", mesh.dimensions,
              ", stride=", mesh.coordinate_stride, ")");
    if (mesh.vertex_count > 0) {
      const Real3 first = mesh.GetPosition(0);
      const Real3 last = mesh.GetPosition(mesh.vertex_count - 1);
      Log::Info("AgentMapping", "  First vertex at (", first[0], ", ", first[1],
                ", ", first[2], ")");
      Log::Info("AgentMapping", "  Last vertex at (", last[0], ", ", last[1],
                ", ", last[2], ")");
    }
    return mesh;
  }

  /// Create one BioDynaMo agent at the center of each volume cell
  /// represented by the provided mesh. The mesh vertices are treated as
  /// cell centroids. Each created agent is initialized at the vertex
  /// coordinates and mapped to the corresponding preCICE vertex ID.
  ///
  /// @param mesh       Volume mesh data containing vertex IDs and coordinates
  /// @param diameter   Initial diameter for the created agents
  /// @return           Mapping from created AgentUid to preCICE vertex ID
  static std::vector<std::pair<AgentUid, precice::VertexID>>
  CreateAgentsFromMesh(const MeshData& mesh, double diameter,
                       size_t sample_stride = 1,
                       const std::function<void(MyCell&, int)>& configure = {},
                       bool force2d = false, double z_plane = 0.0,
                       const Real3* domain_min = nullptr,
                       const Real3* domain_max = nullptr) {
    // Debug log for mesh data
    Log::Info("AgentMapping", "CreateAgentsFromMesh called with mesh vertex_count=",
              mesh.vertex_count, ", diameter=", diameter, ", sample_stride=", sample_stride);
    std::vector<std::pair<AgentUid, precice::VertexID>> mapping;
    if (mesh.vertex_count == 0) {
      Log::Fatal("AgentMapping", "No mesh vertices available from OpenFOAM. Cannot create agents without mesh data.");
      return mapping; // Empty mapping
    }

    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();
    if (sample_stride == 0) {
      sample_stride = 1;
    }
    mapping.reserve((mesh.vertex_count + sample_stride - 1) / sample_stride);

    std::vector<Real3> debug_samples;
    debug_samples.reserve(std::min<size_t>(5, mesh.vertex_count));

    // Create agents at each OpenFOAM volume cell center (represented by mesh vertices)
    Log::Info("AgentMapping", "Creating agents at the center of each OpenFOAM volume cell");

    for (size_t i = 0; i < mesh.vertex_count; i += sample_stride) {
      // Get position from mesh (this is the center of the OpenFOAM volume cell)
      Real3 pos = mesh.GetPosition(i);

      // Handle 2D/3D alignment
      if (force2d) {
        pos[2] = z_plane;
      } else {
        pos = mesh.AlignToPlane(pos);
      }

      // Ensure position is within domain bounds
      if (domain_min && domain_max) {
        for (int axis = 0; axis < 3; ++axis) {
          pos[axis] = std::clamp(pos[axis], (*domain_min)[axis], (*domain_max)[axis]);
        }
      }

      // Create agent at the volume cell center
      auto* cell = new MyCell(pos);

      // Set diameter based on OpenFOAM cell size
      cell->SetDiameter(diameter);

      // Get the vertex ID from the mesh
      precice::VertexID vid = mesh.vertex_ids[i];
      cell->SetMeshVertexId(static_cast<int>(vid));

      // Apply any additional configuration (like motion behavior)
      if (configure) {
        configure(*cell, static_cast<int>(vid));
      }

      // Add to simulation and mapping
      rm->AddAgent(cell);
      mapping.emplace_back(cell->GetUid(), vid);

      // Store sample positions for debugging
      if (debug_samples.size() < 5) {
        debug_samples.push_back(cell->GetPosition());
      }
    }

    Log::Info("AgentMapping", "Created ", mapping.size(),
              " agents from volume mesh centres (stride=", sample_stride, ")");
    if (!debug_samples.empty()) {
      for (size_t i = 0; i < debug_samples.size(); ++i) {
        const auto& sample = debug_samples[i];
        Log::Info("AgentMapping", "  Sample agent[", i,
                  "] position = (", sample[0], ", ", sample[1], ", ",
                  sample[2], ")");
      }
    }
    return mapping;
  }


  // Nearest-cell assignment:
  //   c(i) = argmin_c ||x_i - x_c||
  // where mesh vertices are the exported OpenFOAM cell centers x_c.
  static int
  FindNearestMeshVertex(const Real3& position, const MeshData& mesh) {
    if (mesh.vertex_count == 0)
      return -1;

    double min_dist_sq = std::numeric_limits<double>::max();
    int nearest_idx = -1;

    Real3 query = mesh.AlignToPlane(position);

    for (size_t i = 0; i < mesh.vertex_count; ++i) {
      Real3 vertex = mesh.GetPosition(i);

      double dx = query[0] - vertex[0];
      double dy = query[1] - vertex[1];
      double dz = query[2] - vertex[2];
      double dist_sq = dx * dx + dy * dy + dz * dz;
      if (dist_sq < min_dist_sq) {
        min_dist_sq = dist_sq;
        nearest_idx = static_cast<int>(i);
      }
    }

    return nearest_idx;
  }

  // Agent-to-mesh mapping using simple domain bounds
  static std::vector<std::pair<AgentUid, precice::VertexID>>
  MapAgentsToMesh(
      const MeshData& mesh, const Real3& domain_min, const Real3& domain_max,
      std::unordered_map<uint64_t, precice::VertexID>* mapping_out = nullptr,
      const std::unordered_map<precice::VertexID, size_t>* vertex_lookup =
          nullptr) {
    auto* sim = Simulation::GetActive();
    auto* rm = sim->GetResourceManager();

    std::vector<std::pair<AgentUid, precice::VertexID>> mapping;
    mapping.reserve(rm->GetNumAgents());

    std::unordered_map<precice::VertexID, size_t> local_lookup;
    if (!vertex_lookup) {
      local_lookup.reserve(mesh.vertex_ids.size());
      for (size_t i = 0; i < mesh.vertex_ids.size(); ++i) {
        local_lookup[mesh.vertex_ids[i]] = i;
      }
      vertex_lookup = &local_lookup;
    }

    rm->ForEachAgent([&](Agent* agent) {
      auto* cell = dynamic_cast<MyCell*>(agent);
      if (!cell)
        return;

      // Skip helper agents (e.g., bounding box corners) that do not participate in coupling
      if (cell->GetPreciceVertexId() < 0 || cell->IsStationary()) {
        return;
      }

  precice::VertexID vid = static_cast<precice::VertexID>(cell->GetMeshVertexId());
      Real3 pos = mesh.AlignToPlane(cell->GetPosition());
      for (int axis = 0; axis < 3; ++axis) {
        pos[axis] = std::clamp(pos[axis], domain_min[axis], domain_max[axis]);
      }

      auto ensure_vertex_valid = [&](precice::VertexID vertex_id) {
        if (vertex_id < 0 || !vertex_lookup) {
          return false;
        }
        return vertex_lookup->find(vertex_id) != vertex_lookup->end();
      };

      auto store_mapping = [&](precice::VertexID vertex_id) {
        if (vertex_id < 0) {
          return;
        }
  cell->SetMeshVertexId(static_cast<int>(vertex_id));
        mapping.emplace_back(cell->GetUid(), vertex_id);
        if (mapping_out) {
          (*mapping_out)[cell->GetUid().GetIndex()] = vertex_id;
        }
      };

      auto ensure_vertex_index = [&](precice::VertexID vertex_id) -> const size_t* {
        if (!ensure_vertex_valid(vertex_id)) {
          return nullptr;
        }
        auto it = vertex_lookup->find(vertex_id);
        if (it == vertex_lookup->end()) {
          return nullptr;
        }
        return &it->second;
      };

      auto is_within_neighborhood = [&](precice::VertexID vertex_id) {
        const size_t* idx_ptr = ensure_vertex_index(vertex_id);
        if (!idx_ptr) {
          return false;
        }
        size_t idx = *idx_ptr;
        Real3 vertex = mesh.GetPosition(idx);
        double dx = pos[0] - vertex[0];
        double dy = pos[1] - vertex[1];
        double dz = pos[2] - vertex[2];
        double dist_sq = dx * dx + dy * dy + dz * dz;

        double radius = cell->GetDiameter() * 0.5;
        double tolerance = std::max(radius, 1e-9);
        return dist_sq <= tolerance * tolerance;
      };

      auto attempt_store_with_validation = [&](precice::VertexID vertex_id) {
        if (!ensure_vertex_valid(vertex_id)) {
          return false;
        }

        precice::VertexID resolved_vid = vertex_id;
        if (!is_within_neighborhood(resolved_vid)) {
          int nearest_vertex = FindNearestMeshVertex(pos, mesh);
          if (nearest_vertex >= 0) {
            resolved_vid = mesh.vertex_ids[nearest_vertex];
          }
        }

        if (!ensure_vertex_valid(resolved_vid)) {
          return false;
        }

        store_mapping(resolved_vid);
        return true;
      };

      if (attempt_store_with_validation(vid)) {
        return;
      }

      if (mapping_out) {
        auto it = mapping_out->find(cell->GetUid().GetIndex());
        if (it != mapping_out->end() && attempt_store_with_validation(it->second)) {
          return;
        }
      }

      int nearest_vertex = FindNearestMeshVertex(pos, mesh);
      if (nearest_vertex >= 0) {
        precice::VertexID nearest_vid = mesh.vertex_ids[nearest_vertex];
        store_mapping(nearest_vid);
      } else if (mapping_out) {
        mapping_out->erase(cell->GetUid().GetIndex());
      }
    });

    Log::Info("AgentMapping", "Mapped ", mapping.size(),
              " agents to mesh vertices");
    return mapping;
  }
};

} // namespace bdm
#endif // __CLING__

#endif // AGENT_MAPPING_H_
