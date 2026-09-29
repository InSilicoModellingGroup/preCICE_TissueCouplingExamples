#ifndef BDM_ADAPTER_MESH_MAPPING_H_
#define BDM_ADAPTER_MESH_MAPPING_H_

#include <vector>
#include <string>
#include <map>
#include <array>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <limits>

namespace bdm_adapter {

/**
 * @brief Bounding box for mesh/domain bounds
 *
 * NO DEFAULT VALUES - must be explicitly set from mesh data.
 * Initialized to NaN to catch uninitialized usage.
 */
struct BoundingBox {
    std::array<double, 3> min = {std::numeric_limits<double>::quiet_NaN(),
                                  std::numeric_limits<double>::quiet_NaN(),
                                  std::numeric_limits<double>::quiet_NaN()};
    std::array<double, 3> max = {std::numeric_limits<double>::quiet_NaN(),
                                  std::numeric_limits<double>::quiet_NaN(),
                                  std::numeric_limits<double>::quiet_NaN()};

    bool Contains(const std::array<double, 3>& point) const {
        return point[0] >= min[0] && point[0] <= max[0] &&
               point[1] >= min[1] && point[1] <= max[1] &&
               point[2] >= min[2] && point[2] <= max[2];
    }

    std::array<double, 3> Center() const {
        return {
            (min[0] + max[0]) / 2.0,
            (min[1] + max[1]) / 2.0,
            (min[2] + max[2]) / 2.0
        };
    }

    std::array<double, 3> Size() const {
        return {
            max[0] - min[0],
            max[1] - min[1],
            max[2] - min[2]
        };
    }

    double Volume() const {
        auto s = Size();
        return s[0] * s[1] * s[2];
    }

    void Expand(double margin) {
        min[0] -= margin; min[1] -= margin; min[2] -= margin;
        max[0] += margin; max[1] += margin; max[2] += margin;
    }
};

/**
 * @brief Mesh data container for preCICE meshes
 */
struct MeshData {
    std::string name;
    std::vector<int> vertex_ids;
    std::vector<double> coordinates;  // Interleaved: x0,y0,z0,x1,y1,z1,...
    std::vector<double> mapped_cell_centers;  // Optional owner cell centers for boundary vertices
    std::vector<double> face_areas;  // Optional face area for each boundary vertex
    BoundingBox bounds;
    bool is_initialized = false;

    size_t NumVertices() const {
        return vertex_ids.size();
    }

    std::array<double, 3> GetVertex(size_t index) const {
        size_t i = index * 3;
        return {coordinates[i], coordinates[i + 1], coordinates[i + 2]};
    }

    bool HasMappedCellCenters() const {
        return mapped_cell_centers.size() == coordinates.size() && !mapped_cell_centers.empty();
    }

    std::array<double, 3> GetMappedCellCenter(size_t index) const {
        size_t i = index * 3;
        return {mapped_cell_centers[i], mapped_cell_centers[i + 1], mapped_cell_centers[i + 2]};
    }

    bool HasFaceAreas() const {
        return face_areas.size() == vertex_ids.size() && !face_areas.empty();
    }

    void ComputeBounds() {
        if (coordinates.empty()) return;

        bounds.min = {coordinates[0], coordinates[1], coordinates[2]};
        bounds.max = bounds.min;

        for (size_t i = 0; i < coordinates.size(); i += 3) {
            bounds.min[0] = std::min(bounds.min[0], coordinates[i]);
            bounds.min[1] = std::min(bounds.min[1], coordinates[i + 1]);
            bounds.min[2] = std::min(bounds.min[2], coordinates[i + 2]);
            bounds.max[0] = std::max(bounds.max[0], coordinates[i]);
            bounds.max[1] = std::max(bounds.max[1], coordinates[i + 1]);
            bounds.max[2] = std::max(bounds.max[2], coordinates[i + 2]);
        }
    }
};

/**
 * @brief Wall mesh data with geometry information
 */
struct WallMeshData : public MeshData {
    std::vector<double> normals;    // Face normals (3 per face)
    std::vector<double> areas;      // Face areas (1 per face)
    std::vector<double> geometry;   // Face centers (3 per face)
    std::vector<double> velocity;   // Face velocities (3 per face) for moving walls

    size_t NumFaces() const {
        return areas.size();
    }

    std::array<double, 3> GetNormal(size_t index) const {
        size_t i = index * 3;
        return {normals[i], normals[i + 1], normals[i + 2]};
    }

    std::array<double, 3> GetFaceCenter(size_t index) const {
        size_t i = index * 3;
        return {geometry[i], geometry[i + 1], geometry[i + 2]};
    }

    std::array<double, 3> GetVelocity(size_t index) const {
        if (velocity.empty()) {
            return {0.0, 0.0, 0.0};  // Stationary wall default
        }
        size_t i = index * 3;
        return {velocity[i], velocity[i + 1], velocity[i + 2]};
    }

    bool HasVelocity() const {
        return !velocity.empty();
    }
};

/**
 * @brief Exact fixed boundary triangles received as per-record preCICE data.
 */
struct BoundaryTriangleMeshData : public MeshData {
    std::vector<double> vertex0;
    std::vector<double> vertex1;
    std::vector<double> vertex2;
    std::vector<double> normals;
    std::vector<double> velocity;
    std::vector<double> kinds;

    bool HasCompleteRecords() const {
        const std::size_t n = NumVertices();
        return n > 0
            && vertex0.size() == n * 3
            && vertex1.size() == n * 3
            && vertex2.size() == n * 3
            && normals.size() == n * 3
            && kinds.size() == n;
    }
};

/**
 * @brief Parse OpenFOAM blockMeshDict for domain bounds
 *
 * @param filepath Path to blockMeshDict file
 * @return BoundingBox Domain bounds
 */
inline BoundingBox ParseBlockMeshDict(const std::string& filepath) {
    BoundingBox box;
    std::ifstream file(filepath);

    if (!file.is_open()) {
        return box;
    }

    std::string content((std::istreambuf_iterator<char>(file)),
                        std::istreambuf_iterator<char>());

    // Find vertices section
    auto vertices_pos = content.find("vertices");
    if (vertices_pos == std::string::npos) {
        return box;
    }

    auto start = content.find('(', vertices_pos);
    auto end = content.find(");", start);

    if (start == std::string::npos || end == std::string::npos) {
        return box;
    }

    std::string vertices_section = content.substr(start + 1, end - start - 1);

    // Parse vertex coordinates
    std::vector<std::array<double, 3>> vertices;
    std::istringstream iss(vertices_section);
    std::string line;

    while (std::getline(iss, line)) {
        // Remove comments
        auto comment_pos = line.find("//");
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        // Find coordinate tuple
        auto open_paren = line.find('(');
        auto close_paren = line.find(')');

        if (open_paren != std::string::npos && close_paren != std::string::npos) {
            std::string coords = line.substr(open_paren + 1, close_paren - open_paren - 1);
            std::istringstream coord_stream(coords);
            double x, y, z;
            if (coord_stream >> x >> y >> z) {
                vertices.push_back({x, y, z});
            }
        }
    }

    if (vertices.empty()) {
        return box;
    }

    // Compute bounding box
    box.min = vertices[0];
    box.max = vertices[0];

    for (const auto& v : vertices) {
        box.min[0] = std::min(box.min[0], v[0]);
        box.min[1] = std::min(box.min[1], v[1]);
        box.min[2] = std::min(box.min[2], v[2]);
        box.max[0] = std::max(box.max[0], v[0]);
        box.max[1] = std::max(box.max[1], v[1]);
        box.max[2] = std::max(box.max[2], v[2]);
    }

    return box;
}

/**
 * @brief Load agent positions from CSV file
 *
 * @param filepath Path to CSV file
 * @return Vector of positions as {x, y, z} arrays
 */
inline std::vector<std::array<double, 3>> LoadPositionsFromCSV(const std::string& filepath) {
    std::vector<std::array<double, 3>> positions;
    std::ifstream file(filepath);

    if (!file.is_open()) {
        return positions;
    }

    std::string line;
    bool header_skipped = false;

    while (std::getline(file, line)) {
        // Skip header line
        if (!header_skipped) {
            if (line.find("x") != std::string::npos || line.find("X") != std::string::npos) {
                header_skipped = true;
                continue;
            }
            header_skipped = true;
        }

        // Parse CSV line
        std::istringstream iss(line);
        std::string token;
        std::array<double, 3> pos = {0.0, 0.0, 0.0};
        int col = 0;

        while (std::getline(iss, token, ',') && col < 3) {
            try {
                pos[col++] = std::stod(token);
            } catch (...) {
                break;
            }
        }

        if (col >= 2) {  // At least x, y
            positions.push_back(pos);
        }
    }

    return positions;
}

/**
 * @brief Generate grid positions within bounding box
 *
 * @param box Bounding box
 * @param nx Number of points in x direction
 * @param ny Number of points in y direction
 * @param nz Number of points in z direction
 * @param margin Margin from boundaries
 * @return Vector of positions
 */
inline std::vector<std::array<double, 3>> GenerateGridPositions(
    const BoundingBox& box,
    int nx, int ny, int nz,
    double margin = 0.0)
{
    std::vector<std::array<double, 3>> positions;

    double x_min = box.min[0] + margin;
    double x_max = box.max[0] - margin;
    double y_min = box.min[1] + margin;
    double y_max = box.max[1] - margin;
    double z_min = box.min[2] + margin;
    double z_max = box.max[2] - margin;

    double dx = (nx > 1) ? (x_max - x_min) / (nx - 1) : 0.0;
    double dy = (ny > 1) ? (y_max - y_min) / (ny - 1) : 0.0;
    double dz = (nz > 1) ? (z_max - z_min) / (nz - 1) : 0.0;

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                positions.push_back({
                    x_min + i * dx,
                    y_min + j * dy,
                    z_min + k * dz
                });
            }
        }
    }

    return positions;
}

/**
 * @brief Filter out positions that are too close to wall mesh face centres
 *
 * For each candidate position the minimum distance to every face centre in
 * every supplied WallMeshData is computed.  Positions whose minimum distance
 * is less than @p clearance are discarded.
 *
 * The function uses the @c WallMeshData::geometry field (face centres, 3
 * doubles per face) when available.  If @c geometry is empty it falls back
 * to the generic vertex @c coordinates instead.
 *
 * @param positions  Input candidate positions (x,y,z).
 * @param wall_meshes  Map of wall mesh name -> WallMeshData.
 * @param clearance  Minimum acceptable distance from any wall face [m].
 * @return Filtered position vector (positions farther than clearance from all walls).
 */
inline std::vector<std::array<double, 3>> FilterWallAdjacentPositions(
    const std::vector<std::array<double, 3>>& positions,
    const std::map<std::string, WallMeshData>&  wall_meshes,
    double clearance)
{
    if (clearance <= 0.0 || wall_meshes.empty()) {
        return positions;  // nothing to filter
    }

    // Collect wall sample points from BOTH face centres and vertices.
    // Using both is more robust near edges/corners where centre-only sampling
    // may miss close positions.
    std::vector<std::array<double, 3>> wall_points;
    for (const auto& kv : wall_meshes) {
        const WallMeshData& wm = kv.second;
        if (!wm.is_initialized) continue;

        if (!wm.geometry.empty()) {
            // geometry = face centres laid out as [x0,y0,z0, x1,y1,z1, ...]
            const size_t n_faces = wm.geometry.size() / 3;
            wall_points.reserve(wall_points.size() + n_faces);
            for (size_t i = 0; i < n_faces; ++i) {
                wall_points.push_back({
                    wm.geometry[i * 3],
                    wm.geometry[i * 3 + 1],
                    wm.geometry[i * 3 + 2]
                });
            }
        }

        if (!wm.coordinates.empty()) {
            // Also include wall vertices for tighter clearance checks.
            const size_t n_verts = wm.coordinates.size() / 3;
            wall_points.reserve(wall_points.size() + n_verts);
            for (size_t i = 0; i < n_verts; ++i) {
                wall_points.push_back({
                    wm.coordinates[i * 3],
                    wm.coordinates[i * 3 + 1],
                    wm.coordinates[i * 3 + 2]
                });
            }
        }
    }

    if (wall_points.empty()) {
        return positions;  // no wall data yet – skip filtering
    }

    const double clearance_sq = clearance * clearance;
    std::vector<std::array<double, 3>> filtered;
    filtered.reserve(positions.size());

    for (const auto& pos : positions) {
        double min_dist_sq = std::numeric_limits<double>::max();
        for (const auto& wp : wall_points) {
            double dx = pos[0] - wp[0];
            double dy = pos[1] - wp[1];
            double dz = pos[2] - wp[2];
            double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 < min_dist_sq) min_dist_sq = d2;
        }
        if (min_dist_sq >= clearance_sq) {
            filtered.push_back(pos);
        }
    }

    return filtered;
}

} // namespace bdm_adapter

#endif // BDM_ADAPTER_MESH_MAPPING_H_
