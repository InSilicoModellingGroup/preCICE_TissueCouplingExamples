#ifndef BDM_ADAPTER_BOUNDARY_GEOMETRY_H_
#define BDM_ADAPTER_BOUNDARY_GEOMETRY_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
#include <vector>

namespace bdm_adapter {

using BoundaryVec3 = std::array<double, 3>;

inline BoundaryVec3 BAdd(const BoundaryVec3& a, const BoundaryVec3& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}
inline BoundaryVec3 BSub(const BoundaryVec3& a, const BoundaryVec3& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
inline BoundaryVec3 BScale(const BoundaryVec3& a, double s) {
    return {a[0] * s, a[1] * s, a[2] * s};
}
inline double BDot(const BoundaryVec3& a, const BoundaryVec3& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
inline BoundaryVec3 BCross(const BoundaryVec3& a, const BoundaryVec3& b) {
    return {a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]};
}
inline double BNormSq(const BoundaryVec3& a) { return BDot(a, a); }
inline double BNorm(const BoundaryVec3& a) { return std::sqrt(BNormSq(a)); }
inline BoundaryVec3 BNormalized(const BoundaryVec3& a) {
    const double n = BNorm(a);
    return n > 1e-30 ? BScale(a, 1.0 / n) : BoundaryVec3{0, 0, 0};
}

enum class BoundaryKind : unsigned char {
    kPassThrough = 0,
    kWall = 1
};

struct BoundaryTriangle {
    BoundaryVec3 a{};
    BoundaryVec3 b{};
    BoundaryVec3 c{};
    BoundaryVec3 normal{};
    BoundaryVec3 velocity{};
    BoundaryKind kind = BoundaryKind::kWall;
};

struct BoundaryClosestPoint {
    bool found = false;
    std::size_t triangle = 0;
    BoundaryVec3 point{};
    BoundaryVec3 normal{};
    BoundaryVec3 velocity{};
    double distance = std::numeric_limits<double>::infinity();
};

struct BoundaryContact : BoundaryClosestPoint {
    double fraction = 1.0;
};

enum class BoundaryVelocityMode {
    kNone,
    kSlip,
    kRebound,
    kStick
};

struct BoundaryResponse {
    BoundaryVec3 relativeVelocity{};
    BoundaryVec3 remainingDisplacement{};
};

inline BoundaryResponse ApplyBoundaryResponse(
    BoundaryVelocityMode mode,
    BoundaryVec3 relativeVelocity,
    BoundaryVec3 remainingDisplacement,
    const BoundaryVec3& outwardNormal,
    const BoundaryVec3& wallVelocity,
    double remainingTime,
    double restitution,
    double friction) {
    const double normalVelocity = BDot(relativeVelocity, outwardNormal);
    const double normalRemaining =
        BDot(remainingDisplacement, outwardNormal);

    switch (mode) {
        case BoundaryVelocityMode::kSlip:
            if (normalVelocity > 0.0) {
                relativeVelocity = BSub(
                    relativeVelocity,
                    BScale(outwardNormal, normalVelocity));
            }
            if (normalRemaining > 0.0) {
                remainingDisplacement = BSub(
                    remainingDisplacement,
                    BScale(outwardNormal, normalRemaining));
            }
            break;
        case BoundaryVelocityMode::kRebound:
            if (normalVelocity > 0.0) {
                const BoundaryVec3 tangent = BSub(
                    relativeVelocity,
                    BScale(outwardNormal, normalVelocity));
                relativeVelocity = BSub(
                    BSub(relativeVelocity,
                         BScale(outwardNormal,
                                (1.0 + restitution) * normalVelocity)),
                    BScale(tangent, friction));
            }
            if (normalRemaining > 0.0) {
                remainingDisplacement = BSub(
                    remainingDisplacement,
                    BScale(outwardNormal,
                           (1.0 + restitution) * normalRemaining));
            }
            break;
        case BoundaryVelocityMode::kStick:
            relativeVelocity = {0.0, 0.0, 0.0};
            remainingDisplacement = BScale(wallVelocity, remainingTime);
            break;
        case BoundaryVelocityMode::kNone:
            remainingDisplacement = {0.0, 0.0, 0.0};
            break;
    }
    return {relativeVelocity, remainingDisplacement};
}

class BoundaryGeometry {
 public:
    bool Build(std::vector<BoundaryTriangle> triangles) {
        triangles_.clear();
        triangles_.reserve(triangles.size());
        for (auto& tri : triangles) {
            BoundaryVec3 geometricNormal =
                BCross(BSub(tri.b, tri.a), BSub(tri.c, tri.a));
            if (BNormSq(geometricNormal) <= 1e-28) {
                continue;
            }
            geometricNormal = BNormalized(geometricNormal);
            if (BNormSq(tri.normal) <= 1e-28) {
                tri.normal = geometricNormal;
            } else {
                tri.normal = BNormalized(tri.normal);
                if (BDot(tri.normal, geometricNormal) < 0.0) {
                    std::swap(tri.b, tri.c);
                }
            }
            triangles_.push_back(tri);
        }

        indices_.resize(triangles_.size());
        std::iota(indices_.begin(), indices_.end(), std::size_t{0});
        nodes_.clear();
        if (!indices_.empty()) {
            BuildNode(0, indices_.size());
        }
        using Edge = std::pair<BoundaryVec3, BoundaryVec3>;
        std::map<Edge, int> edgeCounts;
        const auto addEdge = [&](BoundaryVec3 a, BoundaryVec3 b) {
            if (b < a) std::swap(a, b);
            ++edgeCounts[{a, b}];
        };
        for (const auto& tri : triangles_) {
            addEdge(tri.a, tri.b);
            addEdge(tri.b, tri.c);
            addEdge(tri.c, tri.a);
        }
        closed_ = !edgeCounts.empty();
        for (const auto& entry : edgeCounts) {
            if (entry.second != 2) {
                closed_ = false;
                break;
            }
        }
        return !triangles_.empty();
    }

    bool Empty() const { return triangles_.empty(); }
    std::size_t Size() const { return triangles_.size(); }
    bool IsClosed() const { return closed_; }
    const std::vector<BoundaryTriangle>& Triangles() const { return triangles_; }

    BoundaryClosestPoint ClosestPoint(const BoundaryVec3& p,
                                      bool wallsOnly = false) const {
        BoundaryClosestPoint result;
        if (nodes_.empty()) return result;
        ClosestNode(0, p, wallsOnly, result);
        return result;
    }

    // A closed boundary is required for a globally meaningful result.
    bool IsInside(const BoundaryVec3& p) const {
        const BoundaryVec3 direction =
            BNormalized(BoundaryVec3{1.0, 0.371390676, 0.173205081});
        const std::size_t intersections =
            nodes_.empty() ? 0 : CountRayIntersections(0, p, direction);
        return (intersections % 2) == 1;
    }

    // Positive in the closed computational domain, negative outside.
    double SignedDistance(const BoundaryVec3& p) const {
        const auto closest = ClosestPoint(p, false);
        if (!closest.found) return -std::numeric_limits<double>::infinity();
        return IsInside(p) ? closest.distance : -closest.distance;
    }

    BoundaryContact FindFirstWallContact(const BoundaryVec3& start,
                                         const BoundaryVec3& finish,
                                         double radius,
                                         double tolerance = 1e-10) const {
        BoundaryContact contact;
        if (nodes_.empty()) return contact;
        const BoundaryVec3 displacement = BSub(finish, start);
        const double pathLength = BNorm(displacement);
        if (pathLength <= tolerance) {
            const auto cp = ClosestPoint(start, true);
            if (cp.found && cp.distance <= radius + tolerance) {
                static_cast<BoundaryClosestPoint&>(contact) = cp;
                contact.fraction = 0.0;
            }
            return contact;
        }

        Box sweptBox;
        for (int d = 0; d < 3; ++d) {
            sweptBox.min[d] = std::min(start[d], finish[d]) - radius - tolerance;
            sweptBox.max[d] = std::max(start[d], finish[d]) + radius + tolerance;
        }
        std::vector<std::size_t> candidates;
        CollectOverlappingWalls(0, sweptBox, &candidates);
        const double thresholdSq =
            (radius + tolerance) * (radius + tolerance);
        double bestFraction = std::numeric_limits<double>::infinity();

        for (const std::size_t index : candidates) {
            const auto& tri = triangles_[index];
            const auto distanceSqAt = [&](double t) {
                const BoundaryVec3 p =
                    BAdd(start, BScale(displacement, t));
                return BNormSq(BSub(p, ClosestOnTriangle(p, tri)));
            };

            double left = 0.0;
            double right = 1.0;
            for (int i = 0; i < 48; ++i) {
                const double m1 = left + (right - left) / 3.0;
                const double m2 = right - (right - left) / 3.0;
                if (distanceSqAt(m1) < distanceSqAt(m2)) {
                    right = m2;
                } else {
                    left = m1;
                }
            }
            const double minimumT = 0.5 * (left + right);
            if (distanceSqAt(minimumT) > thresholdSq) continue;

            double hitT = 0.0;
            if (distanceSqAt(0.0) > thresholdSq) {
                double lo = 0.0;
                double hi = minimumT;
                for (int i = 0; i < 56; ++i) {
                    const double mid = 0.5 * (lo + hi);
                    if (distanceSqAt(mid) <= thresholdSq) hi = mid;
                    else lo = mid;
                }
                hitT = hi;
            }
            if (hitT >= bestFraction) continue;

            const BoundaryVec3 p =
                BAdd(start, BScale(displacement, hitT));
            const BoundaryVec3 q = ClosestOnTriangle(p, tri);
            contact.found = true;
            contact.triangle = index;
            contact.point = q;
            contact.normal = tri.normal;
            contact.velocity = tri.velocity;
            contact.distance = BNorm(BSub(p, q));
            contact.fraction = hitT;
            bestFraction = hitT;
        }
        return contact;
    }

 private:
    struct Box {
        BoundaryVec3 min{
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity()};
        BoundaryVec3 max{
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity()};
    };
    struct Node {
        Box box;
        std::size_t begin = 0;
        std::size_t end = 0;
        int left = -1;
        int right = -1;
        bool IsLeaf() const { return left < 0; }
    };

    static BoundaryVec3 ClosestOnTriangle(const BoundaryVec3& p,
                                          const BoundaryTriangle& tri) {
        const BoundaryVec3 ab = BSub(tri.b, tri.a);
        const BoundaryVec3 ac = BSub(tri.c, tri.a);
        const BoundaryVec3 ap = BSub(p, tri.a);
        const double d1 = BDot(ab, ap);
        const double d2 = BDot(ac, ap);
        if (d1 <= 0.0 && d2 <= 0.0) return tri.a;

        const BoundaryVec3 bp = BSub(p, tri.b);
        const double d3 = BDot(ab, bp);
        const double d4 = BDot(ac, bp);
        if (d3 >= 0.0 && d4 <= d3) return tri.b;

        const double vc = d1*d4 - d3*d2;
        if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
            return BAdd(tri.a, BScale(ab, d1 / (d1 - d3)));
        }

        const BoundaryVec3 cp = BSub(p, tri.c);
        const double d5 = BDot(ab, cp);
        const double d6 = BDot(ac, cp);
        if (d6 >= 0.0 && d5 <= d6) return tri.c;

        const double vb = d5*d2 - d1*d6;
        if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
            return BAdd(tri.a, BScale(ac, d2 / (d2 - d6)));
        }

        const double va = d3*d6 - d5*d4;
        if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
            const BoundaryVec3 bc = BSub(tri.c, tri.b);
            return BAdd(tri.b, BScale(bc, (d4 - d3) /
                                          ((d4 - d3) + (d5 - d6))));
        }

        const double denom = 1.0 / (va + vb + vc);
        return BAdd(tri.a, BAdd(BScale(ab, vb * denom),
                                BScale(ac, vc * denom)));
    }

    static bool RayTriangle(const BoundaryVec3& origin,
                            const BoundaryVec3& direction,
                            const BoundaryTriangle& tri, double* t) {
        const BoundaryVec3 e1 = BSub(tri.b, tri.a);
        const BoundaryVec3 e2 = BSub(tri.c, tri.a);
        const BoundaryVec3 h = BCross(direction, e2);
        const double det = BDot(e1, h);
        if (std::abs(det) < 1e-13) return false;
        const double invDet = 1.0 / det;
        const BoundaryVec3 s = BSub(origin, tri.a);
        const double u = invDet * BDot(s, h);
        if (u < -1e-12 || u > 1.0 + 1e-12) return false;
        const BoundaryVec3 q = BCross(s, e1);
        const double v = invDet * BDot(direction, q);
        if (v < -1e-12 || u + v > 1.0 + 1e-12) return false;
        *t = invDet * BDot(e2, q);
        return *t >= 0.0;
    }

    static Box TriangleBox(const BoundaryTriangle& tri) {
        Box box;
        for (int d = 0; d < 3; ++d) {
            box.min[d] = std::min({tri.a[d], tri.b[d], tri.c[d]});
            box.max[d] = std::max({tri.a[d], tri.b[d], tri.c[d]});
        }
        return box;
    }
    static void Extend(Box& a, const Box& b) {
        for (int d = 0; d < 3; ++d) {
            a.min[d] = std::min(a.min[d], b.min[d]);
            a.max[d] = std::max(a.max[d], b.max[d]);
        }
    }
    static bool BoxesOverlap(const Box& a, const Box& b) {
        for (int d = 0; d < 3; ++d) {
            if (a.max[d] < b.min[d] || b.max[d] < a.min[d]) return false;
        }
        return true;
    }
    static double BoxDistanceSq(const Box& box, const BoundaryVec3& p) {
        double result = 0.0;
        for (int d = 0; d < 3; ++d) {
            const double delta =
                p[d] < box.min[d] ? box.min[d] - p[d]
              : p[d] > box.max[d] ? p[d] - box.max[d] : 0.0;
            result += delta * delta;
        }
        return result;
    }

    static bool RayHitsBox(const Box& box, const BoundaryVec3& origin,
                           const BoundaryVec3& direction) {
        double tMin = 0.0;
        double tMax = std::numeric_limits<double>::infinity();
        for (int d = 0; d < 3; ++d) {
            if (std::abs(direction[d]) < 1e-15) {
                if (origin[d] < box.min[d] || origin[d] > box.max[d]) {
                    return false;
                }
                continue;
            }
            double t0 = (box.min[d] - origin[d]) / direction[d];
            double t1 = (box.max[d] - origin[d]) / direction[d];
            if (t1 < t0) std::swap(t0, t1);
            tMin = std::max(tMin, t0);
            tMax = std::min(tMax, t1);
            if (tMax < tMin) return false;
        }
        return tMax >= 0.0;
    }

    int BuildNode(std::size_t begin, std::size_t end) {
        Node node;
        node.begin = begin;
        node.end = end;
        for (std::size_t i = begin; i < end; ++i) {
            Extend(node.box, TriangleBox(triangles_[indices_[i]]));
        }
        const int nodeIndex = static_cast<int>(nodes_.size());
        nodes_.push_back(node);
        if (end - begin <= 8) return nodeIndex;

        BoundaryVec3 extent = BSub(node.box.max, node.box.min);
        int axis = extent[1] > extent[0] ? 1 : 0;
        if (extent[2] > extent[axis]) axis = 2;
        const std::size_t middle = begin + (end - begin) / 2;
        std::nth_element(indices_.begin() + begin, indices_.begin() + middle,
                         indices_.begin() + end,
                         [&](std::size_t lhs, std::size_t rhs) {
            const auto centroid = [&](std::size_t index) {
                const auto& t = triangles_[index];
                return (t.a[axis] + t.b[axis] + t.c[axis]) / 3.0;
            };
            return centroid(lhs) < centroid(rhs);
        });
        nodes_[nodeIndex].left = BuildNode(begin, middle);
        nodes_[nodeIndex].right = BuildNode(middle, end);
        return nodeIndex;
    }

    void ClosestNode(int nodeIndex, const BoundaryVec3& p, bool wallsOnly,
                     BoundaryClosestPoint& result) const {
        const Node& node = nodes_[nodeIndex];
        if (BoxDistanceSq(node.box, p) >= result.distance * result.distance) {
            return;
        }
        if (node.IsLeaf()) {
            for (std::size_t i = node.begin; i < node.end; ++i) {
                const std::size_t index = indices_[i];
                const auto& tri = triangles_[index];
                if (wallsOnly && tri.kind != BoundaryKind::kWall) continue;
                const BoundaryVec3 q = ClosestOnTriangle(p, tri);
                const double distance = BNorm(BSub(p, q));
                if (distance < result.distance) {
                    result.found = true;
                    result.triangle = index;
                    result.point = q;
                    result.normal = tri.normal;
                    result.velocity = tri.velocity;
                    result.distance = distance;
                }
            }
            return;
        }
        const double leftDistance = BoxDistanceSq(nodes_[node.left].box, p);
        const double rightDistance = BoxDistanceSq(nodes_[node.right].box, p);
        if (leftDistance < rightDistance) {
            ClosestNode(node.left, p, wallsOnly, result);
            ClosestNode(node.right, p, wallsOnly, result);
        } else {
            ClosestNode(node.right, p, wallsOnly, result);
            ClosestNode(node.left, p, wallsOnly, result);
        }
    }

    std::size_t CountRayIntersections(
        int nodeIndex, const BoundaryVec3& origin,
        const BoundaryVec3& direction) const {
        const Node& node = nodes_[nodeIndex];
        if (!RayHitsBox(node.box, origin, direction)) return 0;
        if (node.IsLeaf()) {
            std::size_t count = 0;
            for (std::size_t i = node.begin; i < node.end; ++i) {
                double t = 0.0;
                if (RayTriangle(origin, direction,
                                triangles_[indices_[i]], &t)
                    && t > 1e-10) {
                    ++count;
                }
            }
            return count;
        }
        return CountRayIntersections(node.left, origin, direction)
             + CountRayIntersections(node.right, origin, direction);
    }

    void CollectOverlappingWalls(
        int nodeIndex, const Box& query,
        std::vector<std::size_t>* result) const {
        const Node& node = nodes_[nodeIndex];
        if (!BoxesOverlap(node.box, query)) return;
        if (node.IsLeaf()) {
            for (std::size_t i = node.begin; i < node.end; ++i) {
                const std::size_t index = indices_[i];
                if (triangles_[index].kind == BoundaryKind::kWall
                    && BoxesOverlap(TriangleBox(triangles_[index]), query)) {
                    result->push_back(index);
                }
            }
            return;
        }
        CollectOverlappingWalls(node.left, query, result);
        CollectOverlappingWalls(node.right, query, result);
    }

    std::vector<BoundaryTriangle> triangles_;
    std::vector<std::size_t> indices_;
    std::vector<Node> nodes_;
    bool closed_ = false;
};

}  // namespace bdm_adapter

#endif  // BDM_ADAPTER_BOUNDARY_GEOMETRY_H_
