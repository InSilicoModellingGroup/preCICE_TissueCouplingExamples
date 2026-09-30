#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "boundary_geometry.h"

using bdm_adapter::BoundaryGeometry;
using bdm_adapter::BoundaryKind;
using bdm_adapter::BoundaryTriangle;
using bdm_adapter::BoundaryVec3;
using bdm_adapter::BoundaryVelocityMode;

namespace {

void AddQuad(std::vector<BoundaryTriangle>* triangles,
             BoundaryVec3 a, BoundaryVec3 b, BoundaryVec3 c, BoundaryVec3 d,
             BoundaryVec3 outward,
             BoundaryKind kind = BoundaryKind::kWall) {
    triangles->push_back({a, b, c, outward, {0, 0, 0}, kind});
    triangles->push_back({a, c, d, outward, {0, 0, 0}, kind});
}

std::vector<BoundaryTriangle> Cube(BoundaryKind rightKind) {
    std::vector<BoundaryTriangle> t;
    AddQuad(&t, {0,0,0}, {0,0,1}, {0,1,1}, {0,1,0}, {-1,0,0});
    AddQuad(&t, {1,0,0}, {1,1,0}, {1,1,1}, {1,0,1}, {1,0,0}, rightKind);
    AddQuad(&t, {0,0,0}, {1,0,0}, {1,0,1}, {0,0,1}, {0,-1,0});
    AddQuad(&t, {0,1,0}, {0,1,1}, {1,1,1}, {1,1,0}, {0,1,0});
    AddQuad(&t, {0,0,0}, {0,1,0}, {1,1,0}, {1,0,0}, {0,0,-1});
    AddQuad(&t, {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}, {0,0,1});
    return t;
}

std::vector<BoundaryTriangle> ConcaveLPrism() {
    const std::vector<BoundaryVec3> lower = {
        {0,0,0}, {2,0,0}, {2,1,0}, {1,1,0}, {1,2,0}, {0,2,0}};
    std::vector<BoundaryTriangle> t;
    for (std::size_t i = 1; i + 1 < lower.size(); ++i) {
        t.push_back({lower[0], lower[i+1], lower[i],
                     {0,0,-1}, {0,0,0}, BoundaryKind::kWall});
        BoundaryVec3 a = lower[0], b = lower[i], c = lower[i+1];
        a[2] = b[2] = c[2] = 1;
        t.push_back({a, b, c, {0,0,1}, {0,0,0}, BoundaryKind::kWall});
    }
    for (std::size_t i = 0; i < lower.size(); ++i) {
        const auto j = (i + 1) % lower.size();
        BoundaryVec3 a = lower[i], b = lower[j], c = lower[j], d = lower[i];
        c[2] = d[2] = 1;
        const BoundaryVec3 edge{
            lower[j][0] - lower[i][0], lower[j][1] - lower[i][1], 0};
        const BoundaryVec3 outward{edge[1], -edge[0], 0};
        AddQuad(&t, a, b, c, d, outward);
    }
    return t;
}

bool Near(double a, double b, double tolerance = 1e-6) {
    return std::abs(a - b) <= tolerance;
}

}  // namespace

int main() {
    BoundaryGeometry cube;
    assert(cube.Build(Cube(BoundaryKind::kWall)));
    assert(cube.Size() == 12);
    assert(cube.IsClosed());
    assert(cube.IsInside({0.5, 0.5, 0.5}));
    assert(!cube.IsInside({1.5, 0.5, 0.5}));
    assert(Near(cube.SignedDistance({0.5, 0.5, 0.5}), 0.5));
    assert(Near(cube.SignedDistance({1.25, 0.5, 0.5}), -0.25));

    const auto closest = cube.ClosestPoint({1.25, 0.5, 0.5}, true);
    assert(closest.found);
    assert(Near(closest.point[0], 1.0));
    assert(Near(closest.normal[0], 1.0));

    const auto contact = cube.FindFirstWallContact(
        {0.5, 0.5, 0.5}, {2.0, 0.5, 0.5}, 0.1);
    assert(contact.found);
    assert(Near(contact.fraction, (0.9 - 0.5) / 1.5, 1e-5));
    assert(Near(contact.normal[0], 1.0));

    const auto edgeContact = cube.FindFirstWallContact(
        {0.5, 0.5, 0.5}, {1.5, 1.5, 0.5}, 0.1);
    assert(edgeContact.found);
    assert(Near(edgeContact.fraction, 0.4, 1e-5));
    const auto cornerContact = cube.FindFirstWallContact(
        {0.5, 0.5, 0.5}, {1.5, 1.5, 1.5}, 0.1);
    assert(cornerContact.found);
    assert(Near(cornerContact.fraction, 0.4, 1e-5));

    const BoundaryVec3 outward{1, 0, 0};
    const BoundaryVec3 wallVelocity{0, 0, 0};
    const auto slip = bdm_adapter::ApplyBoundaryResponse(
        BoundaryVelocityMode::kSlip, {2, 3, 0}, {0.4, 0.6, 0},
        outward, wallVelocity, 0.5, 0.8, 0.0);
    assert(Near(slip.relativeVelocity[0], 0.0));
    assert(Near(slip.relativeVelocity[1], 3.0));
    assert(Near(slip.remainingDisplacement[0], 0.0));

    const auto rebound = bdm_adapter::ApplyBoundaryResponse(
        BoundaryVelocityMode::kRebound, {2, 3, 0}, {0.4, 0.6, 0},
        outward, wallVelocity, 0.5, 0.5, 0.25);
    assert(Near(rebound.relativeVelocity[0], -1.0));
    assert(Near(rebound.relativeVelocity[1], 2.25));
    assert(Near(rebound.remainingDisplacement[0], -0.2));

    const auto stick = bdm_adapter::ApplyBoundaryResponse(
        BoundaryVelocityMode::kStick, {2, 3, 0}, {0.4, 0.6, 0},
        outward, {0, 1, 0}, 0.5, 0.8, 0.0);
    assert(Near(stick.relativeVelocity[0], 0.0));
    assert(Near(stick.remainingDisplacement[1], 0.5));

    const auto positionOnly = bdm_adapter::ApplyBoundaryResponse(
        BoundaryVelocityMode::kNone, {2, 3, 0}, {0.4, 0.6, 0},
        outward, wallVelocity, 0.5, 0.8, 0.0);
    assert(Near(positionOnly.relativeVelocity[0], 2.0));
    assert(Near(positionOnly.remainingDisplacement[0], 0.0));

    BoundaryGeometry openRight;
    assert(openRight.Build(Cube(BoundaryKind::kPassThrough)));
    const auto passThrough = openRight.FindFirstWallContact(
        {0.5, 0.5, 0.5}, {2.0, 0.5, 0.5}, 0.1);
    assert(!passThrough.found);

    BoundaryGeometry concave;
    assert(concave.Build(ConcaveLPrism()));
    assert(concave.IsClosed());
    assert(concave.IsInside({0.5, 1.5, 0.5}));
    assert(!concave.IsInside({1.5, 1.5, 0.5}));
    const auto concaveContact = concave.FindFirstWallContact(
        {0.5, 1.5, 0.5}, {1.5, 1.5, 0.5}, 0.05);
    assert(concaveContact.found);
    assert(Near(concaveContact.point[0], 1.0, 1e-5));

    std::vector<BoundaryTriangle> withDegenerate = Cube(BoundaryKind::kWall);
    withDegenerate.push_back({
        {2,2,2}, {2,2,2}, {2,2,2}, {0,0,0}, {0,0,0},
        BoundaryKind::kWall});
    BoundaryGeometry filtered;
    assert(filtered.Build(std::move(withDegenerate)));
    assert(filtered.Size() == 12);

    std::cout << "boundary_geometry_test passed\n";
    return 0;
}
