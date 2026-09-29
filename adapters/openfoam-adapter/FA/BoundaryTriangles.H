#ifndef BoundaryTriangles_H
#define BoundaryTriangles_H

#include "CouplingDataUser.H"
#include "fvMesh.H"

namespace preciceAdapter
{
namespace FA
{

/**
 * Writes a fixed-size record for every triangle produced by fan-triangulating
 * the configured OpenFOAM boundary faces. The coupling mesh vertex is the
 * triangle centroid; these data fields carry the exact triangle geometry.
 */
class BoundaryTriangles : public CouplingDataUser
{
public:
    BoundaryTriangles(const Foam::fvMesh& mesh,
                      const std::vector<std::string>& patchNames,
                      const std::string& dataName);

    std::size_t write(double* buffer, bool meshConnectivity,
                      const unsigned int dim) override;
    void read(double*, const unsigned int) override {}
    bool isLocationTypeSupported(bool meshConnectivity) const override;
    std::string getDataName() const override;

private:
    const Foam::fvMesh& mesh_;
    std::vector<int> patchIndices_;
};

}  // namespace FA
}  // namespace preciceAdapter

#endif
