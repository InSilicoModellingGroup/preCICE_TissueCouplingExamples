#ifndef WallSurface_H
#define WallSurface_H

#include "CouplingDataUser.H"
#include "fvMesh.H"

namespace preciceAdapter
{
namespace FA
{

/**
 * @brief Wall surface data writer that exports wall mesh geometry to preCICE
 *
 * This class extracts wall/boundary mesh data (face centers) from OpenFOAM
 * and makes it available to particle-based solvers through preCICE.
 */
class WallSurface : public CouplingDataUser
{
private:
    const Foam::fvMesh& mesh_;
    std::vector<int> patchIndices_;
    std::vector<Foam::vector> wallFaceCenters_;
    std::vector<Foam::vector> wallFaceNormals_;
    std::vector<Foam::scalar> wallFaceAreas_;

    bool initialized_ = false;

public:
    //- Constructor
    WallSurface(const Foam::fvMesh& mesh, const std::vector<std::string>& patchNames);

    //- Destructor
    virtual ~WallSurface() = default;

    //- Write function (called by preCICE adapter)
    virtual std::size_t write(double* buffer, bool meshConnectivity, const unsigned int dim) override;

    //- Read function (not used for wall surfaces)
    virtual void read(double* buffer, const unsigned int dim) override {}

    //- Location support helper
    bool isLocationTypeSupported(bool meshConnectivity) const override;

    //- Provide data name used when registering writer
    std::string getDataName() const override;

    //- Initialize wall mesh data extraction
    void initialize();

    //- Get wall face centers
    const std::vector<Foam::point>& getWallFaceCenters() const { return wallFaceCenters_; }

    //- Get wall face normals
    const std::vector<Foam::vector>& getWallFaceNormals() const { return wallFaceNormals_; }

    //- Get wall face areas
    const std::vector<Foam::scalar>& getWallFaceAreas() const { return wallFaceAreas_; }
};

} // namespace FA
} // namespace preciceAdapter

#endif // WallSurface_H
