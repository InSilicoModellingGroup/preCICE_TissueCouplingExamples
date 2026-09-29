#ifndef WallNormals_H
#define WallNormals_H

#include "CouplingDataUser.H"
#include "fvMesh.H"

namespace preciceAdapter
{
namespace FA
{

/**
 * @brief Writer for wall face normal vectors
 *
 * Exports the normal direction of each wall face to preCICE for particle-wall
 * collision detection in particle-based solvers.
 */
class WallNormals : public CouplingDataUser
{
private:
    const Foam::fvMesh& mesh_;
    std::vector<int> patchIndices_;
    std::vector<Foam::vector> wallFaceNormals_;
    bool initialized_ = false;

public:
    //- Constructor
    WallNormals(const Foam::fvMesh& mesh, const std::vector<std::string>& patchNames);

    //- Destructor
    virtual ~WallNormals() = default;

    //- Write function (called by preCICE adapter)
    virtual std::size_t write(double* buffer, bool meshConnectivity, const unsigned int dim) override;

    //- Read function (not used for wall normals)
    virtual void read(double* buffer, const unsigned int dim) override {}

    //- Location support helper
    bool isLocationTypeSupported(bool meshConnectivity) const override;

    //- Provide data name used when registering writer
    std::string getDataName() const override;

    //- Initialize wall normal data extraction
    void initialize();
};

} // namespace FA
} // namespace preciceAdapter

#endif // WallNormals_H
