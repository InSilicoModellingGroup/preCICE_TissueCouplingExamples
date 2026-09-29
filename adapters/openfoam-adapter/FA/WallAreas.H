#ifndef WallAreas_H
#define WallAreas_H

#include "CouplingDataUser.H"
#include "fvMesh.H"

namespace preciceAdapter
{
namespace FA
{

/**
 * @brief Writer for wall face areas
 *
 * Exports the area of each wall face to preCICE for area-weighted particle-wall
 * collision influence in particle-based solvers.
 */
class WallAreas : public CouplingDataUser
{
private:
    const Foam::fvMesh& mesh_;
    std::vector<int> patchIndices_;
    std::vector<Foam::scalar> wallFaceAreas_;
    bool initialized_ = false;

public:
    //- Constructor
    WallAreas(const Foam::fvMesh& mesh, const std::vector<std::string>& patchNames);

    //- Destructor
    virtual ~WallAreas() = default;

    //- Write function (called by preCICE adapter)
    virtual std::size_t write(double* buffer, bool meshConnectivity, const unsigned int dim) override;

    //- Read function (not used for wall areas)
    virtual void read(double* buffer, const unsigned int dim) override {}

    //- Location support helper
    bool isLocationTypeSupported(bool meshConnectivity) const override;

    //- Provide data name used when registering writer
    std::string getDataName() const override;

    //- Initialize wall area data extraction
    void initialize();
};

} // namespace FA
} // namespace preciceAdapter

#endif // WallAreas_H
