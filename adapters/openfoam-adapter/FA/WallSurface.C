#include "WallSurface.H"
#include "Utilities.H"
#include "fvMesh.H"
#include "polyPatch.H"
#include "fvPatch.H"

using namespace Foam;

namespace preciceAdapter
{
namespace FA
{

WallSurface::WallSurface(const Foam::fvMesh& mesh, const std::vector<std::string>& patchNames)
: mesh_(mesh)
{
    DEBUG(adapterInfo("WallSurface: Constructor called", "debug"));

    // Wall geometry uses face centres
    setLocationsType(LocationType::faceCenters);
    setDataName("WallGeometry");
    dataType_ = DataType::vector;

    // Find patch indices from patch names
    for (const auto& patchName : patchNames)
    {
        label patchI = mesh_.boundaryMesh().findPatchID(patchName);
        if (patchI != -1)
        {
            patchIndices_.push_back(patchI);
            DEBUG(adapterInfo("WallSurface: Found patch '" + patchName + "' with index " + std::to_string(patchI), "debug"));
        }
        else
        {
            adapterInfo("WallSurface: WARNING - Patch '" + patchName + "' not found in mesh", "warning");
        }
    }

    if (patchIndices_.empty())
    {
        adapterInfo("WallSurface: ERROR - No valid patches found", "error");
    }
    else
    {
        adapterInfo("WallSurface: Found " + std::to_string(patchIndices_.size()) + " wall patches", "info");
    }

    if (!patchIndices_.empty())
    {
        setPatchIDs(patchIndices_);
    }
}

void WallSurface::initialize()
{
    if (initialized_) return;
    DEBUG(adapterInfo("WallSurface: Initializing wall mesh data extraction", "debug"));

    wallFaceCenters_.clear();
    wallFaceNormals_.clear();
    wallFaceAreas_.clear();

    // Reserve storage to avoid repeated reallocations
    std::size_t totalFaces = 0;
    for (const label patchI : patchIndices_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        totalFaces += patch.size();
    }

    wallFaceCenters_.reserve(totalFaces);
    wallFaceNormals_.reserve(totalFaces);
    wallFaceAreas_.reserve(totalFaces);

    // Set patch IDs and data type
    setPatchIDs(patchIndices_);
    dataType_ = DataType::vector;

    // Extract face centers, normals, and areas from all wall patches
    for (const label patchI : patchIndices_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];

        // Skip if patch has no faces on this processor (parallel decomposition)
        if (patch.size() == 0)
        {
            DEBUG(adapterInfo("WallSurface: Skipping patch " + std::to_string(patchI) +
                             " (empty on this processor)", "debug"));
            continue;
        }

        const vectorField& faceCentres = patch.Cf();
        tmp<vectorField> faceNormalsTmp = patch.nf();
        const vectorField& faceNormals = faceNormalsTmp();
        const scalarField& faceAreas = patch.magSf();

        // CRITICAL: Validate that all field sizes match patch size to prevent memory corruption
        if (faceCentres.size() != patch.size() ||
            faceNormals.size() != patch.size() ||
            faceAreas.size() != patch.size())
        {
            adapterInfo("WallSurface: CRITICAL - Field size mismatch for patch " + std::to_string(patchI) +
                       " (patch size: " + std::to_string(patch.size()) +
                       ", centers: " + std::to_string(faceCentres.size()) +
                       ", normals: " + std::to_string(faceNormals.size()) +
                       ", areas: " + std::to_string(faceAreas.size()) + "). Skipping corrupted patch.", "warning");
            continue; // Skip this corrupted patch to avoid crash
        }

        DEBUG(adapterInfo("WallSurface: Processing patch " + std::to_string(patchI) +
                         " with " + std::to_string(patch.size()) + " faces", "debug"));

        // Store face data
        for (label faceI = 0; faceI < patch.size(); ++faceI)
        {
            wallFaceCenters_.push_back(faceCentres[faceI]);
            wallFaceNormals_.push_back(faceNormals[faceI]);
            wallFaceAreas_.push_back(faceAreas[faceI]);
        }
    }

    adapterInfo("WallSurface: Extracted " + std::to_string(wallFaceCenters_.size()) +
               " wall faces total", "info");

    initialized_ = true;
}

std::size_t WallSurface::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    if (meshConnectivity)
    {
        adapterInfo("WallSurface: meshConnectivity not supported for faceCenters", "warning");
        return 0;
    }

    if (!initialized_)
    {
        initialize();
    }


    DEBUG(adapterInfo("WallSurface: Writing wall facecenters to buffer", "debug"));

    // Write facecenters to buffer (3D coordinates)
    std::size_t bufferIndex = 0;
    for (const auto& faceCenter : wallFaceCenters_)
    {
        buffer[bufferIndex++] = faceCenter.x();
        buffer[bufferIndex++] = faceCenter.y();
        buffer[bufferIndex++] = faceCenter.z();
    }

    DEBUG(adapterInfo("WallSurface: Wrote " + std::to_string(wallFaceCenters_.size()) +
                     " face centers (" + std::to_string(bufferIndex) + " components)", "debug"));

    return bufferIndex;
}

bool WallSurface::isLocationTypeSupported(bool meshConnectivity) const
{
    return !meshConnectivity && (locationType_ == LocationType::faceCenters);
}

std::string WallSurface::getDataName() const
{
    return "WallGeometry";
}

} // namespace FA
} // namespace preciceAdapter
