#include "WallNormals.H"
#include "Utilities.H"
#include "fvMesh.H"
#include "polyPatch.H"
#include "fvPatch.H"
#include "Pstream.H"

using namespace Foam;

namespace preciceAdapter
{
namespace FA
{

WallNormals::WallNormals(const Foam::fvMesh& mesh, const std::vector<std::string>& patchNames)
: mesh_(mesh)
{
    DEBUG(adapterInfo("WallNormals: Constructor called", "debug"));

    // Wall normals use face centres locations
    setLocationsType(LocationType::faceCenters);
    setDataName("WallNormals");
    dataType_ = DataType::vector;

    // Find patch indices from patch names
    for (const auto& patchName : patchNames)
    {
        label patchI = mesh_.boundaryMesh().findPatchID(patchName);
        if (patchI != -1)
        {
            patchIndices_.push_back(patchI);
            DEBUG(adapterInfo("WallNormals: Found patch '" + patchName + "' with index " + std::to_string(patchI), "debug"));
        }
        else
        {
            adapterInfo("WallNormals: WARNING - Patch '" + patchName + "' not found in mesh", "warning");
        }
    }

    if (patchIndices_.empty())
    {
        adapterInfo("WallNormals: ERROR - No valid patches found", "error");
    }
    else
    {
        adapterInfo("WallNormals: Found " + std::to_string(patchIndices_.size()) + " wall patches", "info");
    }

    if (!patchIndices_.empty())
    {
        setPatchIDs(patchIndices_);
    }
}

void WallNormals::initialize()
{
    if (initialized_) return;
    DEBUG(adapterInfo("WallNormals: Initializing wall normal data extraction", "debug"));

    wallFaceNormals_.clear();

    // Reserve storage to avoid repeated reallocations
    std::size_t totalFaces = 0;
    for (const label patchI : patchIndices_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        totalFaces += patch.size();
    }

    wallFaceNormals_.reserve(totalFaces);

    // Set patch IDs and data type
    setPatchIDs(patchIndices_);
    dataType_ = DataType::vector;

    // Extract face normals from all wall patches
    for (const label patchI : patchIndices_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];

        // Skip if patch has no faces on this processor (parallel decomposition)
        if (patch.size() == 0)
        {
            DEBUG(adapterInfo("WallNormals: Skipping patch " + std::to_string(patchI) +
                             " (empty on this processor)", "debug"));
            continue;
        }

        // Get face normals - these should be pre-computed by OpenFOAM
        tmp<vectorField> faceNormalsTmp = patch.nf();
        const vectorField& faceNormals = faceNormalsTmp();

        // CRITICAL: Validate that normal field size matches patch size to prevent memory corruption
        if (faceNormals.size() != patch.size())
        {
            adapterInfo(
                "WallNormals: CRITICAL (rank " + std::to_string(Pstream::myProcNo()) + ") - Normal field size (" +
                std::to_string(faceNormals.size()) + ") does not match patch size (" +
                std::to_string(patch.size()) + ") for patch " + std::to_string(patchI) +
                ". Skipping corrupted patch.",
                "warning");
            continue; // Skip this corrupted patch to avoid crash
        }

        DEBUG(adapterInfo(
            "WallNormals: Processing patch " + std::to_string(patchI) + " with " +
            std::to_string(patch.size()) + " faces on rank " + std::to_string(Pstream::myProcNo()),
            "debug"));

        // Store normal data - safe now that we've validated sizes match
        for (label faceI = 0; faceI < patch.size(); ++faceI)
        {
            wallFaceNormals_.push_back(faceNormals[faceI]);
        }
    }

    adapterInfo("WallNormals: Extracted " + std::to_string(wallFaceNormals_.size()) +
               " wall face normals", "info");

    initialized_ = true;
}

std::size_t WallNormals::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    if (meshConnectivity)
    {
        adapterInfo("WallNormals: meshConnectivity not supported for faceCenters", "warning");
        return 0;
    }

    if (!initialized_)
    {
        initialize();
    }

    DEBUG(adapterInfo("WallNormals: Writing wall face normals to buffer", "debug"));

    // Write normals to buffer (3D vectors)
    std::size_t bufferIndex = 0;
    for (const auto& normal : wallFaceNormals_)
    {
        buffer[bufferIndex++] = normal.x();
        buffer[bufferIndex++] = normal.y();
        buffer[bufferIndex++] = normal.z();
    }

    DEBUG(adapterInfo("WallNormals: Wrote " + std::to_string(wallFaceNormals_.size()) +
                     " face normals (" + std::to_string(bufferIndex) + " components)", "debug"));

    return bufferIndex;
}

bool WallNormals::isLocationTypeSupported(bool meshConnectivity) const
{
    return !meshConnectivity && (locationType_ == LocationType::faceCenters);
}

std::string WallNormals::getDataName() const
{
    return "WallNormals";
}

} // namespace FA
} // namespace preciceAdapter
