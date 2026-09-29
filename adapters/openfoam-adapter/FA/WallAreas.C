#include "WallAreas.H"
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

WallAreas::WallAreas(const Foam::fvMesh& mesh, const std::vector<std::string>& patchNames)
: mesh_(mesh)
{
    DEBUG(adapterInfo("WallAreas: Constructor called", "debug"));

    // Wall areas use face centres locations
    setLocationsType(LocationType::faceCenters);
    setDataName("WallAreas");
    dataType_ = DataType::scalar;

    // Find patch indices from patch names
    for (const auto& patchName : patchNames)
    {
        label patchI = mesh_.boundaryMesh().findPatchID(patchName);
        if (patchI != -1)
        {
            patchIndices_.push_back(patchI);
            DEBUG(adapterInfo("WallAreas: Found patch '" + patchName + "' with index " + std::to_string(patchI), "debug"));
        }
        else
        {
            adapterInfo("WallAreas: WARNING - Patch '" + patchName + "' not found in mesh", "warning");
        }
    }

    if (patchIndices_.empty())
    {
        adapterInfo("WallAreas: ERROR - No valid patches found", "error");
    }
    else
    {
        adapterInfo("WallAreas: Found " + std::to_string(patchIndices_.size()) + " wall patches", "info");
    }

    if (!patchIndices_.empty())
    {
        setPatchIDs(patchIndices_);
    }
}

void WallAreas::initialize()
{
    if (initialized_) return;
    DEBUG(adapterInfo("WallAreas: Initializing wall area data extraction", "debug"));

    wallFaceAreas_.clear();

    // Reserve storage to avoid repeated reallocations
    std::size_t totalFaces = 0;
    for (const label patchI : patchIndices_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        totalFaces += patch.size();
    }

    wallFaceAreas_.reserve(totalFaces);

    // Set patch IDs and data type
    setPatchIDs(patchIndices_);
    dataType_ = DataType::scalar;

    // Extract face areas from all wall patches
    for (const label patchI : patchIndices_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];

        // Skip if patch has no faces on this processor (parallel decomposition)
        if (patch.size() == 0)
        {
            DEBUG(adapterInfo("WallAreas: Skipping patch " + std::to_string(patchI) +
                             " (empty on this processor)", "debug"));
            continue;
        }

        const scalarField& faceAreas = patch.magSf();

        // CRITICAL: Validate that area field size matches patch size to prevent memory corruption
        if (faceAreas.size() != patch.size())
        {
            adapterInfo(
                "WallAreas: CRITICAL (rank " + std::to_string(Pstream::myProcNo()) + ") - Area field size (" +
                std::to_string(faceAreas.size()) + ") does not match patch size (" +
                std::to_string(patch.size()) + ") for patch " + std::to_string(patchI) +
                ". Skipping corrupted patch.",
                "warning");
            continue; // Skip this corrupted patch to avoid crash
        }

        DEBUG(adapterInfo(
            "WallAreas: Processing patch " + std::to_string(patchI) + " with " +
            std::to_string(patch.size()) + " faces on rank " + std::to_string(Pstream::myProcNo()),
            "debug"));

        // Store area data
        for (label faceI = 0; faceI < patch.size(); ++faceI)
        {
            wallFaceAreas_.push_back(faceAreas[faceI]);
        }
    }

    adapterInfo("WallAreas: Extracted " + std::to_string(wallFaceAreas_.size()) +
               " wall face areas", "info");

    initialized_ = true;
}

std::size_t WallAreas::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    if (meshConnectivity)
    {
        adapterInfo("WallAreas: meshConnectivity not supported for faceCenters", "warning");
        return 0;
    }

    if (!initialized_)
    {
        initialize();
    }

    DEBUG(adapterInfo("WallAreas: Writing wall face areas to buffer", "debug"));

    // Write areas to buffer (scalars)
    std::size_t bufferIndex = 0;
    for (const auto& area : wallFaceAreas_)
    {
        buffer[bufferIndex++] = area;
    }

    DEBUG(adapterInfo("WallAreas: Wrote " + std::to_string(wallFaceAreas_.size()) +
                     " face areas", "debug"));

    return bufferIndex;
}

bool WallAreas::isLocationTypeSupported(bool meshConnectivity) const
{
    return !meshConnectivity && (locationType_ == LocationType::faceCenters);
}

std::string WallAreas::getDataName() const
{
    return "WallAreas";
}

} // namespace FA
} // namespace preciceAdapter
