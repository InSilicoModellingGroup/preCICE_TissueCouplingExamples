#include "BoundaryTriangles.H"

#include "Utilities.H"
#include "polyPatch.H"
#include "wallPolyPatch.H"
#include "volFields.H"

namespace preciceAdapter
{
namespace FA
{

BoundaryTriangles::BoundaryTriangles(
    const Foam::fvMesh& mesh,
    const std::vector<std::string>& patchNames,
    const std::string& dataName)
    : mesh_(mesh)
{
    setDataName(dataName);
    setLocationsType(LocationType::triangleCenters);
    dataType_ = (dataName == "BoundaryKind") ? DataType::scalar
                                               : DataType::vector;

    for (const auto& patchName : patchNames)
    {
        const Foam::label patchI = mesh_.boundaryMesh().findPatchID(patchName);
        if (patchI < 0)
        {
            FatalErrorInFunction
                << "Boundary-triangle patch '" << patchName
                << "' does not exist" << Foam::exit(Foam::FatalError);
        }
        patchIndices_.push_back(patchI);
    }
    setPatchIDs(patchIndices_);
}

std::size_t BoundaryTriangles::write(
    double* buffer, bool meshConnectivity, const unsigned int dim)
{
    if (meshConnectivity)
    {
        adapterInfo("BoundaryTriangles uses triangle records, not preCICE connectivity",
                    "error");
        return 0;
    }

    const std::string field = dataName();
    const bool writeVelocity = field == "BoundaryVelocity";
    const Foam::volVectorField* velocity = nullptr;
    if (writeVelocity && mesh_.foundObject<Foam::volVectorField>("U"))
    {
        velocity = &mesh_.lookupObject<Foam::volVectorField>("U");
    }

    std::size_t out = 0;
    for (const int patchI : patchIndices_)
    {
        const Foam::polyPatch& polyPatch = mesh_.boundaryMesh()[patchI];
        const Foam::List<Foam::face>& faces = polyPatch.localFaces();
        const Foam::pointField& points = polyPatch.localPoints();
        const bool isWall = Foam::isA<Foam::wallPolyPatch>(polyPatch);

        forAll(faces, faceI)
        {
            const Foam::face& f = faces[faceI];
            for (Foam::label i = 1; i + 1 < f.size(); ++i)
            {
                const Foam::point vertices[3] = {
                    points[f[0]], points[f[i]], points[f[i + 1]]};

                if (field == "BoundaryKind")
                {
                    buffer[out++] = isWall ? 1.0 : 0.0;
                    continue;
                }

                Foam::vector value = Foam::vector::zero;
                if (field == "BoundaryTriangleVertex0")
                {
                    value = vertices[0];
                }
                else if (field == "BoundaryTriangleVertex1")
                {
                    value = vertices[1];
                }
                else if (field == "BoundaryTriangleVertex2")
                {
                    value = vertices[2];
                }
                else if (field == "BoundaryNormals")
                {
                    value = (vertices[1] - vertices[0])
                          ^ (vertices[2] - vertices[0]);
                    const Foam::scalar magnitude = Foam::mag(value);
                    if (magnitude > Foam::SMALL)
                    {
                        value /= magnitude;
                    }
                }
                else if (writeVelocity && velocity != nullptr)
                {
                    value = velocity->boundaryField()[patchI][faceI];
                }
                else if (!writeVelocity)
                {
                    FatalErrorInFunction
                        << "Unsupported boundary-triangle field '" << field << "'"
                        << Foam::exit(Foam::FatalError);
                }

                for (unsigned int d = 0; d < dim; ++d)
                {
                    buffer[out++] = value[d];
                }
            }
        }
    }
    return out;
}

bool BoundaryTriangles::isLocationTypeSupported(bool meshConnectivity) const
{
    return !meshConnectivity
        && locationType_ == LocationType::triangleCenters;
}

std::string BoundaryTriangles::getDataName() const
{
    return dataName_;
}

}  // namespace FA
}  // namespace preciceAdapter
