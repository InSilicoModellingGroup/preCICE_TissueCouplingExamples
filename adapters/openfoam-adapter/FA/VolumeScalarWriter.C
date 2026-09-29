#include "VolumeScalarWriter.H"

#include "Utilities.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "volFields.H"

namespace preciceAdapter
{
namespace FA
{

VolumeScalarWriter::VolumeScalarWriter(
    const Foam::fvMesh& mesh, const std::string& fieldName)
    : mesh_(mesh), fieldName_(fieldName)
{
  dataType_ = scalar;
  setLocationsType(LocationType::volumeCenters);
  if (mesh_.foundObject<Foam::volScalarField>(fieldName_))
  {
    field_ = &mesh_.lookupObject<Foam::volScalarField>(fieldName_);
  }
  else
  {
    adapterInfo("VolumeScalarWriter: field '" + fieldName_ +
                    "' does not exist",
                "error-deferred");
  }
}

std::size_t VolumeScalarWriter::write(
    double* dataBuffer, bool meshConnectivity, unsigned int)
{
  if (meshConnectivity || field_ == nullptr)
  {
    return 0;
  }

  const auto& values = field_->internalField();
  std::size_t bufferIndex = 0;
  if (cellSetNames_.empty())
  {
    forAll(values, cellI)
    {
      dataBuffer[bufferIndex++] = values[cellI];
    }
  }
  else
  {
    for (const auto& setName : cellSetNames_)
    {
      Foam::cellSet cells(mesh_, setName);
      for (const Foam::label cellI : cells.toc())
      {
        dataBuffer[bufferIndex++] = values[cellI];
      }
    }
  }
  return bufferIndex;
}

bool VolumeScalarWriter::isLocationTypeSupported(bool meshConnectivity) const
{
  return !meshConnectivity && locationType_ == LocationType::volumeCenters;
}

std::string VolumeScalarWriter::getDataName() const { return fieldName_; }

}  // namespace FA
}  // namespace preciceAdapter
