#include "ScalarSink.H"

#include "Utilities.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "volFields.H"

#include <cstdint>
#include <cstring>

namespace preciceAdapter
{
namespace FA
{

namespace
{
double sanitize(double value)
{
  std::uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  constexpr std::uint64_t signMask = 0x8000000000000000ULL;
  constexpr std::uint64_t exponentMask = 0x7FF0000000000000ULL;
  if ((bits & exponentMask) == exponentMask || (bits & signMask) != 0ULL)
  {
    return 0.0;
  }
  return value;
}
}  // namespace

ScalarSink::ScalarSink(Foam::fvMesh& mesh, const std::string& fieldName)
    : mesh_(mesh), fieldName_(fieldName)
{
  dataType_ = scalar;
  setLocationsType(LocationType::volumeCenters);
  if (mesh_.foundObject<Foam::volScalarField>(fieldName_))
  {
    sink_ = const_cast<Foam::volScalarField*>(
        &mesh_.lookupObject<Foam::volScalarField>(fieldName_));
  }
  else
  {
    adapterInfo("ScalarSink: field '" + fieldName_ + "' does not exist",
                "error-deferred");
  }
}

void ScalarSink::read(double* dataBuffer, unsigned int)
{
  if (sink_ == nullptr)
  {
    return;
  }

  auto& field = sink_->primitiveFieldRef();
  Foam::label bufferIndex = 0;
  if (cellSetNames_.empty())
  {
    forAll(field, cellI)
    {
      field[cellI] = sanitize(dataBuffer[bufferIndex++]);
    }
  }
  else
  {
    for (const auto& setName : cellSetNames_)
    {
      Foam::cellSet cells(mesh_, setName);
      for (const Foam::label cellI : cells.toc())
      {
        field[cellI] = sanitize(dataBuffer[bufferIndex++]);
      }
    }
  }
  sink_->correctBoundaryConditions();
}

bool ScalarSink::isLocationTypeSupported(bool meshConnectivity) const
{
  return !meshConnectivity && locationType_ == LocationType::volumeCenters;
}

std::string ScalarSink::getDataName() const { return fieldName_; }

}  // namespace FA
}  // namespace preciceAdapter
