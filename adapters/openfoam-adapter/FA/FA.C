#include "FA.H"

#include "BoundaryTriangles.H"
#include "VolumeScalarWriter.H"
#include "Interface.H"
#include "ScalarSink.H"
#include "Utilities.H"
#include "WallAreas.H"
#include "WallNormals.H"
#include "WallSurface.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"

namespace preciceAdapter
{
namespace FA
{

namespace
{
const Foam::dictionary* findModuleDictionary(
    const Foam::dictionary& adapterConfig)
{
  const Foam::dictionary* module = adapterConfig.findDict("FA");
  return module ? module : adapterConfig.findDict("FP");
}
}  // namespace

FluidParticle::FluidParticle(const Foam::fvMesh& mesh) : mesh_(mesh) {}

bool FluidParticle::configure(const Foam::dictionary& adapterConfig)
{
  if (const Foam::dictionary* module = findModuleDictionary(adapterConfig))
  {
    nameH2O2_ = module->lookupOrDefault<Foam::word>("nameH2O2", "H2O2");
    nameNO2_ = module->lookupOrDefault<Foam::word>("nameNO2", "NO2");
    enableH2O2Field_ = module->lookupOrDefault<bool>("enableH2O2", true);
    enableNO2Field_ = module->lookupOrDefault<bool>("enableNO2", true);
  }

  if (enableH2O2Field_ &&
      !mesh_.foundObject<Foam::volScalarField>(nameH2O2_))
  {
    adapterInfo("FA: configured H2O2 field '" + nameH2O2_ +
                    "' does not exist",
                "error-deferred");
    return false;
  }
  if (enableNO2Field_ && !mesh_.foundObject<Foam::volScalarField>(nameNO2_))
  {
    adapterInfo("FA: configured NO2 field '" + nameNO2_ +
                    "' does not exist",
                "error-deferred");
    return false;
  }

  isConfigured_ = true;
  adapterInfo("FA: configured tissue scalar transport subset", "info");
  return true;
}

bool FluidParticle::addWriters(std::string dataName, Interface* interface)
{
  if (!isConfigured_ || interface == nullptr)
  {
    return false;
  }

  std::string fieldName;
  if (dataName == "H2O2" && enableH2O2Field_)
  {
    fieldName = nameH2O2_;
  }
  else if (dataName == "NO2" && enableNO2Field_)
  {
    fieldName = nameNO2_;
  }
  else
  {
    return false;
  }

  auto* writer = new VolumeScalarWriter(mesh_, fieldName);
  writer->setLocationsType(LocationType::volumeCenters);
  interface->addCouplingDataWriter(dataName, writer);
  return true;
}

bool FluidParticle::addReaders(std::string dataName, Interface* interface)
{
  if (!isConfigured_ || interface == nullptr ||
      (dataName != "H2O2Sink" && dataName != "NO2Sink"))
  {
    return false;
  }

  auto* reader = new ScalarSink(const_cast<Foam::fvMesh&>(mesh_), dataName);
  reader->setLocationsType(LocationType::volumeCenters);
  interface->addCouplingDataReader(dataName, reader);
  return true;
}

bool FluidParticle::addWallSurfaceWriter(
    const std::string&,
    const std::vector<std::string>& patchNames,
    Interface* interface)
{
  if (!isConfigured_ || interface == nullptr)
  {
    return false;
  }

  interface->addCouplingDataWriter(
      "WallGeometry", new WallSurface(mesh_, patchNames));
  interface->addCouplingDataWriter(
      "WallNormals", new WallNormals(mesh_, patchNames));
  interface->addCouplingDataWriter(
      "WallAreas", new WallAreas(mesh_, patchNames));
  return true;
}

bool FluidParticle::addBoundaryTriangleWriters(
    const std::vector<std::string>& patchNames,
    Interface* interface)
{
  if (!isConfigured_ || interface == nullptr)
  {
    return false;
  }

  const std::vector<std::string> fields = {
      "BoundaryTriangleVertex0", "BoundaryTriangleVertex1",
      "BoundaryTriangleVertex2", "BoundaryKind", "BoundaryNormals",
      "BoundaryVelocity"};
  for (const auto& field : fields)
  {
    auto* writer = new BoundaryTriangles(mesh_, patchNames, field);
    writer->setLocationsType(LocationType::triangleCenters);
    interface->addCouplingDataWriter(field, writer);
  }
  return true;
}

}  // namespace FA
}  // namespace preciceAdapter
