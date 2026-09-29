#ifndef FA_H
#define FA_H

#include <string>
#include <vector>

namespace Foam
{
class dictionary;
class fvMesh;
}

namespace preciceAdapter
{
class Interface;

namespace FA
{

/**
 * Tissue-report subset of the fluid-agent module.
 *
 * The module deliberately exposes only the coupling data used by the four
 * tissue examples: extracellular H2O2/NO2 concentrations, their scalar sink
 * feedback, and static boundary geometry. Particle momentum feedback is not
 * part of these examples and is intentionally absent.
 */
class FluidParticle
{
 public:
  explicit FluidParticle(const Foam::fvMesh& mesh);

  bool configure(const Foam::dictionary& adapterConfig);
  bool addWriters(std::string dataName, Interface* interface);
  bool addReaders(std::string dataName, Interface* interface);

  bool addWallSurfaceWriter(const std::string& meshName,
                            const std::vector<std::string>& patchNames,
                            Interface* interface);
  bool addBoundaryTriangleWriters(
      const std::vector<std::string>& patchNames,
      Interface* interface);

 private:
  const Foam::fvMesh& mesh_;
  std::string nameH2O2_ = "H2O2";
  std::string nameNO2_ = "NO2";
  bool enableH2O2Field_ = true;
  bool enableNO2Field_ = true;
  bool isConfigured_ = false;
};

}  // namespace FA
}  // namespace preciceAdapter

#endif
