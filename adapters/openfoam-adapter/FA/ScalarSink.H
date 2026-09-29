#ifndef SCALARSINK_H
#define SCALARSINK_H

#include "CouplingDataUser.H"
#include "volFieldsFwd.H"

namespace Foam
{
class fvMesh;
}

namespace preciceAdapter
{
namespace FA
{

/** Reads a non-negative cellular uptake rate into an OpenFOAM scalar field. */
class ScalarSink : public CouplingDataUser
{
 public:
  ScalarSink(Foam::fvMesh& mesh, const std::string& fieldName);

  std::size_t write(double*, bool, unsigned int) override { return 0; }
  void read(double* dataBuffer, unsigned int dim) override;
  bool isLocationTypeSupported(bool meshConnectivity) const override;
  std::string getDataName() const override;

 private:
  Foam::fvMesh& mesh_;
  std::string fieldName_;
  Foam::volScalarField* sink_ = nullptr;
};

}  // namespace FA
}  // namespace preciceAdapter

#endif
