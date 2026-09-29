#ifndef VOLUMESCALARWRITER_H
#define VOLUMESCALARWRITER_H

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

/** Writes one OpenFOAM volume scalar on finite-volume cell centres. */
class VolumeScalarWriter : public CouplingDataUser
{
 public:
  VolumeScalarWriter(const Foam::fvMesh& mesh, const std::string& fieldName);

  std::size_t write(double* dataBuffer, bool meshConnectivity,
                    unsigned int dim) override;
  void read(double*, unsigned int) override {}
  bool isLocationTypeSupported(bool meshConnectivity) const override;
  std::string getDataName() const override;

 private:
  const Foam::fvMesh& mesh_;
  std::string fieldName_;
  const Foam::volScalarField* field_ = nullptr;
};

}  // namespace FA
}  // namespace preciceAdapter

#endif
