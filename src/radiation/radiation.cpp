// C/C++ headers
#include <sstream>
#include <stdexcept>

// Athena++ header
//#include "../math_funcs.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
#include "radiation.hpp"

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin):
  pmy_block(pmb), dE(pmb->ncells3, pmb->ncells2, pmb->ncells1)
{
  UserRadiationScalingFunc = pmy_block->pmy_mesh->UserRadiationScalingFunc_;
  default_spec_file = pin->GetOrAddString("radiation", "spec_file","");
  default_rt_solver = pin->GetOrAddString("radiation", "rtsolver", "SIMPLE");

  // Count number of Radiation Bands
  int band_num = 0;
  char band_id[80];

  while(true) {
    sprintf(band_id, "b%d", band_num);
    try {
      pin->GetString("radiation",band_id);
    } catch (const std::runtime_error& e) {
      break;
    }
    band_num++;
  } // while

  nbands = band_num;
  my_bands.NewAthenaArray(nbands);

  for (int i = 0; i<nbands; ++i)
  {
    sprintf(band_id, "b%d", i);
    RadiationBand *p = new RadiationBand(this, band_id, pin);
    my_bands(i) = p;
  }
}

Radiation::~Radiation() {
  for (int i = 0; i < nbands; ++i)
  {
    delete my_bands(i);
  }
}

void Radiation::CalculateRadiativeTransfer(AthenaArray<Real> const& prim, AthenaArray<Real> const& cons_scalar, Real time, int k, int j)
{

  for (int i = 0; i < nbands; ++i)
  {
    RadiationBand *p = my_bands(i);
    Real radiation_scaling = UserRadiationScalingFunc(p, prim, time, k, j);

    p->SetSpectralProperties(pmy_block, prim, cons_scalar, k, j);

    p->RadiativeTransfer(pmy_block, radiation_scaling, k, j);
  }
}