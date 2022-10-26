//! \file disort.cpp
//  \brief SIMPLE radiative transfer solver
//=====================================================

// C/C++ header
#include <cstdlib>
#include <sstream>

// Athena++ headers
#include "../../reconstruct/interpolation.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../../utils/utils.hpp" // StringToArray, ReadTabular
#include "../radiation.hpp"
#include "../../coordinates/coordinates.hpp"

// DISORT headers

#ifdef RT_SIMPLE

void RadiationBand::RadtranRadiance(Direction const rin, Direction const *rout, int nrout,
  Real dist, Real ref_dist, int k, int j, int il, int iu)
{
  /* place holder for calculating radiance
  if (!ds->flag.onlyfl) {
    int count = 0;
    for (int j = 0; j < ds->nphi; ++j)
      for (int k = 0; k < ds->ntau; ++k)
        for (int l = 0; l < ds->numu; ++l, ++count)
          uu(n,j,k,l) = ds_out->uu[count];
  }*/
}

// just do e^-tau
void RadiationBand::RadtranFlux(Direction const rin, Real rad_scaling, int k, int j, int il, int iu)
{
  
  AthenaArray<Real> farea(iu+1);
  pmy_rad->pmy_block->pcoord->Face1Area(k, j, il, iu, farea);

  // Real mu  = rin.mu > 1.E-3 ? rin.mu : 1.E-3;
  // Real phi = rin.phi;

  if (rin.mu != 1.0  || rin.phi !=0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in RadiationBand::RadiativeTransfer (simple): " << std::endl;
    msg << "           <radiation> indir (mu = "<< rin.mu <<", phi = "<< rin.phi <<") must be (1,0)" << std::endl;
    ATHENA_ERROR(msg);
  }
  // bflxup(k,j,iu) = bflxdn(k,j,iu) = 0.;

  // // set flux at top boundary
  // for (int n=0; n< nspec; ++n) {
  //   //stellar source function
  //   Real F0;
  //   if (pmy_rad->GetBeam() < 0.)
  //     F0 = pmy_rad->planet->ParentInsolationFlux(spec[n].wav);
  //   else
  //     F0 = pmy_rad->GetBeam();

  //   // Area ratio account for conical geometry
  //   F0 = F0 * rad_scaling * farea(il)/farea(iu) * spec[n].wgt;

  //   // negative sign indicates downwards transfer
  //   boundary_flux[X1DIR](n,k,j,iu) = -F0;
  //   bflxdn(k,j,iu) += F0;
  // }

  // radiation transfer in column
  Real Ftop, Fbot, attenuation;
  for (int i = iu; i >= il; --i) {
    // reset flux for column
    bflxup(k,j,i) = bflxdn(k,j,i) = 0.;

    // loop over wavelength
    for (int n=0; n< nspec; ++n) {
      if (i==iu) {
        if (pmy_rad->GetBeam() < 0.)
          Fbot = pmy_rad->planet->ParentInsolationFlux(spec[n].wav);
        else
          Fbot = pmy_rad->GetBeam();

        // Area ratio account for conical geometry
        // negative sign indicates downwards transfer
        Fbot = -Fbot * rad_scaling * spec[n].wgt;
      }
      else {
        attenuation = exp(-tau_[i][n]);
        Ftop = boundary_flux[X1DIR](n,k,j,i+1);
        Fbot = Ftop * attenuation;
        net_spectral_flux(n,k,j,i) = (Ftop - Fbot) * farea(i);
      }

      // boundary_flux includes directional sign (negative for down)
      // bflxdn is abs value
      boundary_flux[X1DIR](n,k,j,i) = Fbot;
      bflxdn(k,j,i) -= Fbot;
    }
  }
}

#endif // RT_SIMPLE
