// C/C++ headers
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <cassert>  // assert

// Athena++ headers
#include "../../athena_arrays.hpp"
#include "../../globals.hpp"
#include "absorber.hpp"
#include "../radiation.hpp"
#include "../../mesh/mesh.hpp"
#include "../../thermodynamics/thermodynamics.hpp"
#include "hydrogen_absorber.hpp"


// Hydrogen photoionization cross section from Osterbrock & Ferland 
// Astrophysics 0f Gaseous Nebulae and Active Galactic Nuclei
// eq 2.4 page 20, pdfpage 34

HydrogenAbsorber::HydrogenAbsorber(RadiationBand *pband, ParameterInput *pin):
  Absorber(pband)
{
  // defaults to nanometers
  wave_to_meters_conversion = pin->GetOrAddReal("radiation","wave_to_meters",1.e-9);
}

Real __attribute__((weak)) HydrogenAbsorber::AbsorptionCoefficient(Real wave, Real const prim[], int k, int j, int i) const
{
  // convert microns to meters
  Real freq = c / (wave * wave_to_meters_conversion);

  // cover low energy and edge case
  if (freq < nu_0)
    return 0.;
  else if (freq == nu_0)
    return A0;

  Real eps   = sqrt(freq/nu_0 - 1.);
  Real term1 = A0 * pow(nu_0/freq,4.);
  
  Real numerator   = exp(4. - 4.*(atan2(eps,1)/eps));
  Real denominator = 1. - exp(-(2*M_PI)/eps);

  Real sigma = term1 * numerator/denominator; // cross-section m^2

  Real p = prim[IPR];
  Real T = prim[IDN];
  const Absorber *pabs = this;
  Thermodynamics * pthermo = pabs->pmy_band->pmy_rad->pmy_block->pthermo;
  Real R = pthermo->GetRd();

  Real dens   = p/(R * T); // kg/m^3

  Real n = dens/mh; // number density

  return sigma * n; // 1/m
}

Real __attribute__((weak)) HydrogenAbsorber::EnergyAbsorption(Real wave, Real flux, int k, int j, int i)
{
  Real wave_nm = (wave * wave_to_meters_conversion) * 1e9;
  if (wave_nm > nm_0) {
    // energy not absorbed
    // should also be reflected in tau_lambda = 0
    return flux;
  }
  else {
    // fraction of energy turned into heat
    // removes ionizatoin energy cost
    Real energy_fraction = 1 - (wave_nm/nm_0);
    return flux*energy_fraction;
  }
}
