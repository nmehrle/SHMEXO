// C/C++ headers
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <cassert>  // assert

// Athena++ headers
#include "../athena_arrays.hpp"
#include "../globals.hpp"
#include "absorber.hpp"
#include "radiation.hpp"
#include "../mesh/mesh.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "hydrogen_absorber.hpp"


// Hydrogen photoionization cross section from Osterbrock & Ferland 
// Astrophysics 0f Gaseous Nebulae and Active Galactic Nuclei
// eq 2.4 page 20, pdfpage 34

Real A0=6.30431812E-22; // (m^2)
Real nu_0=3.2898419603E15; // (1/s) Rydberg frequency
Real c=2.998E8; // m/s
Real mh=1.674E-27; // kg 

HydrogenAbsorber::HydrogenAbsorber(RadiationBand *pband, ParameterInput *pin):
  Absorber(pband)
{
  // defaults to nanometers
  wave_to_meters_conversion = pin->GetOrAddReal("radiation","wave_to_meters",1e-9);
}

Real HydrogenAbsorber::AbsorptionCoefficient(Real wave, Real const prim[]) const
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
