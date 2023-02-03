// C/C++ headers
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <sstream>    // stringstream
#include <cmath>
#include <cassert>  // assert

// Athena++ headers
#include "../../athena_arrays.hpp"
#include "../../globals.hpp"
#include "../radiation.hpp"
#include "absorber.hpp"
#include "hydrogen_ionization.hpp"

HydrogenIonization::HydrogenIonization(RadiationBand *pband, int my_scalar_number):
  Absorber(pband, my_scalar_number)
{
  Spectrum *spec = pband->spec;
  int nspec = pband->nspec;
  CalculateCrossSections(spec, nspec);
  CalculateEnergyFunctions(spec, nspec);
}

void HydrogenIonization::CalculateCrossSections(Spectrum const *spec, int nspec) {
  Real wave, freq;
  Real eps, term1, numerator, denominator;
  for (int n = 0; n < nspec; ++n)
  {
    wave = spec[n].wave;
    freq = c/(wave * pmy_band->wavelength_coefficient);

    if (freq<nu_0) {
      crossSection(n) = 0.0;
    }
    else if (freq == nu_0) {
      crossSection(n) = A0;
    }
    else {
      eps = sqrt(freq/nu_0 - 1.0);
      term1 = A0 * pow(nu_0/freq, 4.0);
      numerator = exp(4.0 - 4.0*(atan2(eps,1)/eps));
      denominator = 1.0 - exp(-(2*M_PI)/eps);
      crossSection(n) = term1 * numerator/denominator;
    }
  }
}