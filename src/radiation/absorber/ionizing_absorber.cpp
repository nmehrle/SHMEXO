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
#include "../../parameter_input.hpp"
#include "../../globals.hpp"
#include "../radiation.hpp"
#include "absorber.hpp"
#include "ionizing_absorber.hpp"

IonizingAbsorber::IonizingAbsorber(RadiationBand *pband, int my_scalar_number, std::string name, ParameterInput *pin):
  Absorber(pband, my_scalar_number, pin)
{
  ionization_energy = pin->GetReal("reaction", name+"_ENERGY");
  nu_0 = ionization_energy/planck_constant;
  lambda_0 = speed_of_light/nu_0;

  my_name=name;

  Spectrum *spec = pband->spec;
  int nspec = pband->nspec;

  CalculateEnergyFunctions(spec, nspec);
}

void IonizingAbsorber::CalculateEnergyFunctions(Spectrum const *spec, int nspec) {
  Real wave;

  for (int n = 0; n < nspec; ++n) {
    wave = spec[n].wave * pmy_band->wavelength_coefficient;
    if (wave > lambda_0) {
      // what do we do!!
      // region where cross section is zero so should never be relevant
      // default values
      h(n) = 0.0;
      q(n) = 0.0; 
    }
    else {
      h(n) = wave/lambda_0;
      q(n) = 1.0 - h(n);
    }
  }
}