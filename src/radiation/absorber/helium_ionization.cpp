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
#include "ionizing_absorber.hpp"
#include "helium_ionization.hpp"

HeliumIonization::HeliumIonization(RadiationBand *pband, int my_scalar_number, std::string name, ParameterInput *pin):
  IonizingAbsorber(pband, my_scalar_number, name, pin)
{
  Spectrum *spec = pband->spec;
  int nspec = pband->nspec;
  CalculateCrossSections(spec, nspec);
}

void HeliumIonization::CalculateCrossSections(Spectrum const *spec, int nspec) {
  Real wave, freq, energy;
  Real x, y, F1, F2, F3;

  for (int n = 0; n < nspec; ++n)
  {
    wave = spec[n].wave; // AA
    freq = speed_of_light/(wave * pmy_band->wavelength_coefficient); // Hz
    energy = planck_constant * freq / eV_conversion; // eV

    if (energy < Eth) {
      crossSection(n) = 0.0;
    }
    else if (energy > Emax) {
      // energy exceeds maximum ionization energy for our formula
      // throw error and request person find correct formula

      std::stringstream msg;
      msg << "###### FATAL ERROR IN HeliumIonization::CalculateCrossSections." << std::endl
          << " desired energy (" << energy << " eV) exceeds maximum energy for helium"
          << " photoionization formula from Verner1996 (" << Emax << " eV)." << std::endl
          << " Find new valid cross sections for these higher energies to proceed." <<std::endl;
      ATHENA_ERROR(msg);
    }
    else {
      // valid energy
      x = energy/E0 - y0;
      y = sqrt(x*x + y1*y1);
      F1 = (x-1)*(x-1) + yw*yw;
      F2 = pow(y, 0.5*P - 5.5);
      F3 = pow(1 + sqrt(y/ya), -P);
      crossSection(n) = sig0 * F1 * F2 * F3 * mb_conversion;
    }
  }
}