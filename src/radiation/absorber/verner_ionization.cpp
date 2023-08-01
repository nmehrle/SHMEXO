// C/C++ headers
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <sstream>    // stringstream
#include <cmath>
#include <cassert>  // assert

// Athena++ headers
#include "../../math/core.h"
#include "../../math/interpolation.h"
#include "../../athena_arrays.hpp"
#include "../../globals.hpp"
#include "../../utils/utils.hpp"
#include "../../parameter_input.hpp"
#include "../radiation.hpp"
#include "absorber.hpp"
#include "ionizing_absorber.hpp"
#include "verner_ionization.hpp"

VernerIonization::VernerIonization(RadiationBand *pband, int my_scalar_number, std::string name, ParameterInput *pin, Real Eth_in, Real Emax_in, Real E0_in, Real sigma0_in, Real ya_in, Real P_in, Real yw_in, Real y0_in, Real y1_in):
  IonizingAbsorber(pband, my_scalar_number, name, pin)
{
  Spectrum *spec = pband->spec;
  int nspec = pband->nspec;

  Eth = Eth_in;
  Emax = Emax_in;
  E0 = E0_in;
  sigma0 = sigma0_in;
  ya = ya_in;
  P = P_in;
  yw = yw_in;
  y0 = y0_in;
  y1 = y1_in;

  CalculateCrossSections(spec, nspec);
}

void VernerIonization::CalculateCrossSections(Spectrum const *spec, int nspec) {
  Real wave, freq, energy;
  Real x, y, Fy1, Fy2, Fy3;
  


  for (int n = 0; n < nspec; ++n)
  {
    wave = spec[n].wave;
    freq = speed_of_light/(wave * pmy_band->wavelength_coefficient);
    energy = planck_constant * freq / eV_conversion; // eV

    if (energy < Eth) {
      crossSection(n) = 0.0;
    }
    else if (energy > Emax) {
      // energy exceeds maximum ionization energy for our formula
      // throw error and request user find correct formula

      std::stringstream msg;
      msg << "###### FATAL ERROR in VernerIonization::CalculateCrossSections." << std::endl
          << " desired energy (" << energy << " eV) exceeds maximum energy"
          << " (" << Emax << " eV) for Verner Absorber " << my_name << std::endl
          << " Find new valid cross sections for these higher energies to proceed." <<std::endl;
      ATHENA_ERROR(msg);
    }
    else {
      // valid energy
      x = energy/E0 - y0;
      y = sqrt(x*x + y1*y1);
      Fy1 = (x-1)*(x-1) + yw*yw;
      Fy2 = pow(y, 0.5*P - 5.5);
      Fy3 = pow(1 + sqrt(y/ya), -P);
      crossSection(n) = sigma0 * Fy1 * Fy2 * Fy3 * mb_conversion;
    }
  }
}