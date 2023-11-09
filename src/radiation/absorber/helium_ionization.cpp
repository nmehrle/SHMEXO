// C/C++ headers
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <sstream>    // stringstream
#include <cmath>
#include <vector>
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
#include "helium_ionization.hpp"

HeliumIonization::HeliumIonization(RadiationBand *pband, std::string name, int my_scalar_number, int my_ion_number, ParameterInput *pin, std::string my_xc_file):
  IonizingAbsorber(pband, name, my_scalar_number, my_ion_number, pin)
{
  Spectrum *spec = pband->spec;
  int nspec = pband->nspec;

  xc_file = my_xc_file;

  if (!FileExists(xc_file)) {
    std::stringstream msg;
    msg << "##### FATAL ERROR in HeliumIonization Constructor"
        << std::endl << "Cannot open cross sections file "
        << xc_file << std::endl;
    ATHENA_ERROR(msg);
  }

  CalculateCrossSections(spec, nspec);
}

void HeliumIonization::CalculateCrossSections(Spectrum const *spec, int nspec) {
  std::vector<Real> file_x, file_y;
  int n_file;
  ReadDataTableForInterp(xc_file, file_x, file_y, n_file, true);

  Real wave, freq, energy, ry, xc_mb;

  for (int n = 0; n < nspec; ++n)
  { 
    wave = spec[n].wave; // nm
    freq = speed_of_light/(wave * pmy_band->wavelength_coefficient); // Hz
    energy = planck_constant * freq / eV_conversion; // eV
    ry = energy / ry_conversion; // ry

    if (ry < file_x[0]) {
      crossSection(n) = 0.0;
    }
    else if (ry > file_x[n_file-1]) {
      std::stringstream msg;
      msg << "##### FATAL ERROR in HeliumIonization::CalculateCrossSections"
          << std::endl << "Energy " << energy << "eV (" << ry << " Ry) " <<std::endl
          << "exceeds maximum energy in cross_section file."
          << " (" << file_x[n_file-1] << ")" << std::endl
          << "For absorber " << my_name << std::endl;
      ATHENA_ERROR(msg);
    }
    else {
      xc_mb = interp1(ry, file_y.data(), file_x.data(), n_file);
      crossSection(n) = xc_mb * mb_conversion;
    }
  }
}