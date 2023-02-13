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
#include "helium_ionization.hpp"

HeliumIonization::HeliumIonization(RadiationBand *pband, int my_scalar_number, std::string name, ParameterInput *pin, std::string my_xc_file):
  IonizingAbsorber(pband, my_scalar_number, name, pin)
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
  AthenaArray<Real> file_data;
  ReadDataTable(file_data, xc_file);

  int n_file = file_data.GetDim2();
  Real *file_x = new Real[n_file];
  Real *file_y = new Real[n_file];
  for (int i = 0; i < n_file; ++i)
  {
    file_x[i] = file_data(i,0);
    file_y[i] = file_data(i,1);

    if (i > 0 && file_x[i] < file_x[i-1]) {
      std::stringstream msg;
      msg << "###### FATAL ERROR in HeliumIonization::CalculateCrossSections" << std::endl
          << "Cross Sections file " << xc_file << " must be in ascending order." << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  int ii = -1;
  Real dx;
  Real wave, freq, energy, ry, xc_mb;

  for (int n = 0; n < nspec; ++n)
  { 
    wave = spec[n].wave; // nm
    freq = speed_of_light/(wave * pmy_band->wavelength_coefficient); // Hz
    energy = planck_constant * freq / eV_conversion; // eV
    ry = energy / ry_conversion; // ry

    if (energy < file_x[0]) {
      crossSection(n) = 0.0;
    }
    else if (energy > file_x[n_file-1]) {
      std::stringstream msg;
      msg << "##### FATAL ERROR in HeliumIonization::CalculateCrossSections"
          << std::endl << "Energy too high error." << std::endl;
      ATHENA_ERROR(msg);
    }
    else {
      xc_mb = interp1(ry, file_y, file_x, n_file);
      crossSection(n) = xc_mb * mb_conversion;
    }
  }
}