// C/C++ headers
#include <vector>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <type_traits>

// Athena++ header
#include "../math/core.h"
#include "../math/interpolation.h"

#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../utils/utils.hpp"
#include "absorber/absorber.hpp"
#include "rtsolver/rtsolver.hpp"
#include "rtsolver/simple_rtsolver.hpp"
#include "rtsolver/disort_rtsolver.hpp"
#include "radiation.hpp"

void RadiationBand::LoadInputSpectrum(std::string file, Spectrum *dest) {
  AthenaArray<Real> tmp_flx;
  tmp_flx.NewAthenaArray(nspec);
  tmp_flx.ZeroClear();

  if (!file.empty()) {
    ReadFileOntoBand(file, tmp_flx);
  }

  for (int n = 0; n < nspec; ++n)
  {
    dest[n].flux = tmp_flx(n);
  }

  tmp_flx.DeleteAthenaArray();
}

void RadiationBand::ConstructRTSolver(std::string name, ParameterInput *pin) {
  if (strcmp(name.c_str(), "SIMPLE") == 0) {
    delete my_rtsolver;
    my_rtsolver = new SimpleRTSolver(this, pin);
  }
  else if (strcmp(name.c_str(), "DISORT") == 0) {
    delete my_rtsolver;
    my_rtsolver = new DisortRTSolver(this, pin);
  }
  else {
    std::stringstream msg;
    msg << "### FATAL ERROR in RadiationBand::ConstructRTSolver" << std::endl
        << "Radiation Band " << my_id << " found unknown rtsolver " << name
        << std::endl;
    ATHENA_ERROR(msg);
  }
}

void RadiationBand::LoadBandBoundaryConditions(ParameterInput *pin) {
  char key[100];
  std::stringstream msg;

  sprintf(key, "%s.ix1_bc", my_id.c_str());
  ix1_bc = pin->GetOrAddString("radiation", key, "outflow");

  if (strcmp(ix1_bc.c_str(), "outflow") != 0) {
    msg << "### FATAL ERROR in RadiationBand::LoadBandBoundaryConditions" << std::endl
        << "Band " << my_id << " ix1_bc " << ix1_bc << " not valid." << std::endl
        << "Non Outflow Boundary Conditions not yet implemented"
        << std::endl;
    ATHENA_ERROR(msg);
  }

  sprintf(key, "%s.ox1_bc", my_id.c_str());
  ox1_bc = pin->GetOrAddString("radiation", key, "outflow");

  if (strcmp(ox1_bc.c_str(), "outflow") != 0) {
    msg << "### FATAL ERROR in RadiationBand::LoadBandBoundaryConditions" << std::endl
        << "Band " << my_id << " ox1_bc " << ox1_bc << " not valid." << std::endl
        << "Non Outflow Boundary Conditions not yet implemented"
        << std::endl;
    ATHENA_ERROR(msg);
  }
}

int RadiationBand::AssignWavelengthToBin(Real wave) {
  Real half_bin_width = spec_bin_width/2.0;

  Real left = spec[0].wave - half_bin_width;
  Real right;

  // If wave below leftmost edge return -1 (no bin)
  if (wave < left)
    return -1;

  // Loop bins, if wave is in this bin return its value
  for (int i = 0; i < nspec; ++i)
  {
    right = spec[i].wave + half_bin_width;
    if (wave < right)
      return i;
  }

  // If no bin has been returned by now, wave is too high for this band
  return -1;
}

void RadiationBand::ReadFileOntoBand(std::string file, AthenaArray<Real> &output, int ext) {
  std::vector<Real> file_x, file_y;
  int n_file;
  ReadDataTableForInterp(file, file_x, file_y, n_file);
  
  Real half_bin_width = spec_bin_width/2.0;
  Real bin_start;

  Real subbin_dx = spec_bin_width / integration_subbins;
  Real integrated_value, weight, int_x, int_y;

  for (int n = 0; n < nspec; ++n) {
    integrated_value = 0;
    bin_start = spec[n].wave - half_bin_width;

    for (int i = 0; i <= integration_subbins; ++i)
    {
      int_x = bin_start + (i * subbin_dx);

      weight = 2.0;
      if (i == 0 || i == integration_subbins) {
        weight = 1.0;
      }
      int_y = interp1(int_x, file_y.data(), file_x.data(), n_file);

      // Extrapolation conditions
      if (int_x < file_x[0] || int_x > file_x[n_file-1]) {
        if (ext == 0)
          ;
        else if (ext == 1)
          int_y = 0;
        else if (ext == 2) {
          std::stringstream msg;
          msg << "##### FATAL ERROR in RadiationBand::ReadFileOntoBand"
              << "Need to extrapolate when reading file " << file << std::endl
              << "with extrapolate condition set to 2 (error message)" << std::endl
              << "x value " << int_x << " outside file domain ["
              << file_x[0] << ", " << file_x[n_file-1] << "]." << std::endl;
          ATHENA_ERROR(msg);
        }
        else if (ext == 3) {
          if (int_x < file_x[0])
            int_y = file_y[0];
          else if (int_x > file_x[n_file-1])
            int_y = file_y[n_file-1];
        }
      }

      integrated_value += weight * int_y;
    }
    integrated_value = integrated_value * subbin_dx / 2.0;
    output(n) = integrated_value / spec_bin_width;
  }
}