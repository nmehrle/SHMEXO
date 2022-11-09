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

RadiationBand::RadiationBand(Radiation *prad, std::string band_id, ParameterInput *pin):
  my_id(band_id), pmy_rad(prad)
{
  // Gather band name
  my_name = pin->GetString("radiation", my_id);

  char key[100];
  std::string value;
  std::stringstream msg;

  // Gather wavelength details
  sprintf(key, "%s.wavelength", my_id);
  value = pin->GetString("radiation", key);
  std::vector<Real> v = Vectorize(value.c_str());
  if (v.size() != 3) {
    msg << "### FATAL ERROR in RadiationBand::RadiationBand" << std::endl
        << "Length of input file value" << std::endl
        << "<radiation>" << std::endl
        << "'" << my_id << ".wavelength'" << std::endl
        << "must be 3.";
    ATHENA_ERROR(msg);
  }

  nspec = (int)((v[1] - v[0])/v[2]) + 1; // including the last one
  spec = new Spectrum [nspec];
  for (int i = 0; i < nspec; ++i) {
    spec[i].wave = v[0] + v[2]*i;
    spec[i].wgt = (i == 0) || (i == nspec - 1) ? 0.5*v[2] : v[2];
  }
  if (nspec == 1) spec[0].wgt = 1.;

  // Gather input spectrum details
  sprintf(key, "%s.spec_file", my_id);
  value = pin->GetOrAddString("radiation", key, "");
  if (value.empty()) {
    if (pmy_rad->default_spec_file.empty()) {
      msg << "### FATAL ERROR in RadiationBand::RadiationBand" << std::endl
          << "No input spectrum found for band " << my_id << "." << std::endl
          << "Must specify either <radiation> spec_file "
          << "or <radiation> " << my_id << ".spec_file.";
      ATHENA_ERROR(msg);
    }
    value = pmy_rad->default_spec_file;
  }

  LoadInputSpectrum(value);

  // Gather RT Solver Info
  sprintf(key, "%s.rtsolver", my_id);
  value = pin->GetOrAddString("radiation", key, "");
  if (value.empty()) {
    if (pmy_rad->default_rt_solver.empty()) {
      msg << "### FATAL ERROR in RadiationBand::RadiationBand" << std::endl
          << "No input spectrum found for band " << my_id << "." << std::endl
          << "Must specify either <radiation> rtsolver "
          << "or <radiation> " << my_id << ".rtsolver.";
      ATHENA_ERROR(msg);
    }
    value = pmy_rad->default_rt_solver;
  }

  ConstructRTSolver(value, pin);

  // Gather Absorber Info
  int abs_num = 0;
  while(true) {
    sprintf(key, "%s.abs%d", my_id, abs_num);
    try {
      pin->GetString("radiation", key);
    } catch (const std::runtime_error& e) {
      break;
    }
    abs_num++;
  } // while

  nabs = abs_num;
  absorbers.NewAthenaArray(nabs);

  char scalar_key[100];
  int abs_scalar_num = 0;
  std::string abs_name;
  for (int i = 0; i < nabs; ++i)
  {
    sprintf(key, "%s.abs%d", my_id, abs_num);
    abs_name = pin->GetString("radiation", key);
    sprintf(scalar_key, "%s.scalar", key);
    abs_scalar_num = pin->GetOrAddInteger("radiation", scalar_key, -1);

    Absorber* pabs = GetAbsorberByName(abs_name, pin);
    pabs->SetScalar(abs_scalar_num);

    absorbers(i) = pabs;
  }

  MeshBlock *pmb = pmy_rad->pmy_block;
  band_flux.NewAthenaArray(pmb->ncells3, pmb->ncells2, pmb->ncells1+1);
  band_tau.NewAthenaArray(pmb->ncells3, pmb->ncells2, pmb->ncells1);
  band_tau_cell.NewAthenaArray(pmb->ncells3, pmb->ncells2, pmb->ncells1);

  spectral_flux_density.NewAthenaArray(nspec, pmb->ncells3, pmb->ncells2, pmb->ncells1+1);
  tau.NewAthenaArray(nspec, pmb->ncells3, pmb->ncells2, pmb->ncells1);
  tau_cell.NewAthenaArray(nspec, pmb->ncells3, pmb->ncells2, pmb->ncells1);
}

RadiationBand::~RadiationBand() {
  for (int i = 0; i < nabs; ++i)
  {
    delete absorbers(i);
  }
  delete my_rtsolver;
  delete[] spec;
}

void RadiationBand::LoadInputSpectrum(std::string file) {
  AthenaArray<Real> file_data;
  ReadDataTable(file_data, file);

  int n_file = file_data.GetDim2();
  float_triplet *file_spec = new float_triplet[n_file];
  for (int i = 0; i < n_file; ++i)
  {
    file_spec[i].x = file_data(i,0);
    file_spec[i].y = file_data(i,1);
  }

  spline(n_file, file_spec, 0., 0.);

  int ii = -1;
  Real dx;

  for (int n = 0; n < nspec; ++n)
  {
    Real wave = spec[n].wave;
    ii = find_place_in_table(n_file, file_spec, wave, &dx, ii);
    spec[n].flux = splint(wave, file_spec+ii, dx);
  }
}

void RadiationBand::ConstructRTSolver(std::string name, ParameterInput *pin) {
  if (strcmp(name.c_str(), "SIMPLE") == 0) {
    my_rtsolver = new SimpleRTSolver(this, pin);
  }
  else if (strcmp(name.c_str(), "DISORT") == 0) {
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

// To be overridden in pgen file
// Constructs an absorber object based on the name from the problem file
Absorber* __attribute__((weak)) RadiationBand::GetAbsorberByName(
  std::string name, ParameterInput *pin) {}

void RadiationBand::SetSpectralProperties(AthenaArray<Real> const& w, int k, int j, int il, int iu) {
  Coordinates *pcoord = pmy_rad->pmy_block->pcoord;

  tau.ZeroClear();
  tau_cell.ZeroClear();
  band_tau.ZeroClear();
  band_tau_cell.ZeroClear();
  
  // loop collumn
  for (int i = iu; i >= il; --i)
  {
    Real dx = pcoord->dx1f(i);

    // loop wavelength
    for (int n = 0; n < nspec; ++n)
    {

      // loop absorbers
      for (int a = 0; i < nabs; ++a)
      {
        Absorber *pabs = absorbers(a);
        tau_cell(n,k,j,i) += pabs->AbsorptionCoefficient(w, spec[n].wave, k, j, i);
      }
      tau_cell(n,k,j,i) *= dx;

      // calculate tau as sum of tau_cell for this cell and above
      if (i == iu) {
        tau(n,k,j,i) = tau_cell(n,k,j,i);
      }
      else {
        tau(n,k,j,i) = tau_cell(n,k,j,i) + tau(n,k,j,i+1);
      }
      // band values are integrated along wavelength dimension
      band_tau_cell(k,j,i) += tau_cell(n,k,j,i) * spec[n].wgt;
    }

    if (i == iu) {
      band_tau(k,j,i) = band_tau_cell(k,j,i);
    }
    else {
      band_tau(k,j,i) = band_tau_cell(k,j,i) + band_tau(k,j,i+1);
    }
  }
}

// sets spectral_flux_density at the top
// calls my_rtsolver to compute it through the collumn
void RadiationBand::RadiativeTransfer(Real radiation_scaling, int k, int j, int il, int iu) {
  for (int n = 0; n < nspec; ++n) {
    spectral_flux_density(n,k,j,iu+1) = spec[n].flux*radiation_scaling;
    std::cout << "BREAK HERE" << std::endl;
    std::stringstream msg;
    msg << "end of the road, bud";
    ATHENA_ERROR(msg);

    // fills in spectral_flux_density
    my_rtsolver.RadiativeTransfer(int n, int k, int j, int il, int iu);
  }
}