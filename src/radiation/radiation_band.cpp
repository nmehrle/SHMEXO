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
  sprintf(key, "%s.wavelength", my_id.c_str());
  value = pin->GetString("radiation", key);
  std::vector<Real> v = Vectorize<Real>(value.c_str());
  if (v.size() != 3) {
    msg << "### FATAL ERROR in RadiationBand::RadiationBand" << std::endl
        << "Length of input file value" << std::endl
        << "<radiation>" << std::endl
        << "'" << my_id << ".wavelength'" << std::endl
        << "must be 3.";
    ATHENA_ERROR(msg);
  }

  Real nspec_fractional = ((v[1] - v[0])/v[2]);
  nspec = (int) nspec_fractional;
  if (nspec_fractional != nspec ) {
    nspec +=1;
  }

  Real bin_width = (v[1]-v[0]) / nspec;

  spec = new Spectrum [nspec];
  for (int i = 0; i < nspec; ++i) {
    spec[i].wave = v[0] + bin_width*(i+1.0/2.0);
    spec[i].wgt = bin_width;
  }

  // Gather wavelength unit
  sprintf(key, "%s.wavelength_coefficient", my_id.c_str());
  wavelength_coefficient = pin->GetReal("radiation", key);

  // Gather input spectrum details
  sprintf(key, "%s.spec_file", my_id.c_str());
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
  sprintf(key, "%s.rtsolver", my_id.c_str());
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

  // Gather band scaling
  sprintf(key, "%s.scaling", my_id.c_str());
  band_scaling_factor = pin->GetOrAddReal("radiation", key, 1.0);

  my_rtsolver = new RTSolver(this, pin);
  ConstructRTSolver(value, pin);

  // Gather Absorber Info
  int abs_num = 0;
  while(true) {
    sprintf(key, "%s.abs%d", my_id.c_str(), abs_num);
    try {
      pin->GetString("radiation", key);
    } catch (const std::runtime_error& e) {
      break;
    }
    abs_num++;
  } // while

  nabs = abs_num;
  absorbers.NewAthenaArray(nabs);

  std::string abs_name;
  for (int i = 0; i < nabs; ++i)
  {
    sprintf(key, "%s.abs%d", my_id.c_str(), i);
    abs_name = pin->GetString("radiation", key);
    // sprintf(scalar_key, "%s.scalar", key);
    // abs_scalar_num = pin->GetOrAddInteger("radiation", scalar_key, -1);

    Absorber* pabs = GetAbsorberByName(abs_name, my_name, pin);

    absorbers(i) = pabs;
  }

  MeshBlock *pmb = pmy_rad->pmy_block;
  band_flux.NewAthenaArray(pmb->ncells3, pmb->ncells2, pmb->ncells1+1);
  band_tau.NewAthenaArray(pmb->ncells3, pmb->ncells2, pmb->ncells1);
  band_tau_cell.NewAthenaArray(pmb->ncells3, pmb->ncells2, pmb->ncells1);

  flux_density.NewAthenaArray(nspec, pmb->ncells3, pmb->ncells2, pmb->ncells1+1);
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

// To be overridden in pgen file
// Constructs an absorber object based on the name from the problem file
Absorber* __attribute__((weak)) RadiationBand::GetAbsorberByName(
  std::string name, std::string band_name, ParameterInput *pin) {}

void RadiationBand::SetSpectralProperties(MeshBlock *pmb, AthenaArray<Real> const& w, AthenaArray<Real> const& cons_scalar, int k, int j) {
  Coordinates *pcoord = pmy_rad->pmy_block->pcoord;

  tau.ZeroClear();
  tau_cell.ZeroClear();
  band_tau.ZeroClear();
  band_tau_cell.ZeroClear();

  int is = pmb->is; int ie = pmb->ie;
  
  // loop collumn
  for (int i = ie; i >= is; --i)
  {
    Real dx = pcoord->dx1f(i);

    // loop wavelength
    for (int n = 0; n < nspec; ++n)
    {

      // loop absorbers
      for (int a = 0; a < nabs; ++a)
      {
        Absorber *pabs = absorbers(a);
        pabs->CalculateAbsorptionCoefficient(w, cons_scalar, n, k, j, i);
        tau_cell(n,k,j,i) += pabs->absorptionCoefficient(n, k, j, i);
        // tau_cell(n,k,j,i) += pabs->AbsorptionCoefficient(w, spec[n].wave, k, j, i);
      }
      tau_cell(n,k,j,i) *= dx;

      // calculate tau as sum of tau_cell for this cell and above
      if (i == ie) {
        tau(n,k,j,i) = tau_cell(n,k,j,i);
      }
      else {
        tau(n,k,j,i) = tau_cell(n,k,j,i) + tau(n,k,j,i+1);
      }
      // band values are integrated along wavelength dimension
      band_tau_cell(k,j,i) += tau_cell(n,k,j,i) * spec[n].wgt;
    }

    if (i == ie) {
      band_tau(k,j,i) = band_tau_cell(k,j,i);
    }
    else {
      band_tau(k,j,i) = band_tau_cell(k,j,i) + band_tau(k,j,i+1);
    }
  }
}

// sets flux_density at the top
// calls my_rtsolver to compute it through the collumn
void RadiationBand::RadiativeTransfer(MeshBlock *pmb, Real radiation_scaling, int k, int j) {
  int ie = pmb->ie;
  for (int n = 0; n < nspec; ++n) {
    flux_density(n,k,j,ie+1) = spec[n].flux * spec[n].wgt * radiation_scaling;
    // flux_density(n,k,j,ie+1) = spec[n].flux * radiation_scaling;

    // fills in flux_density
    // computes which absorber absorbs the radiation
    my_rtsolver->RadiativeTransfer(pmb, n, k, j);
  }
}