// C/C++ headers
#include <sstream>

// Athena++ headers
#include "../athena.hpp"
#include "../utils/utils.hpp"
#include "../math/core.h"
#include "../math/interpolation.h"

#include "../parameter_input.hpp"
#include "../radiation/radiation.hpp"
#include "reemission.hpp"

HeliumReemisison::HeliumReemisison(ParameterInput *pin, Radiation *prad_, std::string decay_21S_file, Real branch_continuum_, Real branch_21S_, Real branch_21P_):
  Reemission(pin, prad_)
{
  energy_21S = pin->GetReal("reaction", "HE_21S_ENERGY") * eV_conversion;
  energy_21P = pin->GetReal("reaction", "HE_21P_ENERGY") * eV_conversion;

  SetBranchingRatios(branch_continuum_, branch_21S_, branch_21P_);
  AssignDecayFile(decay_21S_file);
}

void HeliumReemisison::SetBranchingRatios(Real branch_continuum_, Real branch_21S_, Real branch_21P_)
{
  // Fraction to ground directly from 2(1)S
  branch_21S = branch_21S_;

  // Fraction from 2(1)P can go to ground or 2(1)S
  branch_21S += branch_21P_ * (state21P_to_21S);
  branch_21P  = branch_21P_ * (1-state21P_to_21S);

  // Fraction from continuum is split to 2(1)S/2(1)P
  // Can go to 2(1)S directly
  branch_21S += (branch_continuum_) * continuum_to_21S;

  // Can go to 2(1)P then either ground or 2(1)S
  branch_21S += (branch_continuum_) * (1-continuum_to_21S) * (state21P_to_21S);
  branch_21P += (branch_continuum_) * (1-continuum_to_21S) * (1-state21P_to_21S);
}

Real HeliumReemisison::ReemissionFunction(int b, int n, Real T, int k, int j, int i)
{
  // T dependent branching
  if (using_branching_file) {
    Real branch_21S_ = interp1(T, branch_21S_array.data(), branching_T.data(), n_branching_file);
    Real branch_21P_ = interp1(T, branch_21P_array.data(), branching_T.data(), n_branching_file);
    Real branch_continuum_ = 1 - branch_21S_ - branch_21P_;
    SetBranchingRatios(branch_continuum_, branch_21S_, branch_21P_);
  }

  Real from_21S = branch_21S * Reemission21S(b, n); 
  Real from_21P = branch_21P * Reemission21P(b, n);

  return from_21S+from_21P;
}

Real HeliumReemisison::Reemission21S(int b, int n) {
  return decay_21S_fraction(b,n);
}

Real HeliumReemisison::Reemission21P(int b, int n) {
  Real critical_wave = (planck_constant * speed_of_light)/energy_21P;
  return FixedWaveReemission(b, n, critical_wave);
}

void HeliumReemisison::AssignDecayFile(std::string decay_21S_file) {
  int max_band_length=0;
  int this_band_length;

  for (int b = 0; b < prad->nbands; ++b)
  {
    this_band_length = prad->my_bands(b)->nspec;
    if (this_band_length > max_band_length)
      max_band_length = this_band_length;
  }

  decay_21S_fraction.NewAthenaArray(prad->nbands, max_band_length);

  for (int b = 0; b < prad->nbands; ++b)
  {
    this_band_length = prad->my_bands(b)->nspec;
    AthenaArray<Real> tmp_flx;

    tmp_flx.NewAthenaArray(this_band_length);
    prad->my_bands(b)->ReadFileOntoBand(decay_21S_file, tmp_flx, 1);

    for (int n = 0; n < this_band_length; ++n)
    {
      decay_21S_fraction(b,n) = tmp_flx(n);
    }

    tmp_flx.DeleteAthenaArray();
  }
}

void HeliumReemisison::AssignBranchingFile(std::string branching_file, int column_21S, int column_21P) {
  ReadDataTableForInterp(branching_file, branching_T, branch_21S_array, n_branching_file, 0, column_21S, true);
  ReadDataTableForInterp(branching_file, branching_T, branch_21P_array, n_branching_file, 0, column_21P, true);
  using_branching_file = true;
}
