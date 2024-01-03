#ifndef RADIATION_HPP
#define RADIATION_HPP

// C/C++ headers
#include <string>

// Athena++ headers
#include "../athena.hpp"

class MeshBlock;
class ParameterInput;
class Absorber;
class RTSolver;
class Radiation;
class RadiationBand;

struct Spectrum {
  Real flux, wave, wgt;
};

struct float_triplet;

class Radiation {
  friend class RadiationBand;
public:
  // data
  MeshBlock *pmy_block;
  AthenaArray<RadiationBand*> my_bands;

  bool downwards_flag, upwards_flag;

  int nbands;

  // change in energy in a substep due to radiation absorption
  AthenaArray<Real> dE;

  // functions
  Radiation(MeshBlock *pmb, ParameterInput *pin);
  ~Radiation();

  void CalculateRadiativeTransfer(AthenaArray<Real> const& prim, AthenaArray<Real> const& cons_scalar, Real time, int k, int j);
  void CalculateFluxDifference();
  void CalculateEnergyAbsorption(const Real dt);
  void AddRadiationSourceTerm(const Real dt, AthenaArray<Real> &du_hydro);

protected:
  // reserved incoming and outgoing rays
  std::string default_spec_file, default_upwards_spec_file;
  std::string default_rt_solver;

  RadiationScalingFunc UserRadiationScalingFunc;
};

class RadiationBand {
  friend class Absorber;
public:
  std::string my_name;
  Radiation *pmy_rad;
  RTSolver *my_rtsolver;

  AthenaArray<Absorber*> absorbers;
  Real wavelength_coefficient;
  int nabs;

  Real band_scaling_factor;
  int integration_subbins;

  std::string ix1_bc, ox1_bc;

  int nspec;
  Spectrum *spec, *spec_up;
  Real spec_bin_width;

  // band-integrated quantities
  AthenaArray<Real> band_flux;
  AthenaArray<Real> band_tau, band_tau_cell;

  // spectral flux density at cell bottoms
  // units -- [energy / area / time / wavelength]
  // length -- nspec, ncells3, ncells2, ncells1+1
  // flux_down[i] -- Downward facing flux at bottom of cell [i]
  // flux__up[i]   -- Upwards facing flux at bottom of cell [i]
  AthenaArray<Real> flux_down;
  AthenaArray<Real> flux_up;

  // Additional flux generated in cells (due to reactions emitting phtons)
  // units -- [energy / volume / time]
  // length -- nspec, ncells3, ncells2, ncells1
  AthenaArray<Real> emission_coefficient;
  AthenaArray<Real> source_flux_down;
  AthenaArray<Real> source_flux_up;

  // tau_cell -- optical depth inside one cell
  // tau -- optical depth of all cells above this cell including this cell
  //        such that F(i) = F_input * exp[-tau(i)]
  //        note F(i) is flux at bottom of cell (i)
  //        F(iu+1) = F_input
  //        F(iu) = F_input * exp[-tau(iu)]
  AthenaArray<Real> tau, tau_cell;

  RadiationBand(Radiation *prad, std::string name, ParameterInput *pin);
  ~RadiationBand();

  void InitAbsorbers(ParameterInput *pin);
  Absorber* GetAbsorberByName(std::string name, std::string band_name, ParameterInput *pin);

  void SetSpectralProperties(MeshBlock *pmb, AthenaArray<Real> const& w, AthenaArray<Real> const& cons_scalar, int k, int j);
  void RadiativeTransfer(MeshBlock *pmb, Real radiation_scaling, int k, int j);

  int AssignWavelengthToBin(Real wave);
  void ReadFileOntoBand(std::string file, AthenaArray<Real> &output, int ext=2);
protected:
  std::string my_id;

  void ConstructRTSolver(std::string name, ParameterInput *pin);
  void LoadInputSpectrum(std::string file, Spectrum *dest);
  void LoadBandBoundaryConditions(ParameterInput *pin);
};

#endif