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
  Real flux, wav, wgt;
};

class Radiation {
  friend class RadiationBand;
public:
  // data
  MeshBlock *pmy_block;
  AthenaArray<RadiationBand*> my_bands;

  int nbands;

  // change in energy in a substep due to radiation absorption
  AthenaArray<Real> dE;

  // functions
  Radiation(MeshBlock *pmb, ParameterInput *pin);
  ~Radiation();

  RadiationBand* GetBand(int n);
  int GetNumBands();

  void CalculateRadiativeTransfer(AthenaArray<Real> const& prim, Real time,
    int k, int j, int il, int iu);
  void CalculateFluxDifference();

protected:
  // reserved incoming and outgoing rays
  std::string default_spec_file, default_rt_solver;

  RadiationScalingFunc UserRadiationScalingFunc;
};

class RadiationBand {
public:
  std::string my_name;
  Radiation *pmy_rad;
  AthenaArray<Absorber*> absorbers;
  int nabs;

  // band-integrated quantities
  AthenaArray<Real> band_flux;
  AthenaArray<Real> band_tau, band_tau_cell;

  RadiationBand(Radiation *prad, std::string name, ParameterInput *pin);
  ~RadiationBand();

  void LoadInputSpectrum(std::string file);
  void ConstructRTSolver(std::string name, ParameterInput *pin);

  Absorber* GetAbsorberByName(std::string name, ParameterInput *pin);

  void SetSpectralProperties(AthenaArray<Real> const& w, int k, int j, int il, int iu);
  void RadiativeTransfer(Real radiation_scaling, int k, int j, int il, int iu);

  void CalculateEnergyAbsorption(const Real dt);
  void AddRadiationSourceTerm(const Real dt, AthenaArray<Real> &du_hydro);
protected:
  std::string my_id;
  int nspec;
  Spectrum *spec;
  RTSolver *my_rtsolver;

  // spectral flux density at bottom of cells
  // units -- [energy / time * area * wavelength]
  // length -- nspec, ncells3, ncells2, ncells1+1
  AthenaArray<Real> spectral_flux_density;

  // tau_cell -- optical depth inside one cell
  // tau -- optical depth of all cells above this cell including this cell
  //        such that F(i) = F_input * exp[-tau(i)]
  //        note F(i) is flux at bottom of cell (i)
  //        F(iu+1) = F_input
  //        F(iu) = F_input * exp[-tau(iu)]
  AthenaArray<Real> tau, tau_cell;
};

#endif