#ifndef RADIATION_HPP
#define RADIATION_HPP

// C/C++ headers
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../astronomy/celestrial_body.hpp"

class MeshBlock;
class ParameterInput;
class Absorber;
class Radiation;

#ifdef RT_DISORT
  struct disort_state;
  struct disort_output;
#endif

struct Spectrum {
  Real rad, wav, wgt;
};

struct Direction {
  Real mu, phi;
};

class RadiationBand {
public:
  // data
  std::string myname;
  Radiation *pmy_rad;
  RadiationBand *prev, *next;
  Absorber *pabs;

  // spectra
  Spectrum *spec;
  int nspec, npmom;   // number of spectra and Legendre moments

  // radiation results
  AthenaArray<Real> net_spectral_flux; // energy/time/wavelength
  AthenaArray<Real> boundary_flux[3]; // net flux on boundaries


  // band radiation results
  AthenaArray<Real> btau, bssa, bpmom;
  AthenaArray<Real> bflxup, bflxdn;
  AthenaArray<Real> btoa;

  // functions
  RadiationBand(Radiation *prad); // delayed initialization
  RadiationBand(Radiation *prad, std::string name, ParameterInput *pin);
  ~RadiationBand();
  void AddAbsorber(std::string name, std::string file, ParameterInput *pin);
  void AddAbsorber(Absorber *pab);
  void SetSpectralProperties(AthenaArray<Real> const& w, int k, int j, int il, int iu);
  void RadtranFlux(Direction const rin, Real rad_scaling,
    int k, int j, int il, int iu);
  void RadtranRadiance(Direction const rin, Direction const *rout, int nrout, Real dist, Real ref_dist,
    int k, int j, int il, int iu);
  void CalculateEnergyAbsorption(AthenaArray<Real> &dflx, int k, int j, int il, int iu);

#ifdef RT_DISORT
  void init_disort(ParameterInput *pin);
  void free_disort();
#endif

protected:
  Real **tau_, **ssa_, ***pmom_, *tem_;
  Real **flxup_, **flxdn_;
  Real **toa_;
  Real alpha_;  // T ~ Ts*(\tau/\tau_s)^\alpha at lower boundary

#ifdef RT_DISORT
  disort_state *ds;
  disort_output *ds_out;
#endif
};

class Radiation {
  friend class RadiationBand;
public:
  // data
  MeshBlock *pmy_block;
  RadiationBand *pband;
  Real cooldown, current;
  CelestrialBody *planet;

  AthenaArray<Real> du;

  // functions
  Radiation(MeshBlock *pmb); // delayed initialization
  Radiation(MeshBlock *pmb, ParameterInput *pin);
  ~Radiation();
  RadiationBand* GetBand(int n);
  int GetNumBands();
  void CalculateRadiativeTransfer(AthenaArray<Real> const& w, Real time,
    int k, int j, int il, int iu);
  void CalculateRadiances(AthenaArray<Real> const& w, Real time,
    int k, int j, int il, int iu);
  void CalculateEnergyAbsorption(const Real dt);
  void AddRadiationSourceTerm(const Real dt, AthenaArray<Real> &du_hydro);
  bool IsDynamic() { return dynamic_; }
  Real GetBeam() { return beam_; }

protected:
  // reserved incoming and outgoing rays
  bool dynamic_;
  Direction *rin_, *rout_;
  int nrin_, nrout_;
  Real dist_;
  Real ref_dist_;
  Real beam_;

  RadiationTimeFunc UserRadiationTimeFunc = nullptr;

  // scratch space to compute change in radflux
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_;
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;
  AthenaArray<Real> cell_volume_;
  AthenaArray<Real> dflx_;
};

#endif
