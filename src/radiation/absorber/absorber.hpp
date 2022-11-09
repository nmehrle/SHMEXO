#ifndef ABSORBER_HPP_
#define ABSORBER_HPP_

// Athena++ headers
#include "../../athena.hpp"

class RadiationBand;

class Absorber {
public:
  // data
  RadiationBand *pmy_band;
  PassiveScalars *ps;

  Absorber(RadiationBand *pband): pmy_band(pband), ps(pmy_band->pmy_rad->pmy_block->pscalars) {};
  ~Absorber() {}

  void SetScalar(int n) {scalar_num = n;};
  virtual Real AbsorptionCoefficient(AthenaArray<Real> const& prim, Real wave, int k, int j, int i);

protected:
  int scalar_num;
};

#endif
