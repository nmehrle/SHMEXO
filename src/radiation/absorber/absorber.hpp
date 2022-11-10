#ifndef ABSORBER_HPP_
#define ABSORBER_HPP_

// Athena++ headers
#include "../../athena.hpp"

class RadiationBand;
class Radiation;
class MeshBlock;
class PassiveScalars;

class Absorber {
public:
  // data
  RadiationBand *pmy_band;
  PassiveScalars *ps;

  Absorber(RadiationBand *pband);
  virtual ~Absorber();

  void SetScalar(int n);
  virtual Real AbsorptionCoefficient(AthenaArray<Real> const& prim, Real wave, int k, int j, int i);

protected:
  int scalar_num;
  Real mass;
};

#endif
