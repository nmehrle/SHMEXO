#ifndef NULL_ABSORBER_HPP
#define NULL_ABSORBER_HPP

// Athena++ header
#include "absorber.hpp"

class NullAbsorber: public Absorber {
public:
  NullAbsorber(RadiationBand *pband):
    Absorber(pband){};
  ~NullAbsorber() {};
  virtual Real AbsorptionCoefficient(AthenaArray<Real> const& prim, Real wave, int k, int j, int i)
    {return 0;};
};

#endif