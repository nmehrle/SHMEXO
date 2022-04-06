#ifndef NULL_ABSORBER_HPP
#define NULL_ABSORBER_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "absorber.hpp"

class NullAbsorber: public Absorber {
public:
  NullAbsorber(RadiationBand *pband):
    Absorber(pband) {}
  virtual ~NullAbsorber() {}
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;
};

#endif