#ifndef HYDROGEN_ABSORBER_HPP
#define HYDROGEN_ABSORBER_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "absorber.hpp"

class HydrogenAbsorber: public Absorber {
public:
  HydrogenAbsorber(RadiationBand *pband):
    Absorber(pband) {}
  HydrogenAbsorber(RadiationBand *pband, ParameterInput *pin);
  virtual ~HydrogenAbsorber() {}
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;
  Real EnergyDeposition(Real wave, Real flux_in, Real flux_out);


protected:
  Real wave_to_meters_conversion;
};

#endif