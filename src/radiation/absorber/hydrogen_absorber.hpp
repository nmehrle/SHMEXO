#ifndef HYDROGEN_ABSORBER_HPP
#define HYDROGEN_ABSORBER_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "absorber.hpp"

class HydrogenAbsorber: public Absorber {
public:
  Real A0=6.30431812E-22; // (m^2)
  Real nu_0=3.2898419603E15; // (1/s) Rydberg frequency
  Real c=2.998E8; // m/s
  Real mh=1.674E-27; // kg
  Real nm_0=91.126705058; // Hydrogen ionization wavelength in nm

  Real wave_to_meters_conversion;

  HydrogenAbsorber(RadiationBand *pband):
    Absorber(pband) {}
  HydrogenAbsorber(RadiationBand *pband, ParameterInput *pin);
  virtual ~HydrogenAbsorber() {}
  Real AbsorptionCoefficient(Real wave, Real const prim[], int k, int j, int i) const;
  Real EnergyAbsorption(Real wave, Real flux, int k, int j, int i);


protected:
};

#endif