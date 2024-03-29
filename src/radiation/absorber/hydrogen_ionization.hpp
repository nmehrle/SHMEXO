#ifndef HYDROGEN_IONIZATION_HPP
#define HYDROGEN_IONIZATION_HPP

// Athena++ header
#include "absorber.hpp"
#include "ionizing_absorber.hpp"

class HydrogenIonization: public IonizingAbsorber {
public:
  HydrogenIonization(RadiationBand *pband, std::string name, int my_scalar_number, int my_ion_number, ParameterInput *pin);
  ~HydrogenIonization() {};

protected:
  // cross section constant
  Real A0=6.30431812E-18; // cm2

  void CalculateCrossSections(Spectrum const *spec, int nspec);
};

#endif