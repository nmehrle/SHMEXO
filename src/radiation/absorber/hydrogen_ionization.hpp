#ifndef HYDROGEN_IONIZATION_HPP
#define HYDROGEN_IONIZATION_HPP

// Athena++ header
#include "absorber.hpp"
#include "ionizing_absorber.hpp"

class HydrogenIonization: public IonizingAbsorber {
public:
  HydrogenIonization(RadiationBand *pband, int my_scalar_number, std::string name, ParameterInput *pin);
  ~HydrogenIonization() {};

protected:
  // cross section constant
  Real A0=6.30431812E-22; // m2

  void CalculateCrossSections(Spectrum const *spec, int nspec);
};

#endif