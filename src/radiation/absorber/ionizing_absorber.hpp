#ifndef IONIZING_ABSORBER_HPP
#define IONIZING_ABSORBER_HPP

#include <string>

#include "../../athena.hpp"
#include "absorber.hpp"

class IonizingAbsorber: public Absorber {
public:
  IonizingAbsorber(RadiationBand *pband, std::string name, int my_scalar_number, int my_ion_number, ParameterInput *pin);

  std::string my_name;
  int ion_num;

  Real ionization_energy;
  Real nu_0, lambda_0;

  // overridden by children
  virtual void CalculateCrossSections(Spectrum const *spec, int nspec) {};

  // calculates how much ionizing radiaion becomes heat vs absorbed
  void CalculateEnergyFunctions(Spectrum const *spec, int nspec);
protected:
};
#endif 