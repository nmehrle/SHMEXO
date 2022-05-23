#ifndef USER_ABSORBER_HPP
#define USER_ABSORBER_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "absorber.hpp"

using AbsorptionCoefficientFunc = void (*)(
  Real wave, Real const prim[], int k, int j, int i) const;

using EnergyAbsorptionFunc = void (*)(Real wave, Real const flux, int k, int j, int i);

class UserDefinedAbsorber: public Absorber {
public:
  UserDefinedAbsorber(RadiationBand *pband);
  virtual ~UserDefinedAbsorber() {}

  Real AbsorptionCoefficient(Real wave, Real const prim[], int k, int j, int i) const;
  Real EnergyAbsorption(Real wave, Real const flux, int k, int j, int i);

  void EnrollUserAbsorptionCoefficientFunc(AbsorptionCoefficientFunc my_func);
  void EnrollUserEnergyAbsorptionFunc(EnergyAbsorptionFunc my_func);

private:
  AbsorptionCoefficientFunc UserAbsorptionCoeffFunc_;
  EnergyAbsorptionFunc UserEnergyAbsorbFunc_;
};

#endif