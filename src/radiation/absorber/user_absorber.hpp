#ifndef USER_ABSORBER_HPP
#define USER_ABSORBER_HPP

// Athena++ header
#include "absorber.hpp"

using AbsorptionCoefficientFunc = Real (*)(
  Absorber const *me, AthenaArray<Real> const& prim, Real wave, int k, int j, int i);

using EnergyAbsorptionFunc = Real (*)(Absorber *me, Real wave, Real const flux, int k, int j, int i);

class UserDefinedAbsorber: public Absorber {
public:
  UserDefinedAbsorber(RadiationBand *pband);
  virtual ~UserDefinedAbsorber() {};

  Real AbsorptionCoefficient(AthenaArray<Real> const& prim, Real wave, int k, int j, int i);
  Real EnergyAbsorption(Real wave, Real const flux, int k, int j, int i);

  void EnrollUserAbsorptionCoefficientFunc(AbsorptionCoefficientFunc my_func);
  void EnrollUserEnergyAbsorptionFunc(EnergyAbsorptionFunc my_func);

private:
  AbsorptionCoefficientFunc UserAbsorptionCoeffFunc_;
  EnergyAbsorptionFunc UserEnergyAbsorbFunc_;
};

#endif