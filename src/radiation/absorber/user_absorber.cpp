// C/C++ headers
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <cassert>  // assert

// Athena++ headers
#include "../athena_arrays.hpp"
#include "../globals.hpp"
#include "absorber.hpp"
#include "user_absorber.hpp"

UserDefinedAbsorber::UserDefinedAbsorber(RadiationBand *pband):
  Absorber(pband),
  UserAbsorptionCoeffFunc_{},
  UserEnergyAbsorbFunc_{}
{
  return;
}

Real UserDefinedAbsorber::AbsorptionCoefficient(Real wave, Real const prim[], int k, int j, int i) const
{
  if (UserAbsorptionCoeffFunc_ != nullptr) {
    return UserAbsorptionCoeffFunc_(this, wave, prim, k, j, i);
  }
  else {
    std::stringstream msg;
    msg << "### FATAL ERROR in UserDefinedAbsorber::AbsorptionCoefficient" << std::endl
        << "Must define an AbsorptionCoefficient funciton in problem file by calling EnrollUserAbsorptionCoefficientFunc()"
        << std::endl;
    ATHENA_ERROR(msg);
  }
}

Real __attribute__((weak)) UserDefinedAbsorber::EnergyAbsorption(Real wave, Real flux, int k, int j, int i)
{
  if (UserEnergyAbsorbFunc_ != nullptr) {
    return UserEnergyAbsorbFunc_(this, wave, flux, k, j, i);
  }
  else {
    std::stringstream msg;
    msg << "### FATAL ERROR in UserDefinedAbsorber::EnergyAbsorption" << std::endl
        << "Must define an EnergyAbsorption funciton in problem file by calling EnrollUserEnergyAbsorptionFunc()"
        << std::endl;
    ATHENA_ERROR(msg);
  }
}

void UserDefinedAbsorber::EnrollUserAbsorptionCoefficientFunc(AbsorptionCoefficientFunc my_func){
  UserAbsorptionCoeffFunc_ = my_func;
}

void UserDefinedAbsorber::EnrollUserEnergyAbsorptionFunc(EnergyAbsorptionFunc my_func){
  UserEnergyAbsorbFunc_ = my_func;
}
