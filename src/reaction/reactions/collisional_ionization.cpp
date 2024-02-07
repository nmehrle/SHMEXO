// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "../../scalars/scalars.hpp"
#include "../reaction_network.hpp"
#include "../reaction.hpp"
#include "collisional_ionization.hpp"

CollisionalIonization::CollisionalIonization(std::string name, std::vector<int> species, std::vector<Real> stoichiometry, int mode):
    ReactionTemplate(name, species, stoichiometry)
{
  // set coefficients based on mode
  // mode == 0: H ionization
  // mode == 1: He ionization
  // mode == 2: HeII ionization
  // mode == 3: He23S ionization

  if (mode == 0) {
    alpha_coeff = 5.85E-11;
    beta_coeff = 1.27E-21;
    exponent = -157809.1;
    Tmin = 10000;
  }
  else if (mode == 1) {
    alpha_coeff = 2.38E-11;
    beta_coeff = 9.38E-22;
    exponent = -285335.4;
    Tmin = 10000;
  }
  else if (mode == 2) {
    alpha_coeff = 5.68E-12;
    beta_coeff = 4.95E-22;
    exponent = -631515;
    Tmin = 10000;
  }
  else if (mode == 3) {
    alpha_coeff = 8.38E-10;
    beta_coeff = 6.41E-21;
    exponent = -55338;
    Tmin = 1000;
  }
  else {
    std::stringstream msg;
    msg << "### FATAL ERROR in CollisionalIonization::Constructor" << std::endl
        << "Input Mode (" << mode << ") invalid." <<std::endl
        << "Must be one of 0,1,2,3 (H,He, He+, He23S)." <<std::endl;
    ATHENA_ERROR(msg);
  }
}

Real CollisionalIonization::alpha(Real T, int k, int j, int i) {
  // if (T < Tmin) {
  //   return 0;
  // }
  return alpha_coeff * pow(T, 1/2.) * exp(exponent/T);
}

Real CollisionalIonization::beta(Real T, int k, int j, int i) {
  // if (T < Tmin) {
  //   return 0;
  // }
  return -1. * beta_coeff * pow(T, 1/2.) * exp(exponent/T);
}