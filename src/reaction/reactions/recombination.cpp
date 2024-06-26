// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "../../scalars/scalars.hpp"
#include "../reaction_network.hpp"
#include "../reaction.hpp"
#include "recombination.hpp"

HydrogenicRecombination::HydrogenicRecombination(std::string name, std::vector<int> species, std::vector<Real> stoichiometry, int Z_, std::string case_):
  ReactionTemplate(name, species, stoichiometry)
{
  Initialize(Z_, case_);
}

void HydrogenicRecombination::Initialize(int Z_, std::string case_) {
  // Set Z
  Z = Z_;

  // Set Case coefficients following Draine "Physics of the Interstellar and Intergalactic Medium", 2011
  // Coefficients for alpha (C1, C2, C3) taken from equations 14.5, 14.6
  // Coefficients for beta (C1_beta, C2_beta) taken from equations 27.22, 27.23
  if (case_ == "A") {
    C1 = 4.13e-13;
    C2 = -0.7131;
    C3 = -0.0115;
  } else if (case_ == "B") {
    C1 = 2.54e-13;
    C2 = -0.8163;
    C3 = -0.0208;
  } else if (case_ == "1S") {
    C1 = 1.58e-13;
    C2 = -0.5334;
    C3 = -0.0071;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in HydrogenicRecombination::Initialize" << std::endl
        << "Input case (" << case_ << ") must be either \"A\", \"B\", or \"1S\"" <<std::endl;
    ATHENA_ERROR(msg);
  }

  C1_beta = 1.5+C2;
  C2_beta = 2*C3;
}

Real HydrogenicRecombination::alpha(Real T, int k, int j, int i) {
  Real Teff = T/1.e4 / (Z*Z);
  Real exp = C2 + C3 * std::log(Teff);
  return C1 * Z*Z * std::pow(Teff, exp);
}

Real HydrogenicRecombination::beta(Real T, int k, int j, int i) {
  Real Teff = T/1.e4 / (Z*Z);
  Real mean_energy = (C1_beta + C2_beta * std::log(Teff)) * (pmy_network->boltzmann * T);
  return -1 * this->alpha(T, k, j, i) * mean_energy;
}

VernerRecombination::VernerRecombination(std::string name, std::vector<int> species, std::vector<Real> stoichiometry, VernerRecombinationParams *params):
    ReactionTemplate(name, species, stoichiometry)
{
  a = params->a;
  b = params->b;
  T0 = params->T0;
  T1 = params->T1;
}

Real VernerRecombination::alpha(Real T, int k, int j, int i) {
  Real alpha_1 = sqrt(T/T0);
  Real alpha_2 = pow(1+sqrt(T/T0), 1-b);
  Real alpha_3 = pow(1+sqrt(T/T1), 1+b);
  Real alpha = a * pow(alpha_1 * alpha_2 * alpha_3, -1.0); // cm3 s-1
  return alpha;  
}

Real VernerRecombination::beta(Real T, int k, int j, int i) {
  std::stringstream msg;
  msg << "Verner Recombination BETA is incorrect. please update. See Draine chp 27 eq 27.14 and on" << std::endl;
  ATHENA_ERROR(msg);
}

BadnellRecombination::BadnellRecombination(std::string name, std::vector<int> species, std::vector<Real> stoichiometry, Real A_, Real B_, Real T0_, Real T1_):
    ReactionTemplate(name, species, stoichiometry)
{
  A = A_;
  B = B_;
  T0 = T0_;
  T1 = T1_;
}

Real BadnellRecombination::alpha(Real T, int k, int j, int i) {
  Real alpha_1 = sqrt(T/T0);
  Real alpha_2 = pow(1+sqrt(T/T0), 1-B);
  Real alpha_3 = pow(1+sqrt(T/T1), 1+B);
  Real alpha = A * pow(alpha_1 * alpha_2 * alpha_3, -1.0); // cm3 s-1
  return alpha;  
}

Real BadnellRecombination::beta(Real T, int k, int j, int i) {
  Real sqrtT = sqrt(T);
  Real gamma1 = (-1+B) * sqrtT / (2 * (sqrtT + sqrt(T0)));
  Real gamma2 = (-1-B) * sqrtT / (2 * (sqrtT + sqrt(T1)));
  Real gamma = -1 + gamma + gamma2;
  Real mean_energy = (pmy_network->boltzmann * T) * (gamma+2);
  return -1 * this->alpha(T,k,j,i) * mean_energy;
}