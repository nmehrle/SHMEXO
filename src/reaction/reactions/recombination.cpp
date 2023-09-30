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

HydrogenicRecombination::HydrogenicRecombination(std::string name, std::vector<int> species, std::vector<Real> stoichiometry, int Z_, std::string case_, ReactionNetwork *pnetwork, AthenaArray<Real> wave, AthenaArray<Real> frac):
  ReactionTemplate(name, species, stoichiometry, pnetwork, wave, frac)
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

Real HeliumRecombination::alpha(Real T, int k, int j, int i) {
  return 1e-11 / (0.00217556 * pow(T, 1.12036526) + 0.32053997 * pow(T, 0.61526097));
}

Real HeliumRecombination::beta(Real T, int k, int j, int i) {
  return -(1e-11 * pmy_network->boltzmann * T) / (0.00218954 * pow(T, 1.18645797) + 0.352778 * pow(T, 0.62645323));
}

Real HeliumTripletRecombination::alpha(Real T, int k, int j, int i) {
  return 1e-11 / (0.00258173 * pow(T, 0.9848205) + 0.10883234 * pow(T, 0.58864659));
}

Real HeliumTripletRecombination::beta(Real T, int k, int j, int i) {
  return -(1e-11 * pmy_network->boltzmann * T) / (0.00427277 * pow(T, 0.99204123) + 0.12332369 * pow(T, 0.57871911));
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

// Real VernerRecombination::beta(Real T, int k, int j, int i) {
//   PassiveScalars *ps = pmy_network->pscalars;

//   Real rxn_energy = ps->energy(ion_num_) - ps->energy(scalar_num_);
//   return rxn_energy * this->alpha(T, k, j, i);
// }