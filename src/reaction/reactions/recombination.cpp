// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "../../scalars/scalars.hpp"
#include "../reaction_network.hpp"
#include "../reaction.hpp"
#include "recombination.hpp"

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

Real HydrogenRecombination::alpha(Real T, int k, int j, int i) {
  return 2.59E-13 * pow(T/1.E4,-0.7);
}

Real HydrogenRecombination::beta(Real T, int k, int j, int i) {
  return -6.11E-10 * pow(T,-0.89) * (pmy_network->boltzmann * T);
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
