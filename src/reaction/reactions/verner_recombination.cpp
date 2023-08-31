// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "../../scalars/scalars.hpp"
#include "../reaction_network.hpp"
#include "../reaction.hpp"
#include "verner_recombination.hpp"

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