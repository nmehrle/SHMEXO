// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header

#include "../../mesh/mesh.hpp"
#include "../../utils/utils.hpp"
#include "../../scalars/scalars.hpp"
#include "../reaction_network.hpp"
#include "../reaction.hpp"
#include "twobody_reactions.hpp"

Real ChargeExchangeHeliumToHyd::alpha(Real T, int k, int j, int i) {
  return 1.25E-15 * pow(300.0/T, -0.25);
}

Real ChargeExchangeHydToHelium::alpha(Real T, int k, int j, int i) {
  Real t1 = 1.75E-11 * pow(300.0/T, 0.75);
  Real t2 = exp(-128000.0/T);
  return t1 * t2;
}

// Real ChargeExchangeHeliumToHyd::beta(Real T, int k, int j, int i) {
//   PassiveScalars *ps = pmy_network->pscalars;

//   Real E_reactant = ps->energy(r1_num_) + ps->energy(r2_num_);
//   Real E_product = ps->energy(p1_num_) + ps->energy(p2_num_);

//   return this->alpha(T, k, j, i) * (E_reactant - E_product);
// }

// Real ChargeExchangeHydToHelium::beta(Real T, int k, int j, int i) {
//   PassiveScalars *ps = pmy_network->pscalars;

//   Real E_reactant = ps->energy(r1_num_) + ps->energy(r2_num_);
//   Real E_product = ps->energy(p1_num_) + ps->energy(p2_num_);

//   return this->alpha(T, k, j, i) * (E_reactant - E_product);
// }