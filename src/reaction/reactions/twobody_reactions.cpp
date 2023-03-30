// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header

#include "../../mesh/mesh.hpp"
#include "../../utils/utils.hpp"
#include "../../scalars/scalars.hpp"
#include "../reaction_network.hpp"
#include "twobody_reactions.hpp"


TwoBodyReaction::TwoBodyReaction(std::string name, int r1_num, int r2_num, int p1_num, int p2_num, Real E_reactant, Real E_product):
  Reaction(name)
{
  r1_num_ = r1_num;
  r2_num_ = r2_num;
  p1_num_ = p1_num;
  p2_num_ = p2_num;

  E_reactant_ = E_reactant;
  E_product_  = E_product;
}

void TwoBodyReaction::react(int k, int j, int i) {
  PassiveScalars *ps = pmy_network->pscalars;
  Real T = pmy_network->temperature_(k,j,i);

  Real n_r1 = ps->s(r1_num_, k, j, i)/ps->mass(r1_num_);
  Real n_r2 = ps->s(r2_num_, k, j, i)/ps->mass(r2_num_);

  Real alpha = this->alpha(T, k, j, i);

  Real n_rxn = alpha * n_r1 * n_r2;
  Real e_rxn = (E_reactant_ - E_product_) * n_rxn;

  pmy_network->dn_rate(my_rxn_num, r1_num_, k, j, i) -= n_rxn;
  pmy_network->dn_rate(my_rxn_num, r2_num_, k, j, i) -= n_rxn;
  pmy_network->dn_rate(my_rxn_num, p1_num_, k, j, i) += n_rxn;
  pmy_network->dn_rate(my_rxn_num, p2_num_, k, j, i) += n_rxn;

  pmy_network->de_rate(my_rxn_num, k, j, i) += e_rxn;
}

Real ChargeExchange_HeToH::alpha(Real T, int k, int j, int i) {
  return 1.25E-15 * pow(300.0/T, -0.25);
}

Real ChargeExchange_HToHe::alpha(Real T, int k, int j, int i) {
  Real t1 = 1.75E-11 * pow(300.0/T, 0.75);
  Real t2 = exp(-128000.0/T);
  return t1 * t2;
}

Real CollisionalRelaxation_HeH::alpha(Real T, int k, int j, int i) {
  return 5.0E-10;
}