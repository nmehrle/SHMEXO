// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "../reaction_network.hpp"
#include "../reaction.hpp"
#include "additional_reactions.hpp"

LyaCooling::LyaCooling(std::string name, int scalar_num, int ion_num):
    ReactionTemplate(name, {scalar_num, ion_num}, {-1, -1})
{
  ;
}

Real LyaCooling::beta(Real T, int k, int j, int i) {
  return -7.5E-19 * exp(-118348./T);
}

Real HeliumTripletDecay::alpha(Real T, int k, int j, int i) {
  return 1.272E-4; // s^-1 Lampon+ 2020, Table 2
}

// is it zero?
// Energy released as 19eV photon
// What does it do lol
Real HeliumTripletDecay::beta(Real T, int k, int j, int i) {
  // PassiveScalars *ps = pmy_network->pscalars;

  // Real E_reactant = ps->energy(r1_num_);
  // Real E_product = ps->energy(p1_num_);

  // return this->alpha(T, k, j, i) * (E_reactant - E_product);
  return 0;
}