// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "../reaction_network.hpp"
#include "../reaction.hpp"
#include "hydrogen_reactions.hpp"

Real HydrogenRecombination::alpha(Real T, int k, int j, int i) {
  return 2.59E-13 * pow(T/1.E4,-0.7);
}

Real HydrogenRecombination::beta(Real T, int k, int j, int i) {
  return -6.11E-10 * pow(T,-0.89) * (pmy_network->boltzmann * T);
}

LyaCooling::LyaCooling(std::string name, int scalar_num, int ion_num):
    ReactionTemplate(name, {scalar_num, ion_num}, {-1, -1})
{
  ;
}

Real LyaCooling::beta(Real T, int k, int j, int i) {
  return -7.5E-19 * exp(-118348./T);
}