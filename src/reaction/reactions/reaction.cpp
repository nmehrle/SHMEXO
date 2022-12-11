// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "../reaction_network.hpp"

Reaction::Reaction(ReactionNetwork *pnetwork, std::string name):
  pmy_network(pnetwork), my_name(name)
{
  next = nullptr;
  prev = nullptr;
  pnetwork->AddReaction(this);
}