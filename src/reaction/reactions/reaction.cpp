// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "../reaction_network.hpp"

Reaction::Reaction(std::string name):
  my_name(name)
{
  pmy_network = nullptr;
  next = nullptr;
  prev = nullptr;
}