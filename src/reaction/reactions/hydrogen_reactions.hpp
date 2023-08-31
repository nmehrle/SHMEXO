// C/C++ headers
#include <string>

// Athena++ headers
#include "../../athena.hpp"
#include "../reaction_network.hpp"
#include "../reaction.hpp"

class HydrogenRecombination: public ReactionTemplate {
public:
  using ReactionTemplate::ReactionTemplate;
  Real alpha(Real T, int k, int j, int i);
  Real beta(Real T, int k, int j, int i);
};

class LyaCooling: public ReactionTemplate {
public:
  LyaCooling(std::string name, int scalar_num, int ion_num);
  Real alpha(Real T, int k, int j, int i) {return 0;};
  Real beta(Real T, int k, int j, int i);
};