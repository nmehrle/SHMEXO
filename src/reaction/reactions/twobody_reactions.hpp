// C/C++ headers
#include <string>

// Athena++ headers
#include "../../athena.hpp"
#include "../reaction_network.hpp"
#include "../reaction.hpp"

class MeshBlock;
class ReactionNetwork;
class Reaction;
class PassiveScalars;

class ChargeExchangeHeliumToHyd: public ReactionTemplate {
public:
  using ReactionTemplate::ReactionTemplate;
  Real alpha(Real T, int k, int j, int i);
  // Real beta(Real T, int k, int j, int i);
};

class ChargeExchangeHydToHelium: public ReactionTemplate {
public:
  using ReactionTemplate::ReactionTemplate;
  Real alpha(Real T, int k, int j, int i);
  // Real beta(Real T, int k, int j, int i);
};