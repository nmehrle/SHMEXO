// C/C++ headers
#include <string>
#include <vector>

// Athena++ headers
#include "../../athena.hpp"
#include "../reaction_network.hpp"
#include "../reaction.hpp"

class TripletHydrogenCollision: public ReactionTemplate {
public:
  using ReactionTemplate::ReactionTemplate;
  Real alpha(Real T, int k, int j, int i);
  // Real beta(Real T, int k, int j, int i);
};

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

class HeElectronCollisions: public ReactionTemplate {
public:
  HeElectronCollisions(std::string name, std::vector<int> species, std::vector<Real> stoichiometry, std::string data_file_, int alpha_column_, int beta_column_);

  Real alpha(Real T, int k, int j, int i);
  Real beta(Real T, int k, int j, int i);
protected:
  Real TempScalingFactor(Real T);

  int n_file;
  std::vector<Real> file_T, alpha_data, beta_data;
};