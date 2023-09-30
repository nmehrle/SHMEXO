// C/C++ headers
#include <string>

// Athena++ headers
#include "../../athena.hpp"
#include "../reaction_network.hpp"
#include "../reaction.hpp"

class HydrogenicRecombination: public ReactionTemplate {
public:
  HydrogenicRecombination(std::string name, std::vector<int> species, std::vector<Real> stoichiometry, int Z_, std::string case_);

  void Initialize(int Z_, std::string case_);

  Real alpha(Real T, int k, int j, int i);
  Real beta(Real T, int k, int j, int i);

protected:
  Real Z;
  Real C1, C2, C3;
  Real C1_beta, C2_beta;
};

class HeliumRecombination: public ReactionTemplate {
public:
  using ReactionTemplate::ReactionTemplate;
  Real alpha(Real T, int k, int j, int i);
  Real beta(Real T, int k, int j, int i);
};

class HeliumTripletRecombination: public ReactionTemplate {
public:
  using ReactionTemplate::ReactionTemplate;
  Real alpha(Real T, int k, int j, int i);
  Real beta(Real T, int k, int j, int i);
};

// Recombination following Verner & Ferland 1996
// https://ui.adsabs.harvard.edu/abs/1996ApJS..103..467V
struct VernerRecombinationParams {
  Real a, b, T0, T1;

  VernerRecombinationParams(Real a_, Real b_, Real T0_, Real T1_):
    a(a_), b(b_), T0(T0_), T1(T1_)
  {
  }
};

class VernerRecombination: public ReactionTemplate {
public:
  // params is [a,b,T0, T1]
  VernerRecombination(std::string name, std::vector<int> species, std::vector<Real> stoichiometry, VernerRecombinationParams *params);
  Real alpha(Real T, int k, int j, int i);
  // Real beta(Real T, int k, int j, int i);

protected:
  Real a,b,T0,T1;
};