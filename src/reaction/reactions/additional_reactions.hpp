// C/C++ headers
#include <string>

// Athena++ headers
#include "../../athena.hpp"
#include "../reaction_network.hpp"
#include "../reaction.hpp"

class LyaCooling: public ReactionTemplate {
public:
  LyaCooling(std::string name, int scalar_num, int elec_num);
  Real alpha(Real T, int k, int j, int i) {return 0;};
  Real beta(Real T, int k, int j, int i);
};

class Bremsstrahlung: public Reaction {
public:
  Bremsstrahlung(std::string name, int elec_num, int HII_num, int HeII_num, int HeIII_num);
  void react(int k, int j, int i);
protected:
  int elec_num_, HII_num_, HeII_num_, HeIII_num_;
};

class HeliumTripletDecay: public ReactionTemplate {
public:
  using ReactionTemplate::ReactionTemplate;
  Real alpha(Real T, int k, int j, int i);
  Real beta(Real T, int k, int j, int i);
};

class HeIIElecCollisions: public ReactionTemplate {
public:
  HeIIElecCollisions(std::string name, int scalar_num, int elec_num);
  Real alpha(Real T, int k, int j, int i) {return 0;};
  Real beta(Real T, int k, int j, int i);
};

class HeIElecCollisions: public ReactionTemplate {
public:
  HeIElecCollisions(std::string name, int scalar_num, int elec_num);
  Real alpha(Real T, int k, int j, int i);
  Real beta(Real T, int k, int j, int i);
};

class UMISTReaction: public ReactionTemplate {
public:
  UMISTReaction(std::string name, std::vector<int> species, std::vector<Real> stoichiometry, Real alpha_, Real beta_, Real gamma_);
  Real alpha(Real T, int k, int j, int i);
  Real beta(Real T, int k, int j, int i);
protected:
  Real param_alpha, param_beta, param_gamma;
};