// C/C++ headers
#include <string>

// Athena++ headers
#include "../../athena.hpp"
#include "../reaction_network.hpp"

class MeshBlock;
class ReactionNetwork;
class Reaction;
class PassiveScalars;

class TwoBodyReaction: public Reaction {
public:
  TwoBodyReaction(std::string name, int r1_num, int r2_num, int p1_num, int p2_num, Real E_reactant, Real E_product);
  void react(int k, int j, int i);

  virtual Real alpha(Real T, int k, int j, int i) {};

protected:
  int r1_num_;
  int r2_num_;
  int p1_num_;
  int p2_num_;

  Real E_reactant_;
  Real E_product_;
};

class ChargeExchange_HeToH: public TwoBodyReaction {
public:
  // ChargeExchange_HeToH(ReactionNetwork *pnetwork, std::string name, int HEplus_num, int H_num, int HE_num, int Hplus_num, Real E_reactant, Real E_product): TwoBodyReaction(pnetwork, name, HEplus_num, H_num, HE_num, Hplus_num, E_reactant, E_product) {};
  using TwoBodyReaction::TwoBodyReaction;

  Real alpha(Real T, int k, int j, int i);
};

class ChargeExchange_HToHe: public TwoBodyReaction {
public:
  // ChargeExchange_HToHe(ReactionNetwork *pnetwork, std::string name, int HE_num, int Hplus_num, int HEplus_num, int H_num, Real E_reactant, Real E_product);

  using TwoBodyReaction::TwoBodyReaction;
  Real alpha(Real T, int k, int j, int i);
};

class CollisionalRelaxation_HeH: public TwoBodyReaction {
public:
  // CollisionalRelaxation_HeH(ReactionNetwork *pnetwork, std::string name, int HEtrip_num, int H_num, int HE_num, Real E_reactant, Real E_product);

  using TwoBodyReaction::TwoBodyReaction;
  Real alpha(Real T, int k, int j, int i);
};