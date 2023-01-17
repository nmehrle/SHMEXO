#ifndef REACTION_NETWORK_HPP
#define REACTION_NETWORK_HPP

// C/C++ headers
#include <string>

// Athena++ headers
#include "../athena.hpp"

class MeshBlock;
class ParameterInput;
class PassiveScalars;
class Absorber;
class ReactionNetwork;
class Reaction;

class ReactionNetwork {
  friend class Reaction;
public:
  // data
  MeshBlock *pmy_block;
  PassiveScalars *pscalars;
  Reaction *pfirst, *plast;

  AthenaArray<Reaction*> my_reactions;

  // int stoichiometry_matrix;
  AthenaArray<Real> dn_rate, de_rate;

  // functions
  ReactionNetwork(MeshBlock *pmb, ParameterInput *pin);
  ~ReactionNetwork();

  Real boltzmann;
  AthenaArray<Real> temperature_;

  // To be included in problem generator
  // Adds user reactions
  void __attribute__((weak)) InitUserReactions(ParameterInput *pin);

  // to be called at end of problem generator.
  // creates athena arrays for output
  void Initialize();
  void AddReaction(Reaction *prxn);

  void ComputeReactionForcing(const Real dt, const AthenaArray<Real> prim, const AthenaArray<Real> cons, const AthenaArray<Real> cons_scalar, AthenaArray<Real> &du, AthenaArray<Real> &ds);

  int num_reactions;
protected:
};

// general type to be overridden by individual reactions
class Reaction {
public:
  ReactionNetwork *pmy_network;
  std::string my_name;

  Reaction *next;
  Reaction *prev;

  Reaction(ReactionNetwork *pnetwork, std::string name);
  virtual ~Reaction() {};
  virtual void react(AthenaArray<Real> &dn_rate, AthenaArray<Real> &de_rate, int k, int j, int i) {};

  int my_rxn_num;
protected:
};

#endif