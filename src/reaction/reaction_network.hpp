#ifndef REACTION_NETWORK_HPP
#define REACTION_NETWORK_HPP

// C/C++ headers
#include <string>

// Math headers
#include "../math/eigen335/Eigen/Core"
#include "../math/eigen335/Eigen/LU"

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

  // All have spatial position (k,j,i) as last three elements
  // dn_rate[NSCALARS]            -- Element i contains dn_i/dt summed over all reactions
  // de_rate[num_reactions]       -- Element i contains dE/dt for reaction i
  // jacobian[NSCALARS, NSCALARS] -- Element i,j contains d (dn_i/dt)/ dn_j, summed over reactions
  AthenaArray<Real> dn_rate, de_rate, jacobian;

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
  bool implicit_reactions;
  typedef Eigen::Matrix<Real,NSCALARS,NSCALARS> MatrixNSR;
  typedef Eigen::Matrix<Real,NSCALARS,1> VectorNSR;

  // called at each timestep to rest dn_rate, etc. to be re-calculated
  void ResetRates();

  // calculates dn inside ComputeReactionForcing
  void ComputeScalarDensityChange(const Real dt, Real drho[NSCALARS], int k, int j, int i);
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
  virtual void react(int k, int j, int i) {};

  int my_rxn_num;
protected:
};

#endif