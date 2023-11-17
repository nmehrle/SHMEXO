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
class ReactionNetwork;
class Reaction;
class Radiation;

class ReactionNetwork {
  friend class Reaction;
public:
  // data
  MeshBlock *pmy_block;
  PassiveScalars *pscalars;
  Radiation *prad;
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

  Real boltzmann, speed_of_light, planck_constant;
  Real eV_conversion, ry_conversion;
  AthenaArray<Real> temperature_;

  // Reads reactions in from pin
  // used in non-radiation network
  void ReadReactionsFromInput(ParameterInput *pin);

  // To be included/overwritten in problem generator
  Reaction* GetReactionByName(std::string name, ParameterInput *pin);
  bool DoRunReactions(int k, int j, int i);

  // to be called at end of problem generator.
  // creates athena arrays for output
  void Initialize();
  void AddReaction(Reaction *prxn);

  void ComputeReactionForcing(const Real dt, const AthenaArray<Real> prim, const AthenaArray<Real> cons, const AthenaArray<Real> cons_scalar, AthenaArray<Real> &du, AthenaArray<Real> &ds);

  int num_reactions;

protected:
  Real consumption_tolerance;
  Real sfloor;

  bool implicit_reactions;
  typedef Eigen::Matrix<Real,NSCALARS,NSCALARS> MatrixNSR;
  typedef Eigen::Matrix<Real,NSCALARS,1> VectorNSR;

  // called at each timestep to rest dn_rate, etc. to be re-calculated
  void ResetRates();

  // calculates dn inside ComputeReactionForcing
  void ComputeScalarDensityChange(const Real dt, Real drho[NSCALARS], int k, int j, int i);
};

#endif