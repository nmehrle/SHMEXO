#ifndef REACTION_HPP
#define REACTION_HPP

// C/C++ headers
#include <string>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "reaction_network.hpp"

struct ReemissionCoordinate {
  int band;
  int spec_bin;
  Real frac;
};

class ReactionNetwork;
class PassiveScalars;

// Generic reaction class to be overridden by individual reactions
class Reaction {
public:
  ReactionNetwork *pmy_network;
  std::string my_name;

  Reaction *next;
  Reaction *prev;

  Reaction(std::string name);
  virtual ~Reaction() {};
  virtual void react(int k, int j, int i) {};

  // Reaction Rate Coefficient
  virtual Real alpha(Real T, int k, int j, int i) {};

  // Heating (Cooling) Rate Coefficient
  virtual Real beta(Real T, int k, int j, int i) {};

  int my_rxn_num;
protected:
};

class ReactionTemplate: public Reaction {
public:
  ReactionTemplate(std::string name, std::vector<int> species, std::vector<Real> stoichiometry);
  ReactionTemplate(std::string name, std::vector<int> species, std::vector<Real> stoichiometry, ReactionNetwork *pnetwork, AthenaArray<Real> reemission_wave, AthenaArray<Real> reemission_frac);

  void react(int k, int j, int i);
  void ProducePhotons(Real n_rxn, int k, int j, int i);

  // Heating (Cooling) Rate Coefficient
  virtual Real beta(Real T, int k, int j, int i);

protected:
  int num_species_;
  std::vector<int> species_;
  std::vector<Real> stoichiometry_;

  std::vector<Real> number_density_;

  int nwave_ = 0;
  AthenaArray<ReemissionCoordinate> reemission_coordinates_;
};

class NullReaction: public Reaction {
public:
  NullReaction();
  void react(int k, int j, int i) {};
};

#endif