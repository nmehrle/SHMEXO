#ifndef ABSORBER_HPP_
#define ABSORBER_HPP_

// Athena++ headers
#include "../../athena.hpp"

class RadiationBand;
class Radiation;
class MeshBlock;
class PassiveScalars;
class Reaction;

class Absorber {
  friend class Reaction;
public:
  // data
  RadiationBand *pmy_band;
  PassiveScalars *ps;
  Reaction *pmy_rxn;

  int scalar_num;

  // set by constructor, fixed
  // dims -- [nspec]
  // h -- dictates how much energy goes into phase change
  // q -- dictates how much energy goes into heat
  AthenaArray<Real> crossSection, h, q;

  // updated in time loop
  // dims -- [nspec, ncells3, ncells2, ncells1]
  AthenaArray<Real> absorptionCoefficient;
  AthenaArray<Real> energyAbsorbed;

  Absorber(RadiationBand *pband, int my_scalar_number);
  Absorber(RadiationBand *pband): Absorber(pband, -1);
  virtual ~Absorber();

  // void AssociateReaction(Reaction *prxn);

  // should take in scalars as well
  void CalculateAbsorptionCoefficient(AthenaArray<Real> const& prim, AthenaArray<Real> const& cons_scalar, int n, int k, int j, int i);

  // used in some children
  Real ionization_energy;
protected:
  void SetScalar(int n);

  Real mass;
};

#endif
