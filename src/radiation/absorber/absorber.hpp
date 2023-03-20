#ifndef ABSORBER_HPP_
#define ABSORBER_HPP_

// Athena++ headers
#include "../../athena.hpp"

class ParameterInput;
class RadiationBand;
class Radiation;
class MeshBlock;
class PassiveScalars;
class Reaction;
struct Spectrum;

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

  Absorber(RadiationBand *pband, int my_scalar_number, ParameterInput *pin);
  Absorber(RadiationBand *pband, ParameterInput *pin): Absorber(pband, -1, pin) {};
  virtual ~Absorber();

  // should take in scalars as well
  void CalculateAbsorptionCoefficient(AthenaArray<Real> const& prim, AthenaArray<Real> const& cons_scalar, int n, int k, int j, int i);

  // constants for unit conversion, from pin
  Real speed_of_light, planck_constant;
  Real eV_conversion;// = 1.602E-12; // erg/eV
  Real mb_conversion;// = 1.E-18; // cm^2/Mb
  Real ry_conversion;// = 13.605693123; // eV/Ry

  // calculate cross sections, to be overridden by children
  virtual void CalculateCrossSections(Spectrum const *spec, int nspec){};

  // calculates how much ionizing radiaion becomes heat vs absorbed
  virtual void CalculateEnergyFunctions(Spectrum const *spec, int nspec){};
protected:
  void SetScalar(int n);

  Real mass;
};

#endif
