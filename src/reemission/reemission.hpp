#ifndef REEMISSION_HPP_
#define REEMISSION_HPP_

// Athena++ headers
#include <vector>
#include "../athena.hpp"

class ParameterInput;
class Radiation;
class MeshBlock;

class Reemission {  
public:
  Reemission(ParameterInput *pin, Radiation *prad_);
  virtual ~Reemission();

  virtual Real ReemissionFunction(int b, int n, Real T, int k, int j, int i){};
  Real FixedWaveReemission(int b, int n, Real wave);

protected:
  Radiation *prad;
  Real planck_constant, speed_of_light;
  Real eV_conversion;
};

class GroundStateReemission: public Reemission {
public:
  GroundStateReemission(ParameterInput *pin, Radiation *prad_, Real ionization_energy_);
  Real ReemissionFunction(int b, int n, Real T, int k, int j, int i);

protected:
  Real ionization_energy;
};

class HeliumReemisison: public Reemission {
public:
  HeliumReemisison(ParameterInput *pin, Radiation *prad_, std::string decay_21S_file, Real branch_continuum_, Real branch_21S_, Real branch_21P_);
  Real ReemissionFunction(int b, int n, Real T, int k, int j, int i);

  void AssignDecayFile(std::string decay_21S_file);
  void AssignBranchingFile(std::string branching_file, int column_21S, int column_21P);

protected:
  void SetBranchingRatios(Real branch_continuum_, Real branch_21S_, Real branch_21P_);
  Real Reemission21S(int b, int n);
  Real Reemission21P(int b, int n);

  AthenaArray<Real> decay_21S_fraction;
  Real energy_21P, energy_21S;

  Real branch_21S, branch_21P;

  Real continuum_to_21S = 1./3;
  Real state21P_to_21S  = 1.098E-3;

  bool using_branching_file = false;
  int n_branching_file;
  std::vector<Real> branching_T, branch_21S_array, branch_21P_array;
};

#endif