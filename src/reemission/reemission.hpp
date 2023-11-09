#ifndef REEMISSION_HPP_
#define REEMISSION_HPP_

// Athena++ headers
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
  Real Reemission21S(int b, int n);
  Real Reemission21P(int b, int n);

protected:
  AthenaArray<Real> decay_21S_fraction;

  Real branch_21S, branch_21P;
  Real continuum_to_21S, continuum_to_21P, state21P_to_21S;

  Real energy_21P, energy_21S;
};

#endif