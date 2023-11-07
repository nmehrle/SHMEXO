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

  // should take in scalars as well
  virtual Real ReemissionFunction(int b, int n, Real T, int k, int j, int i){};

protected:
  Radiation *prad;
};

class HydrogenReemission: public Reemission {
public:
  HydrogenReemission(ParameterInput *pin, Radiation *prad_, Real ionization_energy_);
  Real ReemissionFunction(int b, int n, Real T, int k, int j, int i);

protected:
  Real ionization_energy;
  Real planck_constant, speed_of_light;
};

#endif