#ifndef SIMPLE_RTSOLVER_HPP
#define SIMPLE_RTSOLVER_HPP

// Athena++ headers
#include "../../athena.hpp"
#include "rtsolver.hpp"

class RadiationBand;
class Radiation;
class RTSolver;

class SimpleRTSolver: public RTSolver {
public:
  SimpleRTSolver(RadiationBand *pband, ParameterInput *pin);
  ~SimpleRTSolver() {};

  void RadiativeTransfer(MeshBlock *pmb, int n, int k, int j);
protected:
};

#endif