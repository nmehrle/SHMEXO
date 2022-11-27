#ifndef RTSOLVER_HPP_
#define RTSOLVER_HPP_

// Athena++ headers
#include "../../athena.hpp"

class RadiationBand;
class Radiation;

class RTSolver {
public:
  RadiationBand *pmy_band;

  RTSolver(RadiationBand *pband, ParameterInput *pin);
  ~RTSolver() {};

  virtual void RadiativeTransfer(MeshBlock *pmb, int n, int k, int j){};
protected:
};

#endif