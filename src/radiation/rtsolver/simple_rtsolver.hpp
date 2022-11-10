#ifndef SIMPLE_RTSOLVER_HPP
#define SIMPLE_RTSOLVER_HPP

// Athena++ headers
#include "../../athena.hpp"

class RadiationBand;
class Radiation;

class SimpleRTSolver: public RTSolver {
public:
  RadiationBand *pmy_band;

  SimpleRTSolver(RadiationBand *pband, ParameterInput *pin):
    RTSolver(pband, pin) {};
  ~SimpleRTSolver();

  void RadiativeTransfer(int n, int k, int j, int il, int iu);
protected:
};

#endif