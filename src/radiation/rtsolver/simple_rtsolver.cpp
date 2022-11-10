// Athena++ header
#include "rtsolver.hpp"

SimpleRTSolver::SimpleRTSolver(RadiationBand *pband, ParameterInput *pin):
  RTSolver(pband, pin)
{
  return;
}


// computes F = e^-tau
void SimpleRTSolver::RadiativeTransfer(int n, int k, int j, int il, int iu)
{
  
}