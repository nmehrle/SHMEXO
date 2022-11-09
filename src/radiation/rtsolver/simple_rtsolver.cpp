// Athena++ header
#include "rtsolver.hpp"

SimpleRTSolver::SimpleRTSolver(RadiationBand *pband):
  pmy_band(pband)
{
  return;
}


// computes F = e^-tau
SimpleRTSolver::RadiativeTransfer(int n, int k, int j, int il, int iu)
{
  
}