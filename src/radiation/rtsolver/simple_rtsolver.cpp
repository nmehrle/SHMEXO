// C/C++ headers
#include <sstream>

// Athena++ header
#include "rtsolver.hpp"
#include "simple_rtsolver.hpp"

SimpleRTSolver::SimpleRTSolver(RadiationBand *pband, ParameterInput *pin):
  RTSolver(pband, pin)
{
  return;
}


// computes F = e^-tau
void SimpleRTSolver::RadiativeTransfer(int n, int k, int j, int il, int iu)
{
  
}