// Athena++ header
#include "rtsolver.hpp"
#include "disort_rtsolver.hpp"

DisortRTSolver::DisortRTSolver(RadiationBand *pband, ParameterInput *pin):
  RTSolver(pband, pin)
{
  return;
}

// computes F = e^-tau
void DisortRTSolver::RadiativeTransfer(int n, int k, int j, int il, int iu)
{
  std::stringstream msg;
  msg << "##### FATAL ERROR in DisortRTSolver Constructor" << std::endl
      << "DISORT Radiative Transfer Solver is not yet built.";
  ATHENA_ERROR(msg);      
}