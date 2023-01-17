// C/C++ headers
#include <sstream>

// Athena++ header
#include "rtsolver.hpp"
#include "disort_rtsolver.hpp"
#include "../../athena.hpp"

DisortRTSolver::DisortRTSolver(RadiationBand *pband, ParameterInput *pin):
  RTSolver(pband, pin)
{
  return;
}

DisortRTSolver::~DisortRTSolver() {
  
}

// computes F = e^-tau
void DisortRTSolver::RadiativeTransfer(MeshBlock*pmb, int n, int k, int j)
{
  std::stringstream msg;
  msg << "##### FATAL ERROR in DisortRTSolver Constructor" << std::endl
      << "DISORT Radiative Transfer Solver is not yet built.";
  ATHENA_ERROR(msg);      
}