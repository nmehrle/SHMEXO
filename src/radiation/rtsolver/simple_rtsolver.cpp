// C/C++ headers
#include <sstream>

// Athena++ header
#include "rtsolver.hpp"
#include "simple_rtsolver.hpp"

#include "../../parameter_input.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../radiation.hpp"
#include "../absorber/absorber.hpp"


SimpleRTSolver::SimpleRTSolver(RadiationBand *pband, ParameterInput *pin):
  RTSolver(pband, pin)
{
  return;
}

// computes F = e^-tau
void SimpleRTSolver::RadiativeTransfer(MeshBlock *pmb, int n, int k, int j)
{
  int is = pmb->is; int ie = pmb->ie;
  RadiationBand *pband = pmy_band;

  // Collect 
  // AthenaArray<Real> farea(ie+1);
  // pmb->pcoord->Face1Area(k, j, is, ie+1, farea);

  int nabs = pband->nabs;
  Real Ftop = pband->flux_density(n,k,j,ie+1);

  Real Fbot, attenuation, tau_cell;
  Real delta_F, fraction_absorbed, J_i;
  for (int i = ie; i>=is; --i) {
    tau_cell = pband->tau_cell(n,k,j,i);
    attenuation = exp(-tau_cell);
    Fbot = Ftop * attenuation;

    for (int a = 0; a < nabs; ++a) {
      Absorber *pabs = pband->absorbers(a);
      //energy absorbed by this absorber
      // Ji = (F1-F0)/dx * (a1 * dx / tau)

      delta_F = Ftop - Fbot;
      fraction_absorbed = pabs->absorptionCoefficient(n,k,j,i)/tau_cell;

      J_i = delta_F * fraction_absorbed;

      pabs->energyAbsorbed(n,k,j,i) = J_i;
    }

    pband->flux_density(n,k,j,i) = Fbot;
    Ftop = Fbot;
  }
}