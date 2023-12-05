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

  int nabs = pmy_band->nabs;
  Real Ftop = pmy_band->flux_density(n,k,j,ie+1);

  Real Fsource, Fbot, attenuation, tau_cell, dx;
  Real delta_F, fraction_absorbed, J_i;
  for (int i = ie; i>=is; --i) {
    tau_cell = pmy_band->tau_cell(n,k,j,i);
    attenuation = exp(-tau_cell);

    dx = pmy_band->pmy_rad->pmy_block->pcoord->dx1f(i);
    Fsource = pmy_band->emission_coefficient(n,k,j,i)*dx;

    Fbot = (Ftop) * attenuation;

    // Fsource * (1 - e^(-tau above)) or something
    // half local
    // half with upstream approx
    delta_F = (Ftop - Fbot) + (Fsource * attenuation);

    for (int a = 0; a < nabs; ++a) {
      Absorber *pabs = pmy_band->absorbers(a);
      //energy absorbed by this absorber
      // Ji = (F1-F0)/dx * (a1 * dx / tau)

      fraction_absorbed = pabs->absorptionCoefficient(n,k,j,i)/tau_cell;

      J_i = delta_F * fraction_absorbed;

      pabs->energyAbsorbed(n,k,j,i) = J_i;
    }

    pmy_band->flux_density(n,k,j,i) = Fbot;
    Ftop = Fbot;
  }
}