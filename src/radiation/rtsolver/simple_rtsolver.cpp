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
  RadiationBand *pband = pmy_band;
  int is = pmb->is; int ie = pmb->ie;

  Real tau_cell, dx, emission_coefficient;
  Real Fin, Fout;
  // downwards transfer
  if (pband->pmy_rad->downwards_flag) {
    for (int i = ie; i>=is; --i) {
      tau_cell = pband->tau_cell(n,k,j,i);
      dx = pmb->pcoord->dx1f(i);
      emission_coefficient = pband->emission_coefficient(n,k,j,i);

      Fin = pband->flux_density_down(n,k,j,i+1);
      OneCellRadiativeTransfer(Fin, Fout, tau_cell, emission_coefficient, tau_cell/dx);
      pband->flux_density_down(n,k,j,i) = Fout;
    }
  }
  // upwards transfer
  if (pband->pmy_rad->upwards_flag) {
    for (int i = is; i<=ie; ++i) {
      tau_cell = pband->tau_cell(n,k,j,i);
      dx = pmb->pcoord->dx1f(i);
      emission_coefficient = pband->emission_coefficient(n,k,j,i);

      Fin = pband->flux_density_up(n,k,j,i);
      OneCellRadiativeTransfer(Fin, Fout, tau_cell, emission_coefficient, tau_cell/dx);
      pband->flux_density_up(n,k,j,i+1) = Fout;
    }
  }

  // Assign Energy to Absorbers
  Absorber *pabs;
  Real delta_F, fraction_absorbed, J_i;
  for (int a = 0; a < pband->nabs; ++a) {
    pabs = pband->absorbers(a);

    for (int i = is; i <= ie; ++i)
    {
      dx = pmb->pcoord->dx1f(i);
      Fin  = pband->flux_density_up(n,k,j,i) + pband->flux_density_down(n,k,j,i+1);
      Fout = pband->flux_density_up(n,k,j,i+1) + pband->flux_density_down(n,k,j,i);

      //energy absorbed by this absorber
      // J_i = (Fin-Fout)/dx * (a_i * dx / tau)
      // J_i = (Fin-Fout) * a_i/tau
      tau_cell = pband->tau_cell(n,k,j,i);
      fraction_absorbed = pabs->absorptionCoefficient(n,k,j,i)/tau_cell;
      pabs->energyAbsorbed(n,k,j,i) = (Fin - Fout) * fraction_absorbed;
    }
  }
}

void SimpleRTSolver::OneCellRadiativeTransfer(Real Fin, Real &Fout, Real tau, Real j, Real alpha) 
{
  Real attenuation = exp(-tau);

  // Source Function
  // emission coefficient / absorption coefficient
  Real s = j / alpha;

  // No Emission Case
  Fout = Fin * attenuation;

  // Emission
  // s/2 for isotropic 1D emission
  Fout += s/2. * (1. - attenuation);
}