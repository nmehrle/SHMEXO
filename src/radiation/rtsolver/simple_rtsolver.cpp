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

  Real tau, dx, emission;
  // Fin,Fout -- flux in/out due to input flux
  // Sin,Sout -- flux in/out due to source term
  Real Fin, Fout, Sin, Sout;
  // downwards transfer
  if (pband->pmy_rad->downwards_flag) {
    for (int i = ie; i>=is; --i) {
      dx = pmb->pcoord->dx1f(i);
      tau = pband->tau_cell(n,k,j,i);
      emission = pband->emission_coefficient(n,k,j,i);

      Fin = pband->flux_down(n,k,j,i+1);
      OneCellRadiativeTransfer(Fin, Fout, tau, 0, dx);
      pband->flux_down(n,k,j,i) = Fout;

      Sin = pband->source_flux_down(n,k,j,i+1);
      OneCellRadiativeTransfer(Sin, Sout, tau, emission, dx);
      pband->source_flux_down(n,k,j,i) = Sout;
    }
  }

  // Assign Energy to Absorbers
  Absorber *pabs;
  Real netFup   = 0;
  Real netFdown = 0;
  Real netSup   = 0;
  Real netSdown = 0;
  Real netE_F, netE_S, Eabs;
  Real absorber_fraction;

  for (int i = is; i <= ie; ++i) {
    dx = pmb->pcoord->dx1f(i);
    tau = pband->tau_cell(n,k,j,i);
    emission = pband->emission_coefficient(n,k,j,i);

    // Upwards Transfer
    if (pband->pmy_rad->upwards_flag) {
      Fin = pband->flux_up(n,k,j,i);
      OneCellRadiativeTransfer(Fin, Fout, tau, 0, dx);
      pband->flux_up(n,k,j,i+1) = Fout;

      Sin = pband->source_flux_up(n,k,j,i);
      OneCellRadiativeTransfer(Sin, Sout, tau, emission, dx);
      pband->source_flux_up(n,k,j,i+1) = Sout;

      netFup = Fin-Fout;
      netSup = Sin-Sout;
    }

    if (pband->pmy_rad->downwards_flag) {
      netFdown = pband->flux_down(n,k,j,i+1) - pband->flux_down(n,k,j,i);
      netSdown = pband->source_flux_down(n,k,j,i+1) - pband->source_flux_down(n,k,j,i);
    }

    // for debug
    Real FinDown  = pband->flux_down(n,k,j,i+1);
    Real FoutDown = pband->flux_down(n,k,j,i);
    Real SinDown  = pband->source_flux_down(n,k,j,i+1);
    Real SoutDown = pband->source_flux_down(n,k,j,i);

    netE_F = (netFup + netFdown)/dx;
    netE_S = (netSup + netSdown)/dx + emission;
    Eabs = netE_F + netE_S;

    for (int a = 0; a < pband->nabs; ++a) {
      pabs = pband->absorbers(a);
      absorber_fraction = (pabs->absorptionCoefficient(n,k,j,i) * dx)/tau;
      pabs->energyAbsorbed(n,k,j,i) = Eabs * absorber_fraction;
    }
  }
}

void SimpleRTSolver::OneCellRadiativeTransfer(Real Fin, Real &Fout, Real tau, Real emission, Real dx) 
{
  // No Emission Case
  Fout = Fin * exp(-tau);

  // Emission
  // j/2 for isotropic 1D emission
  Fout += emission/2. * dx * exp(-tau/2.);
}