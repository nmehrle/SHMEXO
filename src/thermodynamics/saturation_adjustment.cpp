// C/C++ header
#include <cstring>
#include <cstdlib>
#include <sstream>

// Athena++ header
#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "thermodynamics.hpp"
#include "thermodynamic_funcs.hpp"

void Thermodynamics::SaturationAdjustment(AthenaArray<Real> &u) const
{
  // return if there's no vapor
  if (NVAPOR*NPHASE == 0) return;

  MeshBlock *pmb = pmy_block_;
  Real gamma = pmb->peos->GetGamma();

  Real q[NHYDRO], q0[NHYDRO];
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        // mass to molar mixing ratio
        Real rho = 0., rho_hat = 0.;
        for (int n = 0; n < NMASS; ++n) {
          rho += u(n,k,j,i);
          rho_hat += u(n,k,j,i)/eps_[n];
          q[n] = u(n,k,j,i)/eps_[n];
        }
        for (int n = 0; n < NMASS; ++n)
          q[n] /= rho_hat;

        // calculate internal energy
        Real KE = 0.5*(u(IM1,k,j,i)*u(IM1,k,j,i) 
                     + u(IM2,k,j,i)*u(IM2,k,j,i)
                     + u(IM3,k,j,i)*u(IM3,k,j,i))/rho;
        Real uhat = (u(IEN,k,j,i) - KE)/(Rd_*rho_hat);
        UpdateTPConservingU(q, rho, uhat);

        // save q for debug purpose
        memcpy(q0, q, NHYDRO*sizeof(Real));

        // check boiling condition
        for (int n = 1; n <= NVAPOR; ++n) {
          int ic = n + NVAPOR, ip = n + 2*NVAPOR;
          Real svp = SatVaporPressure(q[IDN], ic);
          if (svp > q[IPR]) { // boiling
            q[ic] += q[ip];
            q[ip] = 0.;
          }
        }

        Real qsig = 1.;
        for (int n = 1; n < NMASS; ++n)
          qsig += q[n]*(rcv_[n]*eps_[n] - 1.);

        int iter = 0;
        Real t, alpha, rate;
        while (iter++ < max_iter_) {
          // condensation
          for (int n = 1; n <= NVAPOR; ++n) {
            int nc = n + NVAPOR;
            t = q[IDN]/t3_[nc];
            alpha = (gamma - 1.)*(beta_[nc]/t - delta_[nc] - 1.)/qsig;
            rate = GasCloudIdeal(q, n, nc, t3_[nc], p3_[nc], alpha, beta_[nc], delta_[nc]);
            q[n] -= rate;
            q[nc] += rate;
          }
          Real Told = q[IDN];
          UpdateTPConservingU(q, rho, uhat);
          if (fabs(q[IDN] - Told) < ftol_) break;
        }
        if (iter >= max_iter_) {
          std::stringstream msg;
          msg << "### FATAL ERROR in function Thermodynamics::SaturationAdjustment"
              << std::endl << "Iteration reaches maximum."
              << std::endl;
          msg << "Location: x3 = " << pmb->pcoord->x3v(k) << std::endl
              << "          x2 = " << pmb->pcoord->x2v(j) << std::endl
              << "          x1 = " << pmb->pcoord->x1v(i) << std::endl;
          msg << "Variables before iteration u = (";
          for (int n = 0; n < NHYDRO-1; ++n)
            msg << u(n,k,j,i) << ", ";
          msg << u(NHYDRO-1,k,j,i) << ")" << std::endl;
          msg << "Variables before iteration q0 = (";
          for (int n = 0; n < NHYDRO-1; ++n)
            msg << q0[n] << ", ";
          msg << q0[NHYDRO-1] << ")" << std::endl;
          msg << "Variables after iteration q = (";
          for (int n = 0; n < NHYDRO-1; ++n)
            msg << q[n] << ", ";
          msg << q[NHYDRO-1] << ")" << std::endl;
          ATHENA_ERROR(msg);
        }

        // from molar to mass mixing ratios
        Real sum = 1.;
        for (int n = 1; n < NMASS; ++n)
          sum += q[n]*(eps_[n] - 1.);
        for (int n = 1; n < NMASS; ++n)
          u(n,k,j,i) = rho*q[n]*eps_[n]/sum;
      }
}
