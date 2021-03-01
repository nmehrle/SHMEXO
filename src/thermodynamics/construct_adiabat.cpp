#include <iostream>
#include <cstdlib>

#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "thermodynamics.hpp"
#include "thermodynamic_funcs.hpp"

void Thermodynamics::ConstructAdiabat(Real **w, Real Ts, Real Ps,
  Real grav, Real dz, int len, Adiabat method, Real dTdz) const
{
  Real gamma = pmy_block_->peos->GetGamma();

  // mass to molar mixing ratio
  Real q1[NHYDRO];
  put_hat(q1, w[0], eps_, Rd_);
  q1[IDN] = Ts;
  q1[IPR] = Ps;

  Real rcp[NMASS];  // molar cp ratio
  for (int n = 0; n < NMASS;  ++n)
    rcp[n] = rcp_[n]*eps_[n];

  int isat[1+NVAPOR];
  for (int n = 0; n <= NVAPOR; ++n)
    isat[n] = 0;

  // equilibrate with primary condensate, nc = n + NVAPOR
  for (int n = 1; n <= NVAPOR; ++n) {
    int nc = n + NVAPOR;
    Real rate = GasCloudIdeal(q1, n, nc, t3_[nc], p3_[nc], 0., beta_[nc], delta_[nc]);
    q1[n] -= rate;
    q1[nc] += rate;
    if (rate > 0.) isat[n] = 1;
  }

  cut_hat(w[0], q1, eps_, Rd_);
  for (int i = 1; i < len; ++i) {
    // RK4 integration 
    rk4_integrate_z_adaptive(q1, isat, rcp, eps_, beta_, delta_, t3_, p3_, gamma,
      grav/Rd_, dz, ftol_, (int)method, dTdz);
    cut_hat(w[i], q1, eps_, Rd_);
    for (int n = IVX; n <= IVZ; ++n) w[i][n] = w[0][n];
  }
}
