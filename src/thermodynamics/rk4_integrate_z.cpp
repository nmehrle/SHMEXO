#include <iostream>
#include <algorithm>
#include "thermodynamics.hpp"

void rk4_integrate_z(Real q[], int isat[], Real rcp[], Real const eps[],
  Real const beta[], Real const delta[], Real const t3[], Real const p3[], Real gamma,
  Real g_ov_Rd, Real dz, int method, Real user_dTdz)
{
  Real step[] = {0.5, 0.5, 1.};
  Real temp = q[IDN];
  Real pres = q[IPR];
  Real dTdz[4], chi[4];

  Real beta0[NMASS], delta0[NMASS];
  std::fill(beta0, beta0 + NMASS, 0.);
  std::fill(delta0, delta0 + NMASS, 0.);

  for (int rk = 0; rk < 4; ++rk) {
    // reset vapor
    for (int n = 1; n <= NVAPOR; ++n) {
      q[n] += q[n+NVAPOR];
      q[n+NVAPOR] = 0.;
      isat[n] = 0;
    }

    for (int n = 1; n <= NVAPOR; ++n) {
      int nc = n + NVAPOR;
      Real rate = GasCloudIdeal(q, n, nc, t3[nc], p3[nc], 0., beta[nc], delta[nc]);
      //Real rate = 0.;
      if (rate > 0.) isat[n] = 1;
      q[n] -= rate;
      if (method == 0) q[nc] += rate;
    }

    // calculate tendency
    update_gamma(gamma, rcp, q);
    if (method == 0 || method == 1)
      chi[rk] = dlnTdlnP(q, isat, rcp, beta, delta, t3, gamma);
    else if (method == 2)
      chi[rk] = dlnTdlnP(q, isat, rcp, beta0, delta0, t3, gamma);
    else  // isothermal
      chi[rk] = 0.;
    dTdz[rk] = - chi[rk]*g_ov_Rd/qhat_eps(q, eps) + user_dTdz;
    chi[rk] = - qhat_eps(q, eps)/g_ov_Rd*dTdz[rk];
    
    // integrate over dz
    Real chi_avg;
    if (rk < 3) {
      q[IDN] = temp + dTdz[rk]*dz*step[rk];
      chi_avg = chi[rk];
    } else {
      q[IDN] = temp + 1./6.*(dTdz[0] + 2.*dTdz[1] + 2.*dTdz[2] + dTdz[3])*dz;
      chi_avg = 1./6.*(chi[0] + 2.*chi[1] + 2.*chi[2] + chi[3]);
    }
    if (!(q[IDN] > 0.)) q[IDN] = temp;
    if (fabs(q[IDN] - temp) > 0.1) // isothermal limit
      q[IPR] = pres*pow(q[IDN]/temp, 1./chi_avg);
    else
      q[IPR] = pres*exp(-2.*g_ov_Rd*dz/(qhat_eps(q, eps)*(q[IDN] + temp)));
  }

  // recondensation
  for (int n = 1; n <= NVAPOR; ++n) {
    int nc = n + NVAPOR;
    Real rate = GasCloudIdeal(q, n, nc, t3[nc], p3[nc], 0., beta[nc], delta[nc]);
    //Real rate = 0.;
    if (rate > 0.) isat[n] = 1;
    q[n] -= rate;
    if (method == 0) q[nc] += rate;
  }
}

void rk4_integrate_z_adaptive(Real q[], int isat[], Real rcp[], Real const eps[],
  Real const beta[], Real const delta[], Real const t3[], Real const p3[], Real gamma,
  Real g_ov_Rd, Real dz, Real ftol, int method, Real user_dTdz)
{
  Real q1[NHYDRO], q2[NHYDRO];
  int  isat1[1+NVAPOR], isat2[1+NVAPOR];

  for (int n = 0; n < NHYDRO; ++n) {
    q1[n] = q[n];
    q2[n] = q[n];
  }
  for (int n = 0; n <= NVAPOR; ++n) {
    isat1[n] = isat[n];
    isat2[n] = isat[n];
  }

  // trail step
  rk4_integrate_z(q1, isat1, rcp, eps, beta, delta, t3, p3, gamma,
                  g_ov_Rd, dz, method, user_dTdz);

  // refined step
  rk4_integrate_z(q2, isat2, rcp, eps, beta, delta, t3, p3, gamma,
                  g_ov_Rd, dz/2., method, user_dTdz);
  rk4_integrate_z(q2, isat2, rcp, eps, beta, delta, t3, p3, gamma,
                  g_ov_Rd, dz/2., method, user_dTdz);

  // abort if dz is less than 0.001 m
  //std::cout << dz << " ";
  //std::cout << fabs(q2[IDN] - q1[IDN]) << " ";
  //std::cout << ftol << std::endl;
  //assert(dz > 0.001);

  if (fabs(q2[IDN] - q1[IDN]) > ftol) {
    rk4_integrate_z_adaptive(q, isat, rcp, eps, beta, delta, t3, p3, gamma,
                             g_ov_Rd, dz/2., ftol, method, user_dTdz);
    rk4_integrate_z_adaptive(q, isat, rcp, eps, beta, delta, t3, p3, gamma,
                             g_ov_Rd, dz/2., ftol, method, user_dTdz);
  } else {
    for (int n = 0; n < NHYDRO; ++n)
      q[n] = q2[n];
    for (int n = 0; n <= NVAPOR; ++n)
      isat[n] = isat2[n];
  }
}
