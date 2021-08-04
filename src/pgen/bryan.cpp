//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bryan.cpp
//  \brief Problem generator for a moist rising bubble
//
// REFERENCE: Bryan and Fritsch, 2002
//========================================================================================

// Athena++ headers
#include "../math/core.h"
#include "../math/root.h"
#include "../math/interpolation.h"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../utils/utils.hpp"

enum {iH2O = 1, iH2Oc = 2};

Real p0, grav;

struct RootData {
  Real temp_v;
  Real beta;
  Real delta;
  Real p3;
  Real t3;
  Real eps;
  Real *q;
};

Real root_func(Real temp, void *aux) {
  RootData *pd = static_cast<RootData*>(aux);
  Real& temp_v = pd->temp_v;
  Real& beta = pd->beta;
  Real& delta = pd->delta;
  Real& p3 = pd->p3;
  Real& t3 = pd->t3;
  Real& eps = pd->eps;
  Real *q = pd->q;

  q[IDN] = temp;
  Real rate = GasCloudIdeal(q, iH2O, iH2Oc, t3, p3, 0., beta, delta);
  q[iH2O] -= rate;
  q[iH2Oc] += rate;
  Real f1 = 1. - q[iH2Oc];
  Real f2 = 1. + (q[iH2O] + q[iH2Oc])*(eps - 1.);
  return temp*f1/f2 - temp_v;
};

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(6);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "theta_v");
  SetUserOutputVariableName(3, "mse");
  SetUserOutputVariableName(4, "theta_e");
  SetUserOutputVariableName(5, "rh");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  Real dq[1+NVAPOR];
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0,k,j,i) = pthermo->Temp(phydro->w.at(k,j,i));
        user_out_var(1,k,j,i) = pthermo->Theta(phydro->w.at(k,j,i), p0);
        // theta_v
        user_out_var(2,k,j,i) = user_out_var(1,k,j,i)*pthermo->Qeps(phydro->w.at(k,j,i));
        user_out_var(3,k,j,i) = pthermo->MSE(phydro->w.at(k,j,i), grav*pcoord->x1v(i));
        user_out_var(4,k,j,i) = pthermo->ThetaE(phydro->w.at(k,j,i), p0);
        pthermo->SaturationSurplus(dq, phydro->w.at(k,j,i));
        Real rh = phydro->w(iH2O,k,j,i)/(phydro->w(iH2O,k,j,i) - dq[iH2O]);
        user_out_var(5,k,j,i) = rh;
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  grav = - pin->GetReal("hydro", "grav_acc1");
  p0 = pin->GetReal("problem", "p0");
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real Ts = pin->GetReal("problem", "Ts");

  Real xc = pin->GetReal("problem", "xc");
  Real zc = pin->GetReal("problem", "zc");
  Real xr = pin->GetReal("problem", "xr");
  Real zr = pin->GetReal("problem", "zr");
  Real dT = pin->GetReal("problem", "dT");
  Real qt = pin->GetReal("problem", "qt");

  RootData solver;
  solver.p3 = pin->GetReal("thermodynamics", "Ptriple1");
  solver.t3 = pin->GetReal("thermodynamics", "Ttriple1");
  solver.eps = pthermo->GetMassRatio(iH2Oc);
  solver.beta = pthermo->GetBeta(iH2Oc);
  solver.delta = pthermo->GetDelta(iH2Oc);

  Real Rd = pthermo->GetRd();
  Real gamma = peos->GetGamma();
  Real cpd = gamma/(gamma - 1.)*Rd;
  
  // construct a 1D pseudo-moist adiabat with given relative humidity
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;
  int nx1 = 2*pmy_mesh->mesh_size.nx1 + 1;
  Real dz = (x1max - x1min)/(nx1 - 1);
  Real **w1, *z1;
  NewCArray(w1, nx1, NHYDRO);
  z1 = new Real [nx1];

  w1[0][iH2O] = qt;
  Real Ps = p0;
  pthermo->ConstructAdiabat(w1, Ts, Ps, grav, dz, nx1, Adiabat::reversible);
  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i)
    z1[i] = z1[i-1] + dz;

  // setup initial condition
  int kl = block_size.nx3 == 1 ? ks : ks-NGHOST;
  int ku = block_size.nx3 == 1 ? ke : ke+NGHOST;
  int jl = block_size.nx2 == 1 ? js : js-NGHOST;
  int ju = block_size.nx2 == 1 ? je : je+NGHOST;
  for (int i = is; i <= ie; ++i) {
    Real buf[NHYDRO];
    interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, NHYDRO);
    buf[IVX] = buf[IVY] = buf[IVZ] = 0.;
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = kl; k <= ku; ++k)
        for (int j = jl; j <= ju; ++j)
          phydro->w(n,k,j,i) = buf[n];
  }

  // setup open lower boundary
  if (pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j) {
        for (int n = 0; n < NHYDRO; ++n)
          w1[0][n] = phydro->w(n,k,j,is);
        Real P1 = w1[0][IPR];
        Real T1 = pthermo->Temp(w1[0]);
        Real dz = pcoord->dx1f(is);
        // adiabatic extrapolate half a grid down
        pthermo->ConstructAdiabat(w1, T1, P1, grav, -dz/2., 2, Adiabat::reversible);
        for (int n = 0; n < NHYDRO; ++n)
          for (int i = 1; i <= NGHOST; ++i)
            phydro->w(n,k,j,is-i) = w1[1][n];
      }
  }

  // add temperature anomaly
  Real q[NHYDRO]; Real temp;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real L = sqrt(sqr((x2 - xc)/xr) + sqr((x1 - zc)/zr));
        if (L < 1.) {
          pthermo->P2Q(q, phydro->w.at(k,j,i));
          solver.q = q;
          solver.temp_v = phydro->w(IPR,k,j,i)/(phydro->w(IDN,k,j,i)*Rd)*
            (dT*sqr(cos(M_PI*L/2.))/300. + 1.);
          int err = root(q[IDN], q[IDN] + dT, 1.E-8, &temp, root_func, &solver);
          if (err) {
            std::stringstream msg;
            msg << "### TVSolver doesn't converge" << std::endl;
            throw std::runtime_error(msg.str().c_str());
          } else {
            root_func(temp, &solver);
          }
          pthermo->Q2P(phydro->w.at(k,j,i), q);
        }
      }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);

  FreeCArray(w1);
  delete[] z1;
}
