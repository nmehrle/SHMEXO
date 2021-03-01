//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file straka.cpp
//  \brief Problem generator for a dense sinking bubble
//
// REFERENCE: Straka et al., 1993
//========================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
#include "../globals.hpp"
#include "../thermodynamics/thermodynamics.hpp"

namespace math {
  #include "../math/core.h"
}

Real K, p0, cp, Rd;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  for (int j = js; j <= je; ++j)
    for (int i = is; i <= ie; ++i) {
      Real prim[NHYDRO];
      for (int n = 0; n < NHYDRO; ++n)
        prim[n] = phydro->w(n,j,i);
      user_out_var(0,j,i) = pthermo->Temp(phydro->w.at(j,i));
      user_out_var(1,j,i) = pthermo->Theta(phydro->w.at(j,i), p0);
    }
}

void Diffusion(MeshBlock *pmb, Real const time, Real const dt,
  AthenaArray<Real> const& w, AthenaArray<Real> const& bcc, AthenaArray<Real> &u)
{
  Real dx = pmb->pcoord->dx1f(0);
  Thermodynamics *pthermo = pmb->pthermo;
  for (int j = pmb->js; j <= pmb->je; ++j)
    for (int i = pmb->is; i <= pmb->ie; ++i) {
      Real prim[NHYDRO];
      Real temp = pthermo->Temp(w.at(j,i));
      Real theta = pthermo->Theta(w.at(j,i), p0);
      Real theta_ip1_j = pthermo->Theta(w.at(j+1,i), p0);
      Real theta_im1_j = pthermo->Theta(w.at(j-1,i), p0);
      Real theta_i_jp1 = pthermo->Theta(w.at(j,i+1), p0);
      Real theta_i_jm1 = pthermo->Theta(w.at(j,i-1), p0);

      u(IM1,j,i) += dt*K*w(IDN,j,i)/(dx*dx)*(
        w(IV1,j,i-1) + w(IV1,j,i+1) + w(IV1,j-1,i) + w(IV1,j+1,i) - 4.*w(IV1,j,i));
      u(IM2,j,i) += dt*K*w(IDN,j,i)/(dx*dx)*(
        w(IV2,j,i-1) + w(IV2,j,i+1) + w(IV2,j-1,i) + w(IV2,j+1,i) - 4.*w(IV2,j,i));
      u(IEN,j,i) += dt*K*w(IDN,j,i)/(dx*dx)*cp*temp/theta*(theta_ip1_j + theta_im1_j +
        theta_i_jp1 + theta_i_jm1 - 4.*theta);
    }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // forcing parameters
  Real gamma = pin->GetReal("hydro", "gamma");
  K  = pin->GetReal("problem", "K");
  p0 = pin->GetReal("problem", "p0");
  Rd = pin->GetReal("thermodynamics", "Rd");
  cp = gamma/(gamma - 1.)*Rd;

  // forcing function
  EnrollUserExplicitSourceFunction(Diffusion);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real grav = -phydro->hsrc.GetG1();
  Real Ts = pin->GetReal("problem", "Ts");
  Real xc = pin->GetReal("problem", "xc");
  Real xr = pin->GetReal("problem", "xr");
  Real zc = pin->GetReal("problem", "zc");
  Real zr = pin->GetReal("problem", "zr");
  Real dT = pin->GetReal("problem", "dT");

  for (int j = js; j <= je; ++j) {
    for (int i = is; i <= ie; ++i) {
      Real x1 = pcoord->x1v(i);
      Real x2 = pcoord->x2v(j);
      Real L = sqrt(math::sqr((x2 - xc)/xr) + math::sqr((x1 - zc)/zr));
      Real temp = Ts - grav*x1/cp;
      phydro->w(IPR,j,i) = p0*pow(temp/Ts, cp/Rd);
      if (L <= 1.)
        temp += dT*(cos(M_PI*L) + 1.)/2.;
      phydro->w(IDN,j,i) = phydro->w(IPR,j,i)/(Rd*temp);
      phydro->w(IV1,j,i) = phydro->w(IV2,j,i) = 0.;
    }
  }

  if (pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
    for (int j = js-1; j <= je+1; ++j)
      for (int i = 1; i <= NGHOST; ++i) {
        Real temp = Ts - grav*pcoord->x1f(is)/cp;
        phydro->w(IPR,j,is-i) = p0*pow(temp/Ts, cp/Rd);
        phydro->w(IDN,j,is-i) = phydro->w(IPR,j,is-i)/(Rd*temp);
        phydro->w(IVX,j,is-i) = phydro->w(IVY,j,is-i) = 0.;
      }
  }

  /* add tracers if defined
  if (NTRACER > 0) {
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        phydro->w(ITR,j,i) = pcoord->x1v(i);
        if (j > 7./8.*js + 1./8.*je)
          phydro->w(ITR,j,i) += 1000.;
        if (NTRACER > 1) {
          phydro->w(ITR+1,j,i) = pcoord->x2v(j);
          if (i > 7./8.*is + 1./8.*ie)
            phydro->w(ITR+1,j,i) += 1000.;
        }
      }
  }*/

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
