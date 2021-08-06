//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file robert.cpp
//  \brief Problem generator for rising bubble
//
// REFERENCE: Robert et al., 1992
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

Real p0;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  Real gamma = peos->GetGamma();
  for (int j = js; j <= je; ++j)
    for (int i = is; i <= ie; ++i) {
      user_out_var(0,j,i) = pthermo->Temp(phydro->w.at(j,i));
      user_out_var(1,j,i) = pthermo->Theta(phydro->w.at(j,i), p0);
    }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  p0 = pin->GetReal("problem", "p0");
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real gamma = pin->GetReal("hydro", "gamma");
  Real grav = - pin->GetReal("hydro", "grav_acc1");
  Real Ts = pin->GetReal("problem", "Ts");
  Real Rd = pin->GetReal("thermodynamics", "Rd");
  Real cp = gamma/(gamma - 1.)*Rd;

  Real xc = pin->GetReal("problem", "xc");
  Real zc = pin->GetReal("problem", "zc");
  Real s = pin->GetReal("problem", "s");
  Real a = pin->GetReal("problem", "a");
  Real dT = pin->GetReal("problem", "dT");

  bool uniform_bubble = pin->GetOrAddBoolean("problem", "uniform_bubble", false);

  for (int j = js; j <= je; ++j) {
    for (int i = is; i <= ie; ++i) {
      Real x1 = pcoord->x1v(i);
      Real x2 = pcoord->x2v(j);
      Real r = sqrt((x2-xc)*(x2-xc) + (x1-zc)*(x1-zc));
      Real temp = Ts - grav*x1/cp;
      phydro->w(IPR,j,i) = p0*pow(temp/Ts, cp/Rd);
      if (r <= a)
        temp += dT*pow(phydro->w(IPR,j,i)/p0, Rd/cp);
      else if (!uniform_bubble)
        temp += dT*exp(-(r-a)*(r-a)/(s*s))*pow(phydro->w(IPR,j,i)/p0, Rd/cp);
      phydro->w(IDN,j,i) = phydro->w(IPR,j,i)/(Rd*temp);
      phydro->w(IVX,j,i) = phydro->w(IVY,j,i) = 0.;
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

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
