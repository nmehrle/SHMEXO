//! \file sedimentation.cpp
//  \brief Problem generator for testing sedimentation

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
#include "../thermodynamics/thermodynamics.hpp"
#include "../utils/utils.hpp"
#include "../math/core.h"
#include "../math/interpolation.h"

enum {iH2O = 1, iH2Oc = 2, iH2Op = 3};

Real p0, grav;

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
  for (int i = is; i <= ie; ++i) {
    user_out_var(0,i) = pthermo->Temp(phydro->w.at(i));
    user_out_var(1,i) = pthermo->Theta(phydro->w.at(i), p0);
    // theta_v
    user_out_var(2,i) = user_out_var(1,i)*pthermo->Qeps(phydro->w.at(i));
    user_out_var(3,i) = pthermo->MSE(phydro->w.at(i), grav*pcoord->x1v(i));
    user_out_var(4,i) = pthermo->ThetaE(phydro->w.at(i), p0);
    pthermo->SaturationSurplus(dq, phydro->w.at(i));
    Real rh = phydro->w(iH2O,i)/(phydro->w(iH2O,i) - dq[iH2O]);
    user_out_var(5,i) = rh;
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
  Real Qs = pin->GetReal("problem", "Qs");
  Real drho = pin->GetReal("problem", "drho");
  Real zmin = pin->GetReal("problem", "zmin");
  Real zmax = pin->GetReal("problem", "zmax");

  // construct a 1D pseudo-moist adiabat with given relative humidity
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;
  int nx1 = 2*pmy_mesh->mesh_size.nx1 + 1;
  Real dz = (x1max - x1min)/(nx1 - 1);
  Real **w1, *z1;
  NewCArray(w1, nx1, NHYDRO);
  z1 = new Real [nx1];

  w1[0][iH2O] = Qs;
  Real Ps = p0;
  pthermo->ConstructAdiabat(w1, Ts, Ps, grav, dz, nx1, Adiabat::reversible);
  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i)
    z1[i] = z1[i-1] + dz;

  // setup initial condition
  for (int i = is; i <= ie; ++i) {
    Real buf[NHYDRO];
    interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, NHYDRO);
    for (int n = 0; n < NHYDRO; ++n)
        phydro->w(n,i) = buf[n];
  }

  // setup open lower boundary
  if (pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
    for (int n = 0; n < NHYDRO; ++n)
      w1[0][n] = phydro->w(n,is);
    pthermo->ConstructAdiabat(w1, pthermo->Temp(w1[0]), w1[0][IPR],
      grav, -pcoord->dx1f(is), 1+NGHOST, Adiabat::reversible);
    for (int n = 0; n < NHYDRO; ++n)
      for (int i = 1; i <= NGHOST; ++i)
        phydro->w(n,is-i) = w1[i][n];
  }

  // add additional precipitation 
  for (int i = is; i <= ie; ++i) {
    Real c1 = pcoord->x1f(i), c2 = pcoord->x1f(i+1);
    Real dc = c2 - c1;
    if ((c2 > zmin) && (c1 < zmax)) {
      Real f1 = (c2 - zmin)/dc;
      Real f2 = (zmax - c1)/dc;
      phydro->w(iH2Op,i) = std::min(phydro->w(iH2Oc,i), min(1.,f1,f2)*drho/phydro->w(IDN,i));
      phydro->w(iH2Oc,i) -= phydro->w(iH2Op,i);
    }
  }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
  FreeCArray(w1);
  delete[] z1;
}
