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

enum {iH2O = 1, iNH3 = 2, iH2Oc = 3, NH3c = 4, iH2Op = 5, iNH3p = 6};

Real P0, T0, grav;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(7);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "theta_v");
  SetUserOutputVariableName(3, "mse");
  SetUserOutputVariableName(4, "theta_e");
  SetUserOutputVariableName(5, "rh_h2o");
  SetUserOutputVariableName(6, "rh_nh3");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  Real dq[1+NVAPOR], rh;

  for (int i = is; i <= ie; ++i) {
    user_out_var(0,i) = pthermo->Temp(phydro->w.at(i));
    user_out_var(1,i) = pthermo->Theta(phydro->w.at(i), P0);
    // theta_v
    user_out_var(2,i) = user_out_var(1,i)*pthermo->Qeps(phydro->w.at(i));
    user_out_var(3,i) = pthermo->MSE(phydro->w.at(i), grav*pcoord->x1v(i));
    user_out_var(4,i) = pthermo->ThetaE(phydro->w.at(i), P0);
    pthermo->SaturationSurplus(dq, phydro->w.at(i));
    Real rh = phydro->w(iH2O,i)/(phydro->w(iH2O,i) - dq[iH2O]);
    user_out_var(5,i) = rh;
    rh = phydro->w(iNH3,i)/(phydro->w(iNH3,i) - dq[iNH3]);
    user_out_var(6,i) = rh;
  }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  grav = - pin->GetReal("hydro", "grav_acc1");
  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real qH2O = pin->GetReal("problem", "qH2O");
  Real qNH3 = pin->GetReal("problem", "qNH3");
  Real gamma = pin->GetReal("hydro", "gamma");

  // construct a 1D pseudo-moist adiabat with given relative humidity
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;

  int nx1 = 2*pmy_mesh->mesh_size.nx1 + 1;
  Real dz = (x1max - x1min)/(nx1 - 1);
  Real **w1, *z1, *p1, *t1;
  NewCArray(w1, nx1, NHYDRO);
  z1 = new Real [nx1];
  p1 = new Real [nx1];
  t1 = new Real [nx1];

  // 1.1 estimate surface temperature and pressure
  Real Rd = pthermo->GetRd();
  Real cp = gamma/(gamma - 1.)*Rd;
  Real Ts = T0 - grav/cp*x1min;
  Real Ps = P0*pow(Ts/T0, cp/Rd);
  int max_iter = 200, iter = 0;

  w1[0][iH2O] = qH2O;
  w1[0][iNH3] = qNH3;
  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i)
    z1[i] = z1[i-1] + dz;

  Real t0, p0;
  while (iter++ < max_iter) {
    pthermo->ConstructAdiabat(w1, Ts, Ps, grav, dz, nx1, Adiabat::reversible);

    // 1.3 find TP at z = 0
    for (int i = 0; i < nx1; ++i) {
      p1[i] = w1[i][IPR];
      t1[i] = pthermo->Temp(w1[i]);
    }
    p0 = interp1(0., p1, z1, nx1);
    t0 = interp1(0., t1, z1, nx1);

    Ts += T0 - t0;
    Ps *= P0/p0;
    if ((fabs(T0 - t0) < 0.01) && (fabs(P0/p0 - 1.) < 1.E-4)) break;
  }

  if (iter > max_iter) {
    std::stringstream msg;
    msg << "### FATAL ERROR in problem generator"
        << std::endl << "maximum iteration reached."
        << std::endl << "T0 = " << t0
        << std::endl << "P0 = " << p0;
    ATHENA_ERROR(msg);
  }

  // setup initial condition
  for (int i = is; i <= ie; ++i) {
    Real buf[NHYDRO];
    interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, NHYDRO);
    buf[IVX] = buf[IVY] = buf[IVZ] = 0.;
    for (int n = 0; n < NHYDRO; ++n)
      phydro->w(n,i) = buf[n];
  }

  // setup open lower boundary
  if (pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
    for (int n = 0; n < NHYDRO; ++n)
      w1[0][n] = phydro->w(n,is);
    Real P1 = w1[0][IPR];
    Real T1 = pthermo->Temp(w1[0]);
    Real dz = pcoord->dx1f(is);
    // adiabatic extrapolate half a grid down
    pthermo->ConstructAdiabat(w1, T1, P1, grav, -dz/2., 2, Adiabat::reversible);
    for (int n = 0; n < NHYDRO; ++n)
      for (int i = 1; i <= NGHOST; ++i)
        phydro->w(n,is-i) = w1[1][n];
  }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
  FreeCArray(w1);
  delete[] z1;
  delete[] p1;
  delete[] t1;
}
