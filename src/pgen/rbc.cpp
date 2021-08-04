// C/C++ header
#include <ctime>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../globals.hpp"
#include "../math/interpolation.h"
#include "../utils/utils.hpp"
#include "../thermodynamics/thermodynamics.hpp"

// molecules
Real grav, P0, T0, hrate;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  Real dq[1+NVAPOR], rh;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0,k,j,i) = pthermo->Temp(phydro->w.at(k,j,i));
        user_out_var(1,k,j,i) = pthermo->Theta(phydro->w.at(k,j,i), P0);
      }
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
    AthenaArray<Real> const &w, AthenaArray<Real> const &bcc, AthenaArray<Real> &u)
{
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  Thermodynamics *pthermo = pmb->pthermo;

  if (pmb->pbval->block_bcs[outer_x1] != BoundaryFlag::block)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        Real cv = pthermo->Cv(w.at(k,j,ie));
        u(IEN,k,j,ie) += dt*hrate*w(IDN,k,j,ie)*cv;
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  grav = - pin->GetReal("hydro", "grav_acc1");
  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");
  hrate = pin->GetReal("problem", "hrate")/86400.;  // K/day to K/s

  EnrollUserExplicitSourceFunction(Forcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
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

  // estimate surface temperature and pressure
  Real Rd = pthermo->GetRd();
  Real cp = gamma/(gamma - 1.)*Rd;
  Real Ts = T0 - grav/cp*x1min;
  Real Ps = P0*pow(Ts/T0, cp/Rd);
  int max_iter = 20, iter = 0;

  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i)
    z1[i] = z1[i-1] + dz;

  Real t0, p0;
  while (iter++ < max_iter) {
    pthermo->ConstructAdiabat(w1, Ts, Ps, grav, dz, nx1, Adiabat::pseudo);

    // find TP at z = 0
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
  int kl = block_size.nx3 == 1 ? ks : ks-1;
  int ku = block_size.nx3 == 1 ? ke : ke+1;
  int jl = block_size.nx2 == 1 ? js : js-1;
  int ju = block_size.nx2 == 1 ? je : je+1;
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

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);

  FreeCArray(w1);
  delete[] z1;
  delete[] p1;
  delete[] t1;
}
