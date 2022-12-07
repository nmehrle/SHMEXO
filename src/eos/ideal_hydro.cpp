//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file adiabatic_hydro.cpp
//  \brief implements functions in class EquationOfState for adiabatic hydrodynamics`

// C headers

// C++ headers
#include <cmath>   // sqrt()
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../globals.hpp"
#include "eos.hpp"

// EquationOfState constructor

EquationOfState::EquationOfState(MeshBlock *pmb, ParameterInput *pin) :
    pmy_block_(pmb),
    gamma_{pin->GetReal("hydro", "gamma")},
    density_floor_{pin->GetOrAddReal("hydro", "dfloor", std::sqrt(1024*float_min))},
    pressure_floor_{pin->GetOrAddReal("hydro", "pfloor", std::sqrt(1024*float_min))},
    scalar_floor_{pin->GetOrAddReal("hydro", "sfloor", std::sqrt(1024*float_min))},
    boltzmann_{pin->GetOrAddReal("hydro", "boltzmann", 0)} {}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::ConservedToPrimitive(AthenaArray<Real> &cons,
//           const AthenaArray<Real> &prim_old, const FaceField &b,
//           AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Coordinates *pco,
//           int il, int iu, int jl, int ju, int kl, int ku)
// \brief Converts conserved into primitive variables in adiabatic hydro.

void EquationOfState::ConservedToPrimitive(
    AthenaArray<Real> &cons, const AthenaArray<Real> &prim_old, const FaceField &b,
    AthenaArray<Real> &prim, AthenaArray<Real> &bcc,
    Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku) {
  Real gm1 = GetGamma() - 1.0;
  std::stringstream msg;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real& u_d  = cons(IDN,k,j,i);
        Real& u_m1 = cons(IM1,k,j,i);
        Real& u_m2 = cons(IM2,k,j,i);
        Real& u_m3 = cons(IM3,k,j,i);
        Real& u_e  = cons(IEN,k,j,i);

        Real& w_d  = prim(IDN,k,j,i);
        Real& w_vx = prim(IVX,k,j,i);
        Real& w_vy = prim(IVY,k,j,i);
        Real& w_vz = prim(IVZ,k,j,i);
        Real& w_p  = prim(IPR,k,j,i);

        // apply density floor, without changing momentum or energy
        u_d = (u_d > density_floor_) ?  u_d : density_floor_;
        w_d = u_d;

        w_vx = u_m1/u_d;
        w_vy = u_m2/u_d;
        w_vz = u_m3/u_d;

        Real e_k = (0.5/u_d)*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3));
        w_p = gm1*(u_e - e_k);

        // apply pressure floor, correct total energy
        u_e = (w_p > pressure_floor_) ?  u_e : ((pressure_floor_/gm1) + e_k);
        w_p = (w_p > pressure_floor_) ?  w_p : pressure_floor_;
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::PrimitiveToConserved(const AthenaArray<Real> &prim,
//           const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
//           int il, int iu, int jl, int ju, int kl, int ku);
// \brief Converts primitive variables into conservative variables

void EquationOfState::PrimitiveToConserved(
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bc,
    AthenaArray<Real> &cons, Coordinates *pco,
    int il, int iu, int jl, int ju, int kl, int ku) {
  Real igm1 = 1.0/(GetGamma() - 1.0);

  // Force outer-loop vectorization
#pragma omp simd
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      //#pragma omp simd
#pragma novector
      for (int i=il; i<=iu; ++i) {
        Real& u_d  = cons(IDN,k,j,i);
        Real& u_m1 = cons(IM1,k,j,i);
        Real& u_m2 = cons(IM2,k,j,i);
        Real& u_m3 = cons(IM3,k,j,i);
        Real& u_e  = cons(IEN,k,j,i);

        const Real& w_d  = prim(IDN,k,j,i);
        const Real& w_vx = prim(IVX,k,j,i);
        const Real& w_vy = prim(IVY,k,j,i);
        const Real& w_vz = prim(IVZ,k,j,i);
        const Real& w_p  = prim(IPR,k,j,i);

        // density
        u_d = w_d;

        // momentum
        u_m1 = w_vx*w_d;
        u_m2 = w_vy*w_d;
        u_m3 = w_vz*w_d;

        // total energy
        Real KE = 0.5*w_d*(w_vx*w_vx + w_vy*w_vy + w_vz*w_vz);
        u_e = igm1*w_p+ KE;
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
// \!fn Real EquationOfState::SoundSpeed(Real prim[NHYDRO])
// \brief returns adiabatic sound speed given vector of primitive variables
Real EquationOfState::SoundSpeed(const Real prim[NHYDRO]) {
  return std::sqrt(gamma_*prim[IPR]/prim[IDN]);
}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim, int k, int j,
//                                                 =int i)
// \brief Apply density and pressure floors to reconstructed L/R cell interface states
void EquationOfState::ApplyPrimitiveFloors(AthenaArray<Real> &prim, int k, int j, int i) {
  Real& w_d  = prim(IDN,i);
  Real& w_p  = prim(IPR,i);

  // apply (prim) density floor
  w_d = (w_d > density_floor_) ?  w_d : density_floor_;
  // apply pressure floor
  w_p = (w_p > pressure_floor_) ?  w_p : pressure_floor_;

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::ApplyPrimitiveConservedFloors(AthenaArray<Real> &prim,
//!           AthenaArray<Real> &cons, FaceField &b, int k, int j, int i) {
//! \brief Apply pressure (prim) floor and correct energy (cons) (typically after W(U))
void EquationOfState::ApplyPrimitiveConservedFloors(
    AthenaArray<Real> &prim, AthenaArray<Real> &cons, AthenaArray<Real> &bcc,
    int k, int j, int i) {
  Real gm1 = GetGamma() - 1.0;
  Real& w_d  = prim(IDN,k,j,i);
  Real& w_p  = prim(IPR,k,j,i);

  Real& u_d  = cons(IDN,k,j,i);
  Real& u_e  = cons(IEN,k,j,i);
  // apply (prim) density floor, without changing momentum or energy
  w_d = (w_d > density_floor_) ?  w_d : density_floor_;
  // ensure cons density matches
  u_d = w_d;

  Real e_k = 0.5*w_d*(SQR(prim(IVX,k,j,i)) + SQR(prim(IVY,k,j,i))
                      + SQR(prim(IVZ,k,j,i)));
  // apply pressure floor, correct total energy
  u_e = (w_p > pressure_floor_) ?
        u_e : ((pressure_floor_/gm1) + e_k);
  w_p = (w_p > pressure_floor_) ?
        w_p : pressure_floor_;

  return;
}


void EquationOfState::Temperature(const AthenaArray<Real> &w, const AthenaArray<Real> &s, const AthenaArray<Real> &m,
  Real &T, int k, int j, int i)
{
  Real kb = boltzmann_;
  if ( kb <= 0 ) {
    std::stringstream msg;
    msg << "### FATAL ERROR in EOS::Temperature" << std::endl
        << "    Boltzmann constant (<hydro> boltzmann) not defined in problem file" << std::endl;
        << "    Or defined as negative." << std::endl;
    ATHENA_ERROR(msg);   
  }

  // Sum number densities of all species
  Real number_density = 0;
  Real dens, mass;
  for (int n = 0; n < NSCALARS; ++n)
  {
    dens = s(n, k, j, i);
    mass = m(n);
    number_density += dens/mass;
  }

  // P = n K T
  T = w(IPR, k, j, i) / (number_density * kb);
}