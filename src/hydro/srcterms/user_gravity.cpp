//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief source terms due to custom gravity 

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../hydro.hpp"
#include "hydro_srcterms.hpp"

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::MyGravFunc
//  \brief Adds source terms for my grav

void HydroSourceTerms::AddUserGravityTerm(const Real dt,const AthenaArray<Real> *flux,
                                            const AthenaArray<Real> &prim,
                                            AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;

  // acceleration in 1-direction
  if (GravityDefined(1)) {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real src = dt*prim(IDN,k,j,i)*g1(k,j,i);
          cons(IM1,k,j,i) += src;
          if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += src*prim(IVX,k,j,i);
        }
      }
    }
  }

  // acceleration in 2-direction
  if (GravityDefined(2)) {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real src = dt*prim(IDN,k,j,i)*g2(k,j,i);
          cons(IM2,k,j,i) += src;
          if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += src*prim(IVY,k,j,i);
        }
      }
    }
  }

  // acceleration in 3-direction
  if (GravityDefined(3)) {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real src = dt*prim(IDN,k,j,i)*g3(k,j,i);
          cons(IM3,k,j,i) += src;
          if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += src*prim(IVZ,k,j,i);
        }
      }
    }
  }

  return;
}
