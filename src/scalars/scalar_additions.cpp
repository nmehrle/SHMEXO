//========================================================================================
// Athena++ astrophysical MHD code
//========================================================================================
//! \file scalar_additions.cpp
//  \brief additions to scalar class

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "scalars.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn  void PassiveScalars::AddFluxDivergence
//  \brief Adds flux divergence to weighted average of conservative variables from

void PassiveScalars::ImplicitCorrectionFull(AthenaArray<Real> &ds, AthenaArray<Real> const&implicit_correction) {
  // implicit_correction -- NHYDRO, ncells3, ncells2, ncells1
  // ds -- NSCALARS, ncells3, ncells2, ncells1

  MeshBlock *pmb = pmy_block;

  int ks = pmb->ks; int js = pmb->js; int is = pmb->is;
  int ke = pmb->ke; int je = pmb->je; int ie = pmb->ie;

  AthenaArray<Real> scalar_correction;
  scalar_correction.NewAthenaArray(NSCALARS, pmb->ncells3, pmb->ncells2, pmb->ncells1);

  for (int n=0; n<NSCALARS; ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          scalar_correction(n,k,j,i) = r(n,k,j,i) * implicit_correction(IDN, k, j, i);
        }
      }
    }
  }

  Real wghts[3];
  wghts[0] = 1.;
  wghts[1] = 1.;
  wghts[2] = 0.;
  pmb->WeightedAve(ds, scalar_correction, s1, wghts);
}

void PassiveScalars::ImplicitCorrectionReduced(AthenaArray<Real> &ds, AthenaArray<Real> const&implicit_correction) {
  ImplicitCorrectionFull(ds, implicit_correction);
}