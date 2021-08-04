//! \file  noflux.cpp
//  \brief disable hydrodynamic flux

// Athena++ headers
#include "../../hydro.hpp"
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../eos/eos.hpp"

void Hydro::RiemannSolver(int const k, int const j, int const il, int const iu,
                          int const ivx,
                          AthenaArray<Real> &wl, AthenaArray<Real> &wr,
                          AthenaArray<Real> &flx, const AthenaArray<Real> &dxw)
{
  for (int n = 0; n < NHYDRO; ++n)
    for (int i = il; i <= iu; ++i)
      flx(n,k,j,i) = 0.;
}
