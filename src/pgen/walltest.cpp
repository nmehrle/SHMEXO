
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file blast.cpp
//  \brief Problem generator for spherical blast wave problem.  Works in Cartesian,
//         cylindrical, and spherical coordinates.  Contains post-processing code
//         to check whether blast is spherical for regression tests
//
// REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
//   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

Real threshold;

Real channel_xmax, surface_ymin;
Real channel_x, surface_y;

Real pa, da, prat, drat, gamma_, y1_, y2_;

void Mesh::InitUserMeshData(ParameterInput *pin) {
  channel_xmax = pin->GetReal("problem", "channel_xmax");
  surface_ymin = pin->GetReal("problem", "surface_ymin");

  pa   = pin->GetOrAddReal("problem", "pamb", 1.0);
  da   = pin->GetOrAddReal("problem", "damb", 1.0);
  prat = pin->GetReal("problem", "prat");
  drat = pin->GetOrAddReal("problem", "drat", 1.0);

  y1_ = pin->GetOrAddReal("problem", "y1", -1.0);
  y2_ = pin->GetOrAddReal("problem", "y2", -0.8);
  
  gamma_ = pin->GetReal("hydro", "gamma");

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Spherical blast wave test problem generator
//========================================================================================

void Mesh::UserAdjustBlockBCs(LogicalLocation loc, RegionSize &block_size,
                                  BoundaryFlag *block_bcs) {
  // channel
  // if block is lower than surface level
  // and contains channel xmax in it
  // set right bc to reflecting
  if (block_size.x2max <= surface_ymin) {
    if (block_size.x1min < channel_xmax && block_size.x1max >= channel_xmax) {
      block_bcs[BoundaryFace::outer_x1] = BoundaryFlag::reflect;
      channel_x = block_size.x1max;
    }

    if (block_size.x1min <= channel_x) {
      block_bcs[BoundaryFace::inner_x1] = BoundaryFlag::reflect;
    }
  }

  // floor
  // if block is rightward of channel
  // and contains floor level
  // set bottom bc to reflect
  if (block_size.x1min >= channel_xmax) {
    if (block_size.x2min <= surface_ymin && block_size.x2max > surface_ymin) {
      block_bcs[BoundaryFace::inner_x2] = BoundaryFlag::reflect;
      surface_y = block_size.x2min;
    }

    if (block_size.x2max < surface_y) {
      block_bcs[BoundaryFace::outer_x2] = BoundaryFlag::reflect;
    }
  }
 }

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real gamma = gamma_;
  Real gm1 = gamma - 1.0;
  Real y1 = y1_;
  Real y2 = y2_;

  std::cout << "surface y:: " << surface_y << std::endl;
  std::cout << "channel_x:: " << channel_x << std::endl;

  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real rad;

        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);

        Real den = da;
        Real pres = pa;

        if (x < channel_x) {
          if(y < y1) {
            den = drat *da;
            pres = prat * da;
          }
          else if (y > y1 && y < y2) {
            Real f = (y-y1) / (y2 - y1);
            Real log_den = (1.0-f) * std::log(da*drat) + f * std::log(da);
            den = std::exp(log_den);

            Real log_pres = (1.0-f) * std::log(pa*prat) + f * std::log(pa);
            pres = std::exp(log_pres);
          }
        }

        phydro->u(IDN,k,j,i) = den;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = pres/gm1;

      } // i
    } // j
  } // k
}

void MeshBlock::UserWorkInLoop() {
  Real gamma = gamma_;
  Real y1 = y1_;
  Real y2 = y2_;
  Real gm1 = gamma - 1.0;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real rad;

        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);

        Real den = da;
        Real pres = pa;

        if (x < channel_x) {
          if(y < y1) {
            den = drat *da;
            pres = prat * da;

            phydro->u(IDN,k,j,i) = den;
            phydro->u(IM1,k,j,i) = 0.0;
            phydro->u(IM2,k,j,i) = 0.0;
            phydro->u(IM3,k,j,i) = 0.0;
            phydro->u(IEN,k,j,i) = pres/gm1;
          }
        }

      } // i
    } // j
  } // k
}