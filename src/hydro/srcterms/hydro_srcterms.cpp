//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief Class to implement source terms in the hydro equations

// C headers

// C++ headers
#include <cstring>    // strcmp
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../hydro.hpp"
#include "hydro_srcterms.hpp"

// HydroSourceTerms constructor

HydroSourceTerms::HydroSourceTerms(Hydro *phyd, ParameterInput *pin) :
  pmy_hydro_(phyd),
  g1(phyd->pmy_block->ncells3, phyd->pmy_block->ncells2, phyd->pmy_block->ncells1),
  g2(phyd->pmy_block->ncells3, phyd->pmy_block->ncells2, phyd->pmy_block->ncells1),
  g3(phyd->pmy_block->ncells3, phyd->pmy_block->ncells2, phyd->pmy_block->ncells1),
  g1_defined_(false), g2_defined_(false), g3_defined_(false)
{
  hydro_sourceterms_defined = false;

  // read point mass or constant acceleration parameters from input block

  // set the point source only when the coordinate is spherical or 2D
  // cylindrical.
  gm_ = pin->GetOrAddReal("problem","GM",0.0);
  if (gm_ != 0.0) {
    if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0
        || (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0
            && phyd->pmy_block->block_size.nx3==1)) {
      hydro_sourceterms_defined = true;
      g1_defined_ = true;
    }
    else {
      std::stringstream msg;
      msg << "### FATAL ERROR in HydroSourceTerms constructor" << std::endl
          << "The point mass gravity works only in spherical polar coordinates"
          << "or in 2D cylindrical coordinates." << std::endl
          << "Check <problem> GM parameter in the input file." << std::endl;
      ATHENA_ERROR(msg);
    }
  }
  g1_ = pin->GetOrAddReal("hydro","grav_acc1",0.0);
  if (g1_ != 0.0) {
    hydro_sourceterms_defined = true;
    g1_defined_ = true;
  }

  g2_ = pin->GetOrAddReal("hydro","grav_acc2",0.0);
  if (g2_ != 0.0) {
    hydro_sourceterms_defined = true;
    g2_defined_ = true;
  }

  g3_ = pin->GetOrAddReal("hydro","grav_acc3",0.0);
  if (g3_ != 0.0) {
    hydro_sourceterms_defined = true;
    g3_defined_ = true;
  }

  // read shearing box parameters from input block
  Omega_0_ = pin->GetOrAddReal("problem","Omega0",0.0);
  qshear_  = pin->GetOrAddReal("problem","qshear",0.0);
  ShBoxCoord_ = pin->GetOrAddInteger("problem","shboxcoord",1);
  if ((Omega_0_ !=0.0) && (qshear_ != 0.0)) hydro_sourceterms_defined = true;

  if (SELF_GRAVITY_ENABLED) hydro_sourceterms_defined = true;

  UserSourceTerm = phyd->pmy_block->pmy_mesh->UserSourceTerm_;
  if (UserSourceTerm != nullptr) hydro_sourceterms_defined = true;

  UserGravFunc = phyd->pmy_block->pmy_mesh->UserGravFunc_;
  if (UserGravFunc != nullptr) {
    if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0)
      hydro_sourceterms_defined = true;
    else {
      std::stringstream msg;
      msg << "### FATAL ERROR in HydroSourceTerms constructor" << std::endl
          << "User defined gravity only works in cylindrical coordinates at the moment"
          << std::endl
          << "Verify it will work in your geometry before proceeding." << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  SetGravityStores();
}

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::setG()
//  \brief Sets g array on construction. Makes sure it is only set once

void HydroSourceTerms::SetGravityStores() {
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  int is = pmb->is, ie = pmb->ie;

  if (UserGravFunc != nullptr) {
    // Deal with user custom gravity here.

    // Throw error if both user gravity and other gravity.
    if (g1_ != 0 || g2_ != 0 || g3_ !=0 || gm_ != 0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in HydroSourceTerms::SetGravityStores" << std::endl
        << "Cannot define User Gravity and ConstantAcceleration/PointMass gravity"
        << std::endl;
      ATHENA_ERROR(msg);
    }

    // Gravity Stores set in user function
    UserGravFunc(pmb, g1, g2, g3);

    // set g_defined_ variables
    g1_defined_ = GravityDefined(1, pmb->ks, pmb->ke, pmb->js, pmb->je, pmb->is, pmb->ie);
    g2_defined_ = GravityDefined(2, pmb->ks, pmb->ke, pmb->js, pmb->je, pmb->is, pmb->ie);
    g3_defined_ = GravityDefined(3, pmb->ks, pmb->ke, pmb->js, pmb->je, pmb->is, pmb->ie);
  } // if
  else {
    // loop through meshblock, set gravity stores appropriately
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (gm_ != 0.0) {
            Real g = pmb->pcoord->coord_src1_i_(i)*gm_/pmb->pcoord->x1v(i);
            g1(k, j, i) += g;
          }
          if (g1_ != 0.0) {
            g1(k, j, i) += g1_;
          }
          if (g2_ != 0.0) {
            g2(k, j, i) += g2_;
          }
          if (g3_ != 0.0) {
            g3(k, j, i) += g3_;
          }
        } // i loop
      }
    }
  } // else

  // Boundary Conditions
  // only important in i direction at the moment
  // relevant in Hydro::DecomposePressurse()
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=1; i<=NGHOST; ++i) {
        // inner x_1
        if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::reflect)
          g1(k, j, is-i) = -1.0*g1(k, j, is+i);
        else if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::outflow)
          g1(k, j, is-i) = 0.0;
        else // block boundary
          g1(k, j, is-i) = g1(k, j, is);

        // outer_x1
        if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::reflect)
          g1(k, j, ie+i) = -1.0*g1(k, j, ie-i);
        else if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::outflow)
          g1(k, j, ie+i) = 0.0;
        else // block boundary
          g1(k, j, ie+i) = g1(k, j, ie);
      }
    }
  } // k loop
}

//----------------------------------------------------------------------------------------
//! \fn Real HydroSourceTerms::GetG#()
//  \brief Gets grav acc in directions

Real HydroSourceTerms::GetG1(int k, int j, int i) const {
  return g1(k, j, i);
}

Real HydroSourceTerms::GetG2(int k, int j, int i) const {
  return g2(k, j, i);
}

Real HydroSourceTerms::GetG3(int k, int j, int i) const {
  return g3(k, j, i);
}

//----------------------------------------------------------------------------------------
//! \fn bool HydroSourceTerms::GravityDefined()
//  \brief returns true if any vlaue of g in that direction is non zero

// if bounds are given, checks in those bounds
bool HydroSourceTerms::GravityDefined(int direction, int kl, int ku, int jl, int ju, int il, int iu) const {
  // interpret direction
  const AthenaArray<Real> *g;
  if (direction == 1) g = &g1;
  else if (direction == 2) g = &g2;
  else if (direction == 3) g = &g3;
  else {
    std::stringstream msg;
    msg << "### FATAL ERROR in HydroSourceTerms::GravityDefined" << std::endl
        << "int direction (" << direction << ") must be either 1, 2, or 3.";
    ATHENA_ERROR(msg);
  }

  // check for non-zero g values
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j)  {
      for (int i=il; i<=iu; ++i) {
        if ((*g)(k,j,i) != 0.0) return true;
      }
    }
  }
  return false;
}

bool HydroSourceTerms::GravityDefined(int kl, int ku, int jl, int ju, int il, int iu) const {
  // check for non-zero g values
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j)  {
      for (int i=il; i<=iu; ++i) {
        if ((g1(k,j,i) != 0.0) || (g2(k,j,i) != 0.0) || (g3(k,j,i) != 0.0)) return true;
      }
    }
  }
  return false;
}

// if no bounds are given, return g_defined_ variables
bool HydroSourceTerms::GravityDefined(int direction) const {
  if (direction == 1) return g1_defined_;
  else if (direction == 2) return g2_defined_;
  else if (direction == 3) return g3_defined_;
  else {
    std::stringstream msg;
    msg << "### FATAL ERROR in HydroSourceTerms::GravityDefined" << std::endl
        << "int direction (" << direction << ") must be either 1, 2, or 3.";
    ATHENA_ERROR(msg);
  }
}

bool HydroSourceTerms::GravityDefined() const {
  return (g1_defined_ || g2_defined_ || g3_defined_);
}

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::AddHydroSourceTerms
//  \brief Adds source terms to conserved variables

void HydroSourceTerms::AddHydroSourceTerms(const Real time, const Real dt,
                                           const AthenaArray<Real> *flux,
                                           const AthenaArray<Real> &prim,
                                           const AthenaArray<Real> &bcc,
                                           AthenaArray<Real> &cons) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;

  // accleration due to point mass (MUST BE AT ORIGIN)
  if (gm_ != 0.0) PointMass(dt, flux, prim, cons);

  // constant acceleration (e.g. for RT instability)
  if (g1_ != 0.0 || g2_ != 0.0 || g3_ != 0.0) ConstantAcceleration(dt, flux,
                                                                   prim,cons);

  // Add new source terms here
  if (SELF_GRAVITY_ENABLED) SelfGravity(dt, flux, prim, cons);
  // shearing box source terms: tidal and Coriolis forces
  if ((Omega_0_ !=0.0) && (qshear_ != 0.0)) ShearingBoxSourceTerms(dt, flux,
                                                                   prim, cons);
  // MyNewSourceTerms()

  // user-defined gravity term
  if (UserGravFunc != nullptr)
    AddUserGravityTerm(dt, flux, prim, cons);

  //  user-defined source terms
  if (UserSourceTerm != nullptr)
    UserSourceTerm(pmb, time, dt, prim, bcc, cons);

  return;
}
