
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
#include "../radiation/radiation.hpp"
#include "../radiation/hydrogen_cia.hpp"
#include "../radiation/freedman_mean.hpp"
#include "../radiation/user_absorber.hpp"

Real MyMeshSpacingX1(Real x, RegionSize rs) {
  Real min = rs.x1min;
  Real width = rs.x1max - min;
  Real t1 = (2. * x - 1.)*(2. * x - 1.)*(2. * x - 1.) + 1. + 3.*x;
  return t1/5. * width + min;
}

Real MyMeshSpacingX2(Real x, RegionSize rs) {
  Real min = rs.x2min;
  Real width = rs.x2max - min;
  Real t1 = (2. * x - 1.)*(2. * x - 1.)*(2. * x - 1.) + 1. + 3.*x;
  return t1/5. * width + min;
}

Real MyMeshSpacingX3(Real x, RegionSize rs) {
  Real min = rs.x3min;
  Real width = rs.x3max - min;
  Real t1 = (2. * x - 1.)*(2. * x - 1.)*(2. * x - 1.) + 1. + 3.*x;
  return t1/5. * width + min;
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (mesh_size.x1rat == -1.0)
    EnrollUserMeshGenerator(X1DIR, MyMeshSpacingX1);

  if (f2) {
    if (mesh_size.x2rat == -1.0)
      EnrollUserMeshGenerator(X2DIR, MyMeshSpacingX2);
  }
  if (f3) {
    if (mesh_size.x3rat == -1.0)
      EnrollUserMeshGenerator(X3DIR, MyMeshSpacingX3);
  }

  return;
}

// void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
// {
//   AllocateUserOutputVariables(2);
//   SetUserOutputVariableName(0,"rad_eng");
//   SetUserOutputVariableName(1,"rad_eng");
// }

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Spherical blast wave test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rout = pin->GetReal("problem", "radius");
  Real rin  = rout - pin->GetOrAddReal("problem", "ramp", 0.0);
  Real pa   = pin->GetOrAddReal("problem", "pamb", 1.0);
  Real da   = pin->GetOrAddReal("problem", "damb", 1.0);
  Real prat = pin->GetReal("problem", "prat");
  Real drat = pin->GetOrAddReal("problem", "drat", 1.0);
  Real b0, angle;

  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;

  // get coordinates of center of blast, and convert to Cartesian if necessary
  Real x1_0   = pin->GetOrAddReal("problem", "x1_0", 0.0);
  Real x2_0   = pin->GetOrAddReal("problem", "x2_0", 0.0);
  Real x3_0   = pin->GetOrAddReal("problem", "x3_0", 0.0);
  Real x0, y0, z0;

  x0 = x1_0;
  y0 = x2_0;
  z0 = x3_0;

  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real rad;
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));

        Real den = da;
        if (rad < rout) {
          if (rad < rin) {
            den = drat*da;
          } else {   // add smooth ramp in density
            Real f = (rad-rin) / (rout-rin);
            Real log_den = (1.0-f) * std::log(drat*da) + f * std::log(da);
            den = std::exp(log_den);
          }
        }

        phydro->u(IDN,k,j,i) = den;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          Real pres = pa;
          if (rad < rout) {
            if (rad < rin) {
              pres = prat*pa;
            } else {  // add smooth ramp in pressure
              Real f = (rad-rin) / (rout-rin);
              Real log_pres = (1.0-f) * std::log(prat*pa) + f * std::log(pa);
              pres = std::exp(log_pres);
            }
          }
          phydro->u(IEN,k,j,i) = pres/gm1;
        }
      }
    }
  }

  peos->ConservedToPrimitive(phydro->u, phydro->w1, pfield->b, phydro->w, pfield->bcc, pcoord, is, ie, js, je, ks, ke);
  // set spectral properties
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      prad->CalculateRadiativeTransfer(phydro->w, pmy_mesh->time, k, j, is, ie+1);
}

Real AbsorptionCoefficient(Absorber const *pabs, Real wave, Real const prim[], int k, int j, int i)
{
  Hydro *ph = pabs->pmy_band->pmy_rad->pmy_block->phydro;
  Real n = ph->u(IDN,k,j,i);

  // sigma = 1
  return n*10; // 1/m
}

// makes shit slow!
Real EnergyAbsorption(Absorber *pabs, Real wave, Real flux, int k, int j, int i)
{
  Real energy_fraction = 0.15;
  // pabs->pmy_band->pmy_rad->pmy_block->ruser_meshblock_data[0](k,j,i) += (1-energy_fraction) * flux;
  return flux*energy_fraction;
}

void RadiationBand::AddAbsorber(std::string name, std::string file, ParameterInput *pin)
{ 
  std::stringstream msg;
  if (name == "HYDROGEN_IONIZATION") {
    UserDefinedAbsorber hydrogen_ionization(this);

    hydrogen_ionization.EnrollUserAbsorptionCoefficientFunc(AbsorptionCoefficient);
    hydrogen_ionization.EnrollUserEnergyAbsorptionFunc(EnergyAbsorption);
    pabs->AddAbsorber(hydrogen_ionization);
  } else {
    msg << "### FATAL ERROR in RadiationBand::AddAbsorber"
        << std::endl << "unknown absorber: '" << name <<"' ";
    ATHENA_ERROR(msg);
  }
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief Check radius of sphere to make sure it is round
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  if (!pin->GetOrAddBoolean("problem","compute_error",false)) return;

  // analysis - check shape of the spherical blast wave
  int is = pblock->is, ie = pblock->ie;
  int js = pblock->js, je = pblock->je;
  int ks = pblock->ks, ke = pblock->ke;
  AthenaArray<Real> pr;
  pr.InitWithShallowSlice(pblock->phydro->w, 4, IPR, 1);

  // get coordinate location of the center, convert to Cartesian
  Real x1_0 = pin->GetOrAddReal("problem", "x1_0", 0.0);
  Real x2_0 = pin->GetOrAddReal("problem", "x2_0", 0.0);
  Real x3_0 = pin->GetOrAddReal("problem", "x3_0", 0.0);
  Real x0, y0, z0;
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    x0 = x1_0*std::cos(x2_0);
    y0 = x1_0*std::sin(x2_0);
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z0 = x1_0*std::cos(x2_0);
  } else {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::stringstream msg;
    msg << "### FATAL ERROR in blast.cpp ParameterInput" << std::endl
        << "Unrecognized COORDINATE_SYSTEM= " << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }

  // find indices of the center
  int ic, jc, kc;
  for (ic=is; ic<=ie; ic++)
    if (pblock->pcoord->x1f(ic) > x1_0) break;
  ic--;
  for (jc=pblock->js; jc<=pblock->je; jc++)
    if (pblock->pcoord->x2f(jc) > x2_0) break;
  jc--;
  for (kc=pblock->ks; kc<=pblock->ke; kc++)
    if (pblock->pcoord->x3f(kc) > x3_0) break;
  kc--;

  // search pressure maximum in each direction
  Real rmax = 0.0, rmin = 100.0, rave = 0.0;
  int nr = 0;
  for (int o=0; o<=6; o++) {
    int ios = 0, jos = 0, kos = 0;
    if (o == 1) ios=-10;
    else if (o == 2) ios =  10;
    else if (o == 3) jos = -10;
    else if (o == 4) jos =  10;
    else if (o == 5) kos = -10;
    else if (o == 6) kos =  10;
    for (int d=0; d<6; d++) {
      Real pmax = 0.0;
      int imax(0), jmax(0), kmax(0);
      if (d == 0) {
        if (ios != 0) continue;
        jmax = jc+jos, kmax = kc+kos;
        for (int i=ic; i>=is; i--) {
          if (pr(kmax,jmax,i)>pmax) {
            pmax = pr(kmax,jmax,i);
            imax = i;
          }
        }
      } else if (d == 1) {
        if (ios != 0) continue;
        jmax = jc+jos, kmax = kc+kos;
        for (int i=ic; i<=ie; i++) {
          if (pr(kmax,jmax,i)>pmax) {
            pmax = pr(kmax,jmax,i);
            imax = i;
          }
        }
      } else if (d == 2) {
        if (jos != 0) continue;
        imax = ic+ios, kmax = kc+kos;
        for (int j=jc; j>=js; j--) {
          if (pr(kmax,j,imax)>pmax) {
            pmax = pr(kmax,j,imax);
            jmax = j;
          }
        }
      } else if (d == 3) {
        if (jos != 0) continue;
        imax = ic+ios, kmax = kc+kos;
        for (int j=jc; j<=je; j++) {
          if (pr(kmax,j,imax)>pmax) {
            pmax = pr(kmax,j,imax);
            jmax = j;
          }
        }
      } else if (d == 4) {
        if (kos != 0) continue;
        imax = ic+ios, jmax = jc+jos;
        for (int k=kc; k>=ks; k--) {
          if (pr(k,jmax,imax)>pmax) {
            pmax = pr(k,jmax,imax);
            kmax = k;
          }
        }
      } else { // if (d == 5) {
        if (kos != 0) continue;
        imax = ic+ios, jmax = jc+jos;
        for (int k=kc; k<=ke; k++) {
          if (pr(k,jmax,imax)>pmax) {
            pmax = pr(k,jmax,imax);
            kmax = k;
          }
        }
      }

      Real xm, ym, zm;
      Real x1m = pblock->pcoord->x1v(imax);
      Real x2m = pblock->pcoord->x2v(jmax);
      Real x3m = pblock->pcoord->x3v(kmax);
      if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
        xm = x1m;
        ym = x2m;
        zm = x3m;
      } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
        xm = x1m*std::cos(x2m);
        ym = x1m*std::sin(x2m);
        zm = x3m;
      } else {  // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
        xm = x1m*std::sin(x2m)*std::cos(x3m);
        ym = x1m*std::sin(x2m)*std::sin(x3m);
        zm = x1m*std::cos(x2m);
      }
      Real rad = std::sqrt(SQR(xm-x0)+SQR(ym-y0)+SQR(zm-z0));
      if (rad > rmax) rmax = rad;
      if (rad < rmin) rmin = rad;
      rave += rad;
      nr++;
    }
  }
  rave /= static_cast<Real>(nr);

  // use physical grid spacing at center of blast
  Real dr_max;
  Real  x1c = pblock->pcoord->x1v(ic);
  Real dx1c = pblock->pcoord->dx1f(ic);
  Real  x2c = pblock->pcoord->x2v(jc);
  Real dx2c = pblock->pcoord->dx2f(jc);
  Real dx3c = pblock->pcoord->dx3f(kc);
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    dr_max = std::max(std::max(dx1c, dx2c), dx3c);
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    dr_max = std::max(std::max(dx1c, x1c*dx2c), dx3c);
  } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    dr_max = std::max(std::max(dx1c, x1c*dx2c), x1c*std::sin(x2c)*dx3c);
  }
  Real deform=(rmax-rmin)/dr_max;

  // only the root process outputs the data
  if (Globals::my_rank == 0) {
    std::string fname;
    fname.assign("blastwave-shape.dat");
    std::stringstream msg;
    FILE *pfile;

    // The file exists -- reopen the file in append mode
    if ((pfile = std::fopen(fname.c_str(),"r")) != nullptr) {
      if ((pfile = std::freopen(fname.c_str(),"a",pfile)) == nullptr) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Blast shape output file could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
      }

      // The file does not exist -- open the file in write mode and add headers
    } else {
      if ((pfile = std::fopen(fname.c_str(),"w")) == nullptr) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Blast shape output file could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
      }
    }
    std::fprintf(pfile,"# Offset blast wave test in %s coordinates:\n",COORDINATE_SYSTEM);
    std::fprintf(pfile,"# Rmax       Rmin       Rave        Deformation\n");
    std::fprintf(pfile,"%e  %e  %e  %e \n",rmax,rmin,rave,deform);
    std::fclose(pfile);
  }
  return;
}
