#ifndef HYDRO_SRCTERMS_HYDRO_SRCTERMS_HPP_
#define HYDRO_SRCTERMS_HYDRO_SRCTERMS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydro_srcterms.hpp
//  \brief defines class HydroSourceTerms
//  Contains data and functions that implement physical (not coordinate) source terms

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"

// Forward declarations
class Hydro;
class ParameterInput;

//! \class HydroSourceTerms
//  \brief data and functions for physical source terms in the hydro

class HydroSourceTerms {
 public:
  HydroSourceTerms(Hydro *phyd, ParameterInput *pin);

  // Stores of gravitational acceleration
  AthenaArray<Real> g1, g2, g3;

  // accessors
  Real GetGM() const {return gm_;}
  Real GetG1() const {return g1_;}
  Real GetG2() const {return g2_;}
  Real GetG3() const {return g3_;}

  Real GetG1(int k, int j, int i) const;
  Real GetG2(int k, int j, int i) const;
  Real GetG3(int k, int j, int i) const;

  // returns true if (directional) gravity non-zero anywhere within these bounds
  bool GravityDefined(int direction, int kl, int ku, int jl, int ju, int il, int iu) const;
  bool GravityDefined(int kl, int ku, int jl, int ju, int il, int iu) const;
  bool GravityDefined(int direction) const;
  bool GravityDefined() const;

  // data
  bool hydro_sourceterms_defined;

  // functions
  void AddSourceTerms(const Real time, const Real dt, const AthenaArray<Real> *flx,
                           const AthenaArray<Real> &prim,
                           const AthenaArray<Real> &prim_scalar,
                           const AthenaArray<Real> &b, AthenaArray<Real> &cons,
                           AthenaArray<Real> &cons_scalar);
  void SetGravityStores();
  void PointMass(const Real dt, const AthenaArray<Real> *flx,const AthenaArray<Real> &p,
                 AthenaArray<Real> &c);
  void ConstantAcceleration(const Real dt, const AthenaArray<Real> *flx,
                            const AthenaArray<Real> &p, AthenaArray<Real> &c);
  void AddUserGravityTerm(const Real dt, const AthenaArray<Real> *flx,
                  const AthenaArray<Real> &p, AthenaArray<Real> &c);
  // shearing box src terms
  void ShearingBoxSourceTerms(const Real dt, const AthenaArray<Real> *flx,
                              const AthenaArray<Real> &p, AthenaArray<Real> &c);
  Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3);

  void SelfGravity(const Real dt, const AthenaArray<Real> *flx,
                   const AthenaArray<Real> &p, AthenaArray<Real> &c);
  void EnrollSrcTermFunction(SrcTermFunc my_func);
  SrcTermFunc UserSourceTerm;
  GravFunc UserGravFunc;

 private:
  Hydro *pmy_hydro_;  // ptr to Hydro containing this HydroSourceTerms
  Real gm_;           // GM for point mass MUST BE LOCATED AT ORIGIN
  Real g1_, g2_, g3_; // constant acc'n in each direction
  bool g1_defined_, g2_defined_, g3_defined_; // bool for any gravity set in direction
  Real Omega_0_, qshear_; // Orbital freq and shear rate
  int  ShBoxCoord_;       // ShearCoordinate type: 1=xy (default), 2=xz
};
#endif // HYDRO_SRCTERMS_HYDRO_SRCTERMS_HPP_
