// C/C++ header
#include <vector>  // fill
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <stdexcept>

// athena header
#include "thermodynamics.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
//#include "../math_funcs.hpp"

std::ostream& operator<<(std::ostream& os, Thermodynamics const& my)
{
  os << "Rd [J/kg]: " << my.Rd_ << std::endl;
  os << "eps: ";
  for (int i = 0; i < NMASS; ++i)
    os << my.eps_[i] << " ";
  os << std::endl;
  os << "rcp: ";
  for (int i = 0; i < NMASS; ++i)
    os << my.rcp_[i] << " ";
  os << std::endl;
  os << "rcv: ";
  for (int i = 0; i < NMASS; ++i)
    os << my.rcv_[i] << " ";
  os << std::endl;
  os << "beta: ";
  for (int i = 0; i < NMASS; ++i)
    os << my.beta_[i] << " ";
  os << std::endl;
  os << "delta: ";
  for (int i = 0; i < NMASS; ++i)
    os << my.delta_[i] << " ";
  os << std::endl;
  os << "Ttriple [K]: ";
  for (int i = 0; i < NMASS; ++i)
    os << my.t3_[i] << " ";
  os << std::endl;
  os << "Ptriple [pa]: ";
  for (int i = 0; i < NMASS; ++i)
    os << my.p3_[i] << " ";
  os << std::endl;
  os << "Latent [J/kg]: ";
  for (int i = 0; i < NMASS; ++i)
    os << my.latent_[i] << " ";

  return os;
}

void ReadThermoProperty(Real var[], char const name[], int len, Real v0, ParameterInput *pin)
{
  std::stringstream msg;
  char buf[80], cstr[1024], *p;

  var[0] = v0;
  for (int n = 1; n <= NVAPOR; ++n) {
    sprintf(buf, "%s%d", name, n);
    std::string str = pin->GetString("thermodynamics", buf);
    std::strcpy(cstr, str.c_str());
    p = std::strtok(cstr," ,");
    int m = 0;
    while ((p != NULL) && (m++ < len)) {
      var[n+(m-1)*NVAPOR] = std::stod(p);
      p = std::strtok(NULL," ,");
    }
    if (m != len) {
      msg << "### FATAL ERROR in function ReadThermoProperty"
          << std::endl << "Length of '" << name << "' "
          << "doesn't equal to " << len;
      ATHENA_ERROR(msg);
    }
  }
}

Thermodynamics::Thermodynamics(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;
  Rd_ = pin->GetOrAddReal("thermodynamics", "Rd", 1.);

  Real gamma = pin->GetReal("hydro", "gamma");

  ReadThermoProperty(eps_, "eps", NPHASE, 1., pin);
  ReadThermoProperty(rcp_, "rcp", NPHASE, 1., pin);
  ReadThermoProperty(beta_, "beta", NPHASE, 0., pin);
  ReadThermoProperty(t3_, "Ttriple", NPHASE, 0., pin);
  ReadThermoProperty(p3_, "Ptriple", NPHASE, 0., pin);

  // calculate latent
  for (int n = 0; n <= NVAPOR; ++n)
    latent_[n] = 0.;
  for (int n = 1 + NVAPOR; n < NMASS; ++n)
    latent_[n] = beta_[n]*Rd_/eps_[n]*t3_[n];

  // calculate delta
  for (int n = 0; n <= NVAPOR; ++n)
    delta_[n] = 0.;
  for (int n = 1 + NVAPOR; n < NMASS; ++n)
    delta_[n] = (rcp_[n] - rcp_[1+(n-1)%NVAPOR])*eps_[n]/(1. - 1./gamma);

  // calculate rcv
  for (int n = 0; n <= NVAPOR; ++n)
    rcv_[n] = gamma*rcp_[n] + (1. - gamma)/eps_[n];
  for (int n = 1 + NVAPOR; n < NMASS; ++n)
    rcv_[n] = rcp_[n]*gamma;

  // allocate temperature
  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
  if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);

  //temp_.NewAthenaArray(ncells3, ncells2, ncells1);
  //pres_.NewAthenaArray(ncells3, ncells2, ncells1);

  // last time that executes SaturationAdjustment
  last_time = 0.;

  // minimum time between two SaturationAdjustment
  dt = pin->GetOrAddReal("thermodynamics", "dt", 0.1);
  ftol_ = pin->GetOrAddReal("thermodynamics", "ftol", 1.0E-4);
  max_iter_ = pin->GetOrAddInteger("thermodynamics", "max_iter", 10);
}

Thermodynamics::~Thermodynamics() {}

void Thermodynamics::UpdateTPConservingU(Real q[], Real rho, Real uhat) const
{
  Real gamma = pmy_block_->peos->GetGamma();
  Real cv = 1., qtol = 1., qeps = 1.;
  for (int n = 1 + NVAPOR; n < NMASS; ++n) {
    uhat += beta_[n]*t3_[n]*q[n];
    qtol -= q[n];
  }
  for (int n = 1; n < NMASS; ++n) {
    cv += (rcv_[n]*eps_[n] - 1.)*q[n];
    qeps += q[n]*(eps_[n] - 1.);
  }
  q[IDN] = (gamma - 1.)*uhat/cv;
  q[IPR] = rho*Rd_*q[IDN]*qtol/qeps;
}

void Thermodynamics::ConservedToThermodynamic(AthenaArray<Real> &q, 
  AthenaArray<Real> const& u, int il, int iu, int jl, int ju, int kl, int ku) const
{
  Real gamma = pmy_block_->peos->GetGamma();

  for (int k=kl; k<=ku; ++k)
    for (int j=jl; j<=ju; ++j)
      for (int i=il; i<=iu; ++i) {
        // density to molar mixing ratio
        Real rho = 0., rho_hat = 0.;
        for (int n = 0; n < NMASS; ++n) {
          rho += u(n,k,j,i);
          rho_hat += u(n,k,j,i)/eps_[n];
          q(n,k,j,i) = u(n,k,j,i)/eps_[n];
        }
        for (int n = 0; n < NMASS; ++n)
          q(n,k,j,i) /= rho_hat;

        // calculate internal energy
        Real KE = 0.5*(u(IM1,k,j,i)*u(IM1,k,j,i)
                     + u(IM2,k,j,i)*u(IM2,k,j,i)
                     + u(IM3,k,j,i)*u(IM3,k,j,i))/rho;
        Real uhat = (u(IEN,k,j,i) - KE)/(Rd_*rho_hat);

        // calculate temperature and pressure
        Real cv = 1., qtol = 1., qeps = 1.;
        for (int n = 1 + NVAPOR; n < NMASS; ++n) {
          uhat += beta_[n]*t3_[n]*q(n,k,j,i);
          qtol -= q(n,k,j,i);
        }
        for (int n = 1; n < NMASS; ++n) {
          cv += (rcv_[n]*eps_[n] - 1.)*q(n,k,j,i);
          qeps += q(n,k,j,i)*(eps_[n] - 1.);
        }
        q(IDN,k,j,i) = (gamma - 1.)*uhat/cv;
        q(IPR,k,j,i) = rho*Rd_*q(IDN,k,j,i)*qtol/qeps;
      }
}

void Thermodynamics::ThermodynamicToConserved(AthenaArray<Real> &u,
  AthenaArray<Real> const& q, int il, int iu, int jl, int ju, int kl, int ku) const
{
  for (int k=kl; k<=ku; ++k)
    for (int j=jl; j<=ju; ++j)
      for (int i=il; i<=iu; ++i) {
        // calculate total density
        Real qtol = 1., qeps = 1.;
        for (int n = 1 + NVAPOR; n < NMASS; ++n)
          qtol -= q(n,k,j,i);
        for (int n = 1; n < NMASS; ++n)
          qeps += q(n,k,j,i)*(eps_[n] - 1.);
        Real rho = (q(IPR,k,j,i)*qeps)/(Rd_*q(IDN,k,j,i)*qtol);

        // molar mixing ratio to density
        Real sum = 1.;
        for (int n = 1; n < NMASS; ++n)
          sum += q(n,k,j,i)*(eps_[n] - 1.);
        for (int n = 1; n < NMASS; ++n)
          u(n,k,j,i) = rho*q(n,k,j,i)*eps_[n]/sum;
      }
}

void Thermodynamics::PolytropicIndex(AthenaArray<Real> &gm, AthenaArray<Real> &w,
  int kl, int ku, int jl, int ju, int il, int iu) const
{
  Real gamma = pmy_block_->peos->GetGamma();

  // calculate local polytropic index
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i) {
        Real fsig = 1., feps = 1.;
        for (int n = 1 + NVAPOR; n < NMASS; ++n) {
          fsig += w(n,k,j,i)*(rcv_[n] - 1.);
          feps -= w(n,k,j,i);
        }
        for (int n = 1; n <= NVAPOR; ++n) {
          fsig += w(n,k,j,i)*(rcv_[n] - 1.);
          feps += w(n,k,j,i)*(1./eps_[n] - 1.);
        }
        gm(k,j,i) = 1. + (gamma - 1.)*feps/fsig;
      }
}
