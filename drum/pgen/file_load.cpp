// C/C++ header
#include <ctime>
#include <cmath>

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
#include "../radiation/radiation.hpp"
#include "../radiation/hydrogen_cia.hpp"
#include "../radiation/freedman_mean.hpp"

// global parameters
Real user_grav;
//ics
Real T_iso, rho_ambient, p_ambient;
//gravity
Real G, Mp, Ms, Rp, period, a;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "g1");
  // SetUserOutputVariableName(2, "ghost_rho");
  // SetUserOutputVariableName(3, "ghost_press");
  // SetUserOutputVariableName(4, "ghost_temp");
  // SetUserOutputVariableName(5, "ghost_g1");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  Real dq[1+NVAPOR], rh;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is-NGHOST; i <= ie+NGHOST; ++i) {
        user_out_var(0,k,j,i) = pthermo->Temp(phydro->w.at(k,j,i));
        user_out_var(1,k,j,i) = phydro->hsrc.g1(k,j,i);

      //   user_out_var(2,k,j,i) = 0.0;
      //   user_out_var(3,k,j,i) = 0.0;
      //   user_out_var(4,k,j,i) = 0.0;
      //   user_out_var(5,k,j,i) = 0.0;

      //   // low side ghosts
      //   if (i - NGHOST <= is) {
      //     user_out_var(2,k,j,i) = phydro->w(IDN,k,j,i-NGHOST);
      //     user_out_var(3,k,j,i) = phydro->w(IPR,k,j,i-NGHOST);
      //     user_out_var(4,k,j,i) = pthermo->Temp(phydro->w.at(k,j,i-NGHOST));
      //     user_out_var(5,k,j,i) = phydro->hsrc.g1(k,j,i-NGHOST);
      //   }

      //   if (i + NGHOST >= ie) {
      //     user_out_var(2,k,j,i) = phydro->w(IDN,k,j,i+NGHOST);
      //     user_out_var(3,k,j,i) = phydro->w(IPR,k,j,i+NGHOST);
      //     user_out_var(4,k,j,i) = pthermo->Temp(phydro->w.at(k,j,i+NGHOST));
      //     user_out_var(5,k,j,i) = phydro->hsrc.g1(k,j,i+NGHOST);
        // }
      }
}

void RadiationBand::AddAbsorber(std::string name, std::string file, ParameterInput *pin)
{
  Real xHe = pin->GetOrAddReal("radiation", "xHe", 0.136);
  Real xH2 = 1. - xHe;

  std::stringstream msg;

  if (name == "H2-H2") {
    pabs->AddAbsorber(XizH2H2CIA(this, 0, xH2))
        ->LoadCoefficient(file);
  } else if (name == "H2-He") {
    pabs->AddAbsorber(XizH2HeCIA(this, 0, xH2, xHe))
        ->LoadCoefficient(file);
  } else if (name == "FREEDMAN") {
    pabs->AddAbsorber(FreedmanMean(this));
  } else {
    msg << "### FATAL ERROR in RadiationBand::AddAbsorber"
        << std::endl << "unknow absorber: '" << name <<"' ";
    ATHENA_ERROR(msg);
  }
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
    AthenaArray<Real> const &w, AthenaArray<Real> const &bcc, AthenaArray<Real> &u)
{ 
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        // damp kinetic energy
        u(IM1,k,j,i) -= w(IDN,k,j,i)*w(IM1,k,j,i)/5.;
        u(IM2,k,j,i) -= w(IDN,k,j,i)*w(IM2,k,j,i)/5.;
        u(IM3,k,j,i) -= w(IDN,k,j,i)*w(IM3,k,j,i)/5.;
      }
}

void Gravity(MeshBlock *pmb, AthenaArray<Real> &g1, AthenaArray<Real> &g2, AthenaArray<Real> &g3) {
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real x,y,z, rad, rad_sq;
  Real gc1, gc2, gc3;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is-NGHOST; i<=ie+NGHOST; ++i) {
        x = pmb->pcoord->x1v(i);
        y = pmb->pcoord->x2v(j);
        z = pmb->pcoord->x3v(k);

        rad_sq = x*x + y*y + z*z;
        rad = sqrt(rad_sq);

        gc1 = - (G * Mp) / (rad_sq);
        gc2 = (G * Ms) / pow((a-rad),2);
        gc3 = - (G * Ms) * (a-rad) / pow(a,3);

        g1(k, j, i) = gc1+gc2+gc3;
      }
    }
  }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // EnrollUserExplicitSourceFunction(Forcing);

  Real hydro_grav = pin->GetOrAddReal("hydro", "grav_acc1",0.0);
  if (hydro_grav != 0.0) {
    std::cout << std::endl << std::endl << "----------------------" << std::endl << "Using constant gravity" << std::endl << "----------------------" <<std::endl << std::endl << std::endl;
  }
  else {
    EnrollUserExplicitGravityFunction(Gravity);
  }

  if (pin->GetOrAddBoolean("problem","forcing_flag",true))
    EnrollUserExplicitSourceFunction(Forcing);
  // Real G, Mp, Ms, Rp, period, a;
  G = pin->GetReal("problem","G");
  Mp = pin->GetReal("problem","Mp");
  Ms = pin->GetReal("problem","Ms");
  Rp = pin->GetReal("problem","Rp");
  period = pin->GetReal("problem","period");

  Real x = 4 * pow(M_PI,2) / (G * Ms);
  a = pow( pow(period*86400,2)/x ,(1./3));

  rho_ambient = pin->GetOrAddReal("problem","rho_ambient", 0.0);
  T_iso       = pin->GetOrAddReal("problem","T_iso",0.0);

  if (rho_ambient != 0.0) {
    Real R = pin->GetReal("thermodynamics","Rd");
    p_ambient = rho_ambient * R * T_iso;
  }
  else {
    p_ambient = 0.0;
  }

}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  int nx1 = pmy_mesh->mesh_size.nx1;

  if (block_size.nx1 != nx1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in problem generator" << std::endl
        << "Currently file_load is not configured for multiple meshblocks." << std::endl;
    ATHENA_ERROR(msg);
  }

  if (pmy_mesh->f2 || pmy_mesh->f3) {
    std::stringstream msg;
    msg << "### FATAL ERROR in problem generator" << std::endl
        << "File Load is not configured for multiple dimensions." << std::endl;
    ATHENA_ERROR(msg);
  }

  std::string filename;
  filename = pin->GetString("problem", "initial_conditions_file");

  int numRows = GetNumRows(filename.c_str());
  int il, iu;
  bool include_ghost_zones;

  // check size of input file is compatible
  if (numRows == nx1+1) {
    il = is;
    iu = ie;
    include_ghost_zones = false;

    if (pbval->block_bcs[inner_x1] == BoundaryFlag::outflow ||
        pbval->block_bcs[outer_x1] == BoundaryFlag::outflow) {
      std::stringstream msg;
      msg << "### FATAL ERROR in problem generator" << std::endl
          << "    outflow boundary conditions require initial conditions to be set in ghost zones" << std::endl;
      ATHENA_ERROR(msg);
    }

  }
  else if (numRows == nx1 + 2*NGHOST +1) {
    il = is - NGHOST;
    iu = ie + NGHOST;
    include_ghost_zones = true;
  }
  else {
    std::stringstream msg;
    msg << "### FATAL ERROR in problem generator" << std::endl
        << "Input file dimension mismatch with nx1." << std::endl
        << "Input file number of rows must match nx1+1 (no ghost cells)" << std::endl
        << "Or nx1 + 2*NGHOST+1 (with ghost cells)" << std::endl
        << "with one row for header" << std::endl;
    ATHENA_ERROR(msg);
  }

  std::ifstream inp(filename.c_str(), std::ifstream::in);
  std::string line;

  // discard header
  std::getline(inp, line);

  Real *V1 = new Real [numRows-1];
  Real *V2 = new Real [numRows-1];
  Real *V3 = new Real [numRows-1];

  Real *DN = new Real [numRows-1];
  Real *PR = new Real [numRows-1];

  int i = 0;
  while (std::getline(inp, line)) {
    Real dn_i, pr_i, v1_i, v2_i, v3_i;
    sscanf(line.c_str(), "%le %le %le %le %le\n", &dn_i, &pr_i, &v1_i, &v2_i, &v3_i);
    V1[i] = v1_i;
    V2[i] = v2_i;
    V3[i] = v3_i;

    DN[i] = dn_i;
    PR[i] = pr_i;

    i+=1;
  }

  // setup initial condition
  int kl = block_size.nx3 == 1 ? ks : ks-NGHOST;
  int ku = block_size.nx3 == 1 ? ke : ke+NGHOST;
  int jl = block_size.nx2 == 1 ? js : js-NGHOST;
  int ju = block_size.nx2 == 1 ? je : je+NGHOST;

  for (int i = il; i <= iu; ++i) {
    for (int k = kl; k <= ku; ++k) {
      for (int j = jl; j <= ju; ++j) {
        int ii;
        // offset i to match start index depending on ghost include or not
        if (include_ghost_zones) {
          ii = i;
        }
        else {
          ii = i - NGHOST;
        }

        phydro->w(IPR,k,j,i) = PR[ii] + p_ambient;
        phydro->w(IDN,k,j,i) = DN[ii] + rho_ambient;
        phydro->w(IV1,k,j,i) = V1[ii];
        phydro->w(IV2,k,j,i) = V2[ii];
        phydro->w(IV3,k,j,i) = V3[ii];
      }
    }
  }

  // set spectral properties
  RadiationBand *p = prad->pband;
  while (p != NULL) {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        p->SetSpectralProperties(phydro->w, k, j, il, iu);
    p = p->next;
  }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);

  delete[] V1;
  delete[] V2;
  delete[] V3;
  delete[] PR;
  delete[] DN;
}
