// C/C++ header
#include <ctime>
#include <cmath>
#include <string>

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
#include "../scalars/scalars.hpp"
#include "../mesh_generator.hpp"

#include "../radiation/radiation.hpp"
#include "../reaction/reaction_network.hpp"
#include "../radiation/absorber/hydrogen_ionization.hpp"
#include "../reaction/reactions/photoionization.hpp"
#include "../reaction/reactions/hydrogen_reactions.hpp"



namespace {
  // user set variables
  Real G, Mp, Ms, Rp, period, a;
  Real dfloor, pfloor, sfloor;
  Real gas_gamma;
  Real wave_to_meters_conversion;

  // radiation variables
  Real dist_, ref_dist_;

  // tripathi conditions variables
  Real rho_p, cs, space_density_factor;
  Real r_0, r_e, rho_0, rho_e, P_0, P_e;
  Real r_replenish;

  //Physical Constants -- for problem file?
  Real A0=6.30431812E-22; // (m^2)
  Real nu_0=3.2898419603E15; // (1/s) Rydberg frequency
  Real c=2.998E8; // m/s
  Real mh=1.674E-27; // kg
  Real nm_0=91.126705058; // Hydrogen ionization wavelength in nm
  Real Ry=2.1798723611E-18;// J (joules)

  enum species {ELEC = 0, HYD = 1, HPLUS = 2};
} //namespace

  // Function declarations
  Real getRad(Coordinates *pcoord, int i, int j, int k);
  Real rho_func(Real r);
  Real press_func(Real rho);
  void getMBBounds(MeshBlock *pmb, int &il, int &iu, int &jl, int &ju, int &kl, int &ku);
  void gravity_func(MeshBlock *pmb, AthenaArray<Real> &g1, AthenaArray<Real> &g2, AthenaArray<Real> &g3);

  void SourceTerms(MeshBlock *pmb, const Real time, const Real dt,
    const AthenaArray<Real> &w, const AthenaArray<Real> &r,
    const AthenaArray<Real> &bcc,
    AthenaArray<Real> &du, AthenaArray<Real> &ds);
  Real RadiationTime(AthenaArray<Real> const &prim, Real time, int k, int j);

  // mesh generators
  MeshGenerator meshgen_x1, meshgen_x2, meshgen_x3;
  Real MeshSpacingX1(Real x, RegionSize rs);
  Real MeshSpacingX2(Real x, RegionSize rs);
  Real MeshSpacingX3(Real x, RegionSize rs);


//----------------------------------------------------------------------------------------
// User Setup
//----------------------------------------------------------------------------------------
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserExplicitGravityFunction(gravity_func);
  EnrollUserExplicitSourceFunction(SourceTerms);

  // Radiation parameters
  dist_ = pin->GetOrAddReal("radiation", "distance", 1.);
  ref_dist_ = pin->GetOrAddReal("radiation", "reference_distance", 1.);

  wave_to_meters_conversion = pin->GetOrAddReal("radiation","wave_to_meters",1.e-9);
  EnrollUserRadiationScalingFunction(RadiationTime);

  // Gravity/System Parameters
  G = pin->GetReal("problem","G");
  Mp = pin->GetReal("problem","Mp");
  Ms = pin->GetReal("problem","Ms");
  Rp = pin->GetReal("problem","Rp");
  period = pin->GetReal("problem","period");

  Real x = 4. * pow(M_PI,2.) / (G * Ms);
  a = pow( pow(period*86400.,2.)/x ,(1./3));

  gas_gamma = pin->GetReal("hydro","gamma");

  dfloor = pin->GetOrAddReal("hydro", "dfloor", 0);
  pfloor = pin->GetOrAddReal("hydro", "pfloor", 0);
  sfloor = pin->GetOrAddReal("hydro", "sfloor", 0);

  rho_p = pin->GetOrAddReal("problem","rho_p",1.0e-12);
  cs = pin->GetOrAddReal("problem","cs",3.0e3);
  space_density_factor = pin->GetOrAddReal("problem", "space_density_factor", 1.e-4);
  r_replenish = pin->GetOrAddReal("problem", "r_replenish_Rp", 0.75)*Rp;


  // Tripathi initial conditions variables
  r_0 = 0.5*Rp;
  r_e = 1.02*Rp;

  rho_0 = rho_func(r_0);
  rho_e = rho_func(r_e);
  P_0   = press_func(rho_0);
  P_e   = press_func(rho_e);

  // Domain logic
  if (mesh_size.x1rat == -1.0) {
    meshgen_x1 = MeshGenerator(mesh_size.x1min, mesh_size.x1max, mesh_size.nx1, pin);
    EnrollUserMeshGenerator(X1DIR, MeshSpacingX1);
  }
  if (f2) {
    if (mesh_size.x2rat == -1.0) {
      meshgen_x2 = MeshGenerator(mesh_size.x2min, mesh_size.x2max, mesh_size.nx2, pin);  
      EnrollUserMeshGenerator(X2DIR, MeshSpacingX2);
    }
  }
  if (f3) {
    if (mesh_size.x3rat == -1.0) {
      meshgen_x3 = MeshGenerator(mesh_size.x3min, mesh_size.x3max, mesh_size.nx3, pin);
      EnrollUserMeshGenerator(X3DIR, MeshSpacingX3);
    }
  }
}

Real MeshSpacingX1(Real x, RegionSize rs) {
  return meshgen_x1.MeshSpacing(x);
}

Real MeshSpacingX2(Real x, RegionSize rs) {
  return meshgen_x2.MeshSpacing(x);
}

Real MeshSpacingX3(Real x, RegionSize rs) {
  return meshgen_x3.MeshSpacing(x);
}

// outputs and stored values
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  // User outputs
  AllocateUserOutputVariables(3);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "t2");
  SetUserOutputVariableName(2, "flux");
  // SetUserOutputVariableName(2, "ionization_rate");
  // SetUserOutputVariableName(3, "recombination_rate");
  // SetUserOutputVariableName(4, "recombinative_cooling");
  // SetUserOutputVariableName(5, "lya_cooling");
  // SetUserOutputVariableName(6, "du");

  // User mesh data
  // 0 -- energy/time absorbed into absorber -- turns into ions
  // 1 -- radiation energy absorbed prad->du
  // 2 -- ionization rate
  // 3 -- recombination rate
  // 4 -- recombinative cooling rate
  // 5 -- lyman alpha cooling rate
  // 6 -- output du?
  AllocateRealUserMeshBlockDataField(1);
  ruser_meshblock_data[0].NewAthenaArray(ncells3, ncells2, ncells1);
  // ruser_meshblock_data[1].NewAthenaArray(ncells3, ncells2, ncells1);
  // ruser_meshblock_data[2].NewAthenaArray(ncells3, ncells2, ncells1);
  // ruser_meshblock_data[3].NewAthenaArray(ncells3, ncells2, ncells1);
  // ruser_meshblock_data[4].NewAthenaArray(ncells3, ncells2, ncells1);
  // ruser_meshblock_data[5].NewAthenaArray(ncells3, ncells2, ncells1);
  // ruser_meshblock_data[6].NewAthenaArray(ncells3, ncells2, ncells1);
}

//----------------------------------------------------------------------------------------
// Output Variables
//----------------------------------------------------------------------------------------

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  int il, iu, jl, ju, kl, ku;
  getMBBounds(this, il, iu, jl, ju, kl, ku);

  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        Real R = pthermo->GetRd();
        Real T = phydro->w(IPR,k,j,i)/(R * phydro->w(IDN,k,j,i));
        Real ion_f = pscalars->r(HPLUS,k,j,i);
        T = T * (1.-ion_f/2.);

        Real t2;
        peos->Temperature(phydro->w, pscalars->s, pscalars->m, t2, k, j, i);

        user_out_var(0,k,j,i) = T;
        user_out_var(1,k,j,i) = t2;
        // user_out_var(1,k,j,i) = phydro->hsrc.g1(k,j,i);
        // user_out_var(2,k,j,i) = phydro->hsrc.g2(k,j,i);
        // user_out_var(3,k,j,i) = phydro->hsrc.g3(k,j,i);
        user_out_var(2,k,j,i) = ruser_meshblock_data[0](k,j,i);

        // for (int o=0; o<6; ++o) {
        //   user_out_var(o+4,k,j,i) = ruser_meshblock_data[o+1](k,j,i);
        // }
      }
    }
  }
}

void ReactionNetwork::InitUserReactions(ParameterInput *pin) {
  Reaction *rxn = new H_recombination(this, "H Recombination", HYD, HPLUS, ELEC);
  // AddReaction(rxn);
  return;
}

Absorber* RadiationBand::GetAbsorberByName(std::string name, ParameterInput *pin)
{ 
  std::stringstream msg;
  if (name == "HYDROGEN_IONIZATION") {
    HydrogenIonization *a = new HydrogenIonization(this, HYD);

    ReactionNetwork *pnetwork = pmy_rad->pmy_block->pnetwork;
    std::string rxn_name = "Hydrogen Ionization";
    Reaction *rxn = new Photoionization(pnetwork, rxn_name, a, HPLUS, ELEC, a->ionization_energy);

    return a;
  }
  else {
    msg << "### FATAL ERROR in RadiationBand::AddAbsorber"
        << std::endl << "unknown absorber: '" << name <<"' ";
    ATHENA_ERROR(msg);
  }
}

//----------------------------------------------------------------------------------------
// Additional Physics
//----------------------------------------------------------------------------------------


// sources in ghost cells?
void SourceTerms(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &w, const AthenaArray<Real> &r,
  const AthenaArray<Real> &bcc,
  AthenaArray<Real> &du, AthenaArray<Real> &ds)
{
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  AthenaArray<Real> vol;
  vol.NewAthenaArray(pmb->ncells1);

  PassiveScalars *ps = pmb->pscalars;

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      pmb->pcoord->CellVolume(k,j,is,ie,vol);

      for (int i = is; i <= ie; ++i) {
        // number densities
        Real n_ion = ps->s(HPLUS,k,j,i)/mh;
        Real n_neu = ps->s(HYD,k,j,i)/mh;

        // //ionization
        // // negative sign bc downward energy transfer
        // Real energy_ion = -dt * pmb->ruser_meshblock_data[0](k,j,i);
        // Real n_ion_gain = energy_ion/Ry/vol(i);

        // // Real n = ps->s(0,k,j,i)/mh;
        // // Real Fx = pmb->prad->pband->bflxdn(k,j,i)/2.563e-18;// Divide by 16eV per photon to get photon number flux
        // // Real I = Fx * (n_neu * A0);
        // // Real n_ion_gain = dt * I;

        // if (n_ion_gain > n_neu) {
        //   std::stringstream msg;
        //   msg << "### FATAL ERROR in SourceTerms" << std::endl
        //       << "    n_ion_gain ("<< n_ion_gain <<") is greater than n_neutral ("<< n_neu <<")." << std::endl
        //       << "    Re-run with lower timestep." << std::endl;
        // }

        Real n_ion_gain = 0.0;

        // recombination
        Real R = pmb->pthermo->GetRd();
        Real T = w(IPR,k,j,i)/(R * w(IDN,k,j,i));
        Real ion_f = ps->r(HPLUS,k,j,i);
        T = T * (1.- ion_f/2.);

        // Real T = pmb->pthermo->Temp(pmb->phydro->w.at(k,j,i));
        Real alpha_B  = 2.59E-19 * pow(T/1.E4,-0.7); // m3 s-1
        Real n_recomb = dt * alpha_B * (n_ion * n_ion); //m-3

        // cooling processes
        // -- recomb cooling
        Real kb = 1.3806504E-23; // J/K
        Real recomb_cooling_const = 6.11E-16 * pow(T,-0.89); // m3 s-1
        Real recomb_cooling_rate = recomb_cooling_const * (kb * T) * (n_ion * n_ion); //J s-1 m-3

        // -- lya cooling
        Real lya_cooling_const = 7.5E-32 * exp(-118348./T); // J m3 s-1
        Real lya_cooling_rate = lya_cooling_const * n_ion * n_neu; // J m-3 s-1

        // du(IEN,k,j,i) -= (recomb_cooling_rate + lya_cooling_rate) * dt; // J m-3
        du(IEN,k,j,i) -= (lya_cooling_rate) * dt; // J m-3
        
        n_recomb = 0;
        ds(HYD,k,j,i) += ( n_recomb - n_ion_gain) * mh;
        ds(HPLUS,k,j,i) += (-n_recomb + n_ion_gain) * mh;
        ds(ELEC,k,j,i) += (-n_recomb + n_ion_gain) * ps->m(0);

        pmb->ruser_meshblock_data[0](k,j,i) = pmb->prad->my_bands(0)->spectral_flux_density(0,k,j,i);

        // // Outputs
        // // 1 -- radiation energy absorbed prad->du
        // pmb->ruser_meshblock_data[1](k,j,i) = pmb->prad->dE(k,j,i);
        // // 2 -- ionization rate
        // // pmb->ruser_meshblock_data[2](k,j,i) = n_ion_gain/dt;
        // // 3 -- recombination rate
        // // pmb->ruser_meshblock_data[3](k,j,i) = n_recomb/dt;
        // // 4 -- recombinative cooling rate
        // pmb->ruser_meshblock_data[4](k,j,i) = recomb_cooling_rate;
        // // 5 -- lyman alpha cooling rate
        // pmb->ruser_meshblock_data[5](k,j,i) = lya_cooling_rate;
        // // 6 -- output du?
        // pmb->ruser_meshblock_data[6](k,j,i) = du(IEN,k,j,i);
      }
    }
  }

  // clear eng deposited for next iteration
  // pmb->ruser_meshblock_data[0].ZeroClear();
}

// needs validation
void gravity_func(MeshBlock *pmb, AthenaArray<Real> &g1, AthenaArray<Real> &g2, AthenaArray<Real> &g3) {
  Real x,y,z;
  Real r, rsq, rs, rs_sq;
  Real gp, gs, gc;

  bool spherical_coords;
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    spherical_coords = false;
  }
  else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    spherical_coords = true;
  }

  int il, iu, jl, ju, kl, ku;
  getMBBounds(pmb, il, iu, jl, ju, kl, ku);

  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        // simplify in 1D spherical
        if (spherical_coords) {
          x = pmb->pcoord->x1v(i);
          y = z = 0;

          r = x;
          rs = a-x;
        }
        else {
          x = pmb->pcoord->x1v(i);
          y = pmb->pcoord->x2v(j);
          z = pmb->pcoord->x3v(k);

          // planet placed at (x,y,z) = (0,0,0)
          rsq = x*x + y*y + z*z;
          r = sqrt(rsq);

          // star placed at (x,y,z) = (a,0,0)
          rs_sq = (a-x)*(a-x) + y*y + z*z;
          rs = sqrt(rs_sq);
          // rs = (a-x, -y, -z)
        }

        // enforce rmin for center of planet
        // so long as less than reset radius, this shouldn't matter
        if (r < 0.5 * Rp) {
          r = 0.5*Rp;
        }

        // magnitudes (div by r)
        // planet -- in (-r) direction
        gp = (G * Mp) / pow(r,3);

        // star -- in (a-r) direction
        gs = (G * Ms) / pow(rs,3);
        // centrifugal -- in ( -(a-r) ) direction (note negative sign)
        gc = (G * Ms) / pow(a,3);

        // missing negative signs?
        g1(k, j, i) = gp * (-x) + (gs - gc) * (a-x);
        g2(k, j, i) = gp * (-y) + (gs - gc) * (-y);
        g3(k, j, i) = gp * (-z) + (gs - gc) * (-z);
      }
    }
  }
}

Real RadiationTime(AthenaArray<Real> const &prim, Real time, int k, int j) {
  Real rad_scaling = (ref_dist_*ref_dist_)/(dist_*dist_);

  // tripathi
  Real time_factor = 5 * erf(time/8.e4 - 1.5)+5.1;
  // Real time_factor = (std::tanh(time/2000. - 5.)+1.)/2.;

  return rad_scaling*time_factor;
}

//----------------------------------------------------------------------------------------
// Initial Conditions
//----------------------------------------------------------------------------------------
void SetInitialConditions(Real rad, Real &dens, Real &press, Real &ion_f, Real &v1, Real &v2, Real &v3) {
  v1 = v2 = v3 = 0;
  
  if (rad <= r_0) {
    dens  = rho_0;
    press = P_0;
    ion_f = sfloor;
  }
  else if (rad <= r_e) {
    dens  = rho_func(rad);
    press = press_func(dens);
    ion_f = sfloor;
  }
  else {
    dens  = rho_e * space_density_factor;
    press = P_e;
    ion_f = 1-sfloor;
  }
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  if (NSCALARS != 3) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Problem Generator" << std::endl
        << "    NSCALARS ("<< NSCALARS <<") must be exactly 3." << std::endl;
    ATHENA_ERROR(msg);
  }

  // setup initial condition
  int il, iu, jl, ju, kl, ku;
  getMBBounds(this, il, iu, jl, ju, kl, ku);

  Real rad;
  Real dens, press, ion_f, v1, v2, v3;

  v1 = v2 = v3 = 0;

  for (int i = il; i <= iu; ++i) {
    for (int j = jl; j <= ju; ++j) {
      for (int k = kl; k <= ku; ++k) {
        rad = getRad(pcoord, i, j, k);

        SetInitialConditions(rad, dens, press, ion_f, v1, v2, v3);

        phydro->w(IPR,k,j,i) = press;
        phydro->w(IDN,k,j,i) = dens;
        phydro->w(IV1,k,j,i) = v1;
        phydro->w(IV2,k,j,i) = v2;
        phydro->w(IV3,k,j,i) = v3;

        // s0 -- neutral hydrogen
        pscalars->s(HYD,k,j,i) = (1-ion_f) * dens;
        // s1 -- ionized hydrogen
        pscalars->s(HPLUS,k,j,i) = ion_f * dens;
        pscalars->s(ELEC, k,j,i) = ion_f * dens;
      }
    }
  }
  //-- end file loading

  // set spectral properties
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      prad->CalculateRadiativeTransfer(phydro->w, pscalars->s, pmy_mesh->time, k, j);

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);

  pnetwork->Initialize();
}

// resets things inside r_replenish
void MeshBlock::UserWorkInLoop() {
  // happens at last stage, at end

  Real gm1 = gas_gamma - 1.0;
  Real rad, dens, press, ion_f, v1, v2, v3;

  bool nancheck = false;
  std::stringstream nan_msg;
  nan_msg << "### FATAL ERROR" << std::endl;
  nan_msg << "    nan value detected in (";

  // calculate energy that goes into ionizing
  int il, iu, jl, ju, kl, ku;
  getMBBounds(this, il, iu, jl, ju, kl, ku);

  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        rad = getRad(pcoord, i, j, k);

        if (std::isnan(phydro->w(IDN,k,j,i))) {
          nancheck = true;
          nan_msg << "dens";
        } else if (std::isnan(phydro->w(IPR,k,j,i))) {
          nancheck = true;
          nan_msg << "press";
        } else if (std::isnan(phydro->w(IV1,k,j,i))) {
          nancheck = true;
          nan_msg << "vel1";
        } else if (std::isnan(phydro->w(IV2,k,j,i))) {
          nancheck = true;
          nan_msg << "vel2";
        } else if (std::isnan(phydro->w(IV3,k,j,i))) {
          nancheck = true;
          nan_msg << "vel3";
        }


        if (nancheck) {
          nan_msg << ") at (k=" << k << ", j=" << j << ", i=" << i << ")." << std::endl;
          nan_msg << " at time: " << pmy_mesh->time << " cycle: " << pmy_mesh->ncycle << std::endl;
          std::cout << nan_msg.str() << std::endl;
          ATHENA_ERROR(nan_msg);
        }

        // reset conditions interior to r_replenish
        // probably best to do both cons and prim
        if (rad <= r_replenish) {
          SetInitialConditions(rad, dens, press, ion_f, v1, v2, v3);

          // prim
          phydro->w(IPR,k,j,i) = press;
          phydro->w(IDN,k,j,i) = dens;
          phydro->w(IV1,k,j,i) = v1;
          phydro->w(IV2,k,j,i) = v2;
          phydro->w(IV3,k,j,i) = v3;

          //cons
          phydro->u(IEN,k,j,i) = press/gm1;
          phydro->u(IDN,k,j,i) = dens;
          phydro->u(IM1,k,j,i) = v1*dens;
          phydro->u(IM2,k,j,i) = v2*dens;
          phydro->u(IM3,k,j,i) = v3*dens;

          // scalars
          // s0 -- neutral hydrogen
          pscalars->s(HYD,k,j,i) = (1.0-ion_f) * dens;
          pscalars->r(HYD,k,j,i) = (1.0-ion_f);
          // s1 -- ionized hydrogen
          pscalars->s(HPLUS,k,j,i) = ion_f * dens;
          pscalars->r(HPLUS,k,j,i) = ion_f;

          // s1 -- ionized hydrogen
          pscalars->s(ELEC,k,j,i) = ion_f * dens;
          pscalars->r(ELEC,k,j,i) = ion_f;
        }
      } // i
    } // j
  } // k
  return;
}


//----------------------------------------------------------------------------------------
// Utility
//----------------------------------------------------------------------------------------

Real rho_func(Real r) {  
  Real gm1   = gas_gamma - 1.0;

  Real K = pow(rho_p,1.0-gas_gamma) * pow(cs,2);
  Real rho_term_1 = gm1/gas_gamma * G * Mp/K;
  Real rho_term_2 = (1.0/r - 1.0/Rp);
  Real rho_term_3 = pow(rho_p,gm1);

  return pow(rho_term_1*rho_term_2 + rho_term_3, 1.0/gm1);
}

Real press_func(Real rho) {
  Real K = pow(rho_p,1.0-gas_gamma) * pow(cs,2);
  return K * pow(rho, gas_gamma);
}

Real getRad(Coordinates *pcoord, int i, int j, int k) {
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    Real x = pcoord->x1v(i);
    Real y = pcoord->x2v(j);
    Real z = pcoord->x3v(k);

    return  sqrt(x*x + y*y + z*z);
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    return pcoord->x1v(i);
  }
  else {
    std::stringstream msg;
    msg << "### FATAL ERROR" << std::endl
        << "    Coordinate System (" << COORDINATE_SYSTEM << ") must be either cartesian or spherical_polar." << std::endl;
    ATHENA_ERROR(msg);
  }
}

void getMBBounds(MeshBlock *pmb, int &il, int &iu, int &jl, int &ju, int &kl, int &ku) {
  il = pmb->is-NGHOST;
  iu = pmb->ie+NGHOST;

  jl = pmb->block_size.nx2 == 1 ? pmb->js : pmb->js-NGHOST;
  ju = pmb->block_size.nx2 == 1 ? pmb->je : pmb->je+NGHOST;
  
  kl = pmb->block_size.nx3 == 1 ? pmb->ks : pmb->ks-NGHOST;
  ku = pmb->block_size.nx3 == 1 ? pmb->ke : pmb->ke+NGHOST;
}