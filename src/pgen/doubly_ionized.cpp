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

#include "../radiation/absorber/hydrogen_ionization.hpp"
#include "../radiation/absorber/helium_ionization.hpp"
#include "../radiation/absorber/verner_ionization.hpp"


#include "../reaction/reaction_network.hpp"
#include "../reaction/reaction.hpp"

#include "../reaction/reactions/photoionization.hpp"
#include "../reaction/reactions/additional_reactions.hpp"
#include "../reaction/reactions/recombination.hpp"
#include "../reaction/reactions/collisions.hpp"

void FixedBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
                   FaceField &bb, Real time, Real dt,
                   int il, int iu, int jl, int ju, int kl, int ku, int ngh);

namespace {
  // user set variables
  Real G, Mp, Ms, Rp, period, semi_major_axis;
  Real dfloor, pfloor, sfloor;
  Real gas_gamma, boltzmann;
  Real wave_to_meters_conversion;

  // radiation variables
  Real global_rad_scaling;

  // tripathi conditions variables
  Real rho_0, Teq;
  Real r_e, r_replenish, r_inner;

  // species variables
  Real H_He_ratio, H_He_mass_ratio, mu;
  enum speciesList {ELEC = 0,
                H = 1, HII = 2,
                He = 3, He23S = 4,
                HeII = 5, HeIII = 6};
  Real epsilon_concentration;

  AthenaArray<Real> initial_abundances;

  bool check_convergence;
  Real convergence_level;
  AthenaArray<Real> last_step_cons;
  AthenaArray<Real> last_step_cons_scalar;
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
  Real RadiationTime(RadiationBand *band, AthenaArray<Real> const &prim, Real time, int k, int j);

  // mesh generators
  MeshGenerator meshgen_x1, meshgen_x2, meshgen_x3;
  Real MeshSpacingX1(Real x, RegionSize rs);
  Real MeshSpacingX2(Real x, RegionSize rs);
  Real MeshSpacingX3(Real x, RegionSize rs);
  void SetInitialAbundances(MeshBlock *pmb, PassiveScalars *ps);

//----------------------------------------------------------------------------------------
// User Setup
//----------------------------------------------------------------------------------------
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserExplicitGravityFunction(gravity_func);
  // EnrollUserExplicitSourceFunction(SourceTerms);

  check_convergence = pin->GetOrAddBoolean("problem", "check_convergence", false);
  convergence_level = pin->GetOrAddReal("problem", "convergence_level", 1.e-5);

  // Radiation parameters
  global_rad_scaling = pin->GetOrAddReal("radiation", "global_radiation_scaling", 1.);

  wave_to_meters_conversion = pin->GetOrAddReal("radiation","wave_to_meters",1.e-7);
  EnrollUserRadiationScalingFunction(RadiationTime);

  // Gravity/System Parameters
  G = pin->GetReal("problem","G");
  Mp = pin->GetReal("problem","Mp");
  Ms = pin->GetReal("problem","Ms");
  Rp = pin->GetReal("problem","Rp");
  period = pin->GetReal("problem","period");

  Real x = 4. * pow(M_PI,2.) / (G * Ms);
  semi_major_axis = pow( pow(period*86400.,2.)/x ,(1./3));

  gas_gamma = pin->GetReal("hydro","gamma");
  boltzmann = pin->GetReal("hydro", "boltzmann");

  dfloor = pin->GetOrAddReal("hydro", "dfloor", 0);
  pfloor = pin->GetOrAddReal("hydro", "pfloor", 0);
  sfloor = pin->GetOrAddReal("hydro", "sfloor", 0);

  // rho_p = pin->GetOrAddReal("problem","rho_p",1.0e-15);
  // cs = pin->GetOrAddReal("problem","cs",3.0e5);
  // space_density_factor = pin->GetOrAddReal("problem", "space_density_factor", 1.e-4);
  r_replenish = pin->GetOrAddReal("problem", "r_replenish_Rp", 1.0)*Rp;
  r_e = pin->GetOrAddReal("problem", "r_exterior_Rp", 2.0)*Rp;
  r_inner = pin->GetOrAddReal("problem", "r_inner_Rp", 1.0)*Rp;
  
  Teq = pin->GetReal("problem", "Teq");
  rho_0 = pin->GetReal("problem","rho_0");

  // 92.3% H by number == 75% H by mass
  H_He_ratio      = pin->GetOrAddReal("problem", "H_He_ratio", 0);
  H_He_mass_ratio = pin->GetOrAddReal("problem", "H_He_mass_ratio", 0);
  epsilon_concentration = pin->GetOrAddReal("reaction", "epsilon_concentration", 0);


  // Tripathi initial conditions variables

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

  if (GetBoundaryFlag(pin->GetString("mesh", "ix1_bc")) == BoundaryFlag::user) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, FixedBoundary);  
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
  AllocateUserOutputVariables(1);
  SetUserOutputVariableName(0, "temp");
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

        Real T;
        peos->Temperature(phydro->w, pscalars->s, pscalars->mass, T, k, j, i);
        user_out_var(0,k,j,i) = T;
      }
    }
  }
}

Reaction* ReactionNetwork::GetReactionByName(std::string name, ParameterInput *pin)
{
  if (name == "NULL_REACTION") {
    return new NullReaction();

  } else if (name == "H_RECOMBINATION_CASE_A") {
    return new HydrogenicRecombination(name, {H, HII, ELEC}, {+1, -1, -1}, 1, "A");

  } else if (name == "H_RECOMBINATION_CASE_B") {
    return new HydrogenicRecombination(name, {H, HII, ELEC}, {+1, -1, -1}, 1, "B");

  } else if (name == "H_RECOMBINATION_CASE_1S") {
    AthenaArray<Real> wave_temp, frac_temp;
    wave_temp.NewAthenaArray(1);
    frac_temp.NewAthenaArray(1);
    wave_temp(0) = 91.16;
    frac_temp(0) = 1;
    return new HydrogenicRecombination(name, {H, HII, ELEC}, {+1, -1, -1}, 1, "1S", this, wave_temp, frac_temp);

  } else if (name == "LYA_COOLING") {
    return new LyaCooling(name, H, HII);

  } else if (name == "HE_RECOMBINATION") {
    return new HeliumRecombination(name, {He, HeII, ELEC}, {+1, -1, -1});

  } else if (name == "HE23S_RECOMBINATION") {
    return new HeliumTripletRecombination(name, {He23S, HeII, ELEC}, {+1, -1, -1});

  } else if (name == "HE_E_COLLISIONS") {
    std::string collisions_file_name = pin->GetString("reaction", "HE_COLLISIONS_FILE");

    return new He_e_collisions(collisions_file_name, name, He, He23S, ELEC);

  } else if (name == "HE23S_RADIO_DECAY") {
    return new HeliumTripletDecay(name, {He23S, He}, {-1, +1});

  } else if (name == "CHARGE_EXCHANGE_HE_H") {
    return new ChargeExchangeHeliumToHyd(name, {HeII, H, He, HII}, {-1, -1, +1, +1});

  } else if (name == "CHARGE_EXCHANGE_H_HE") {
    return new ChargeExchangeHydToHelium(name,  {He, HII, HeII, H}, {-1, -1, +1, +1});

  } else if (name == "TRIPLET_COLLISIONAL_RELAXATION") {
    return new TripletHydrogenCollision(name, {He23S, H, He, HII, ELEC}, {-1, -1, +1, +1, +1});

  } else if (name == "HEIII_RECOMBINATION_CASE_A") {
    // a, b, T0, T1
    // VernerRecombinationParams *vr_params = new VernerRecombinationParams(1.891e-10, 0.7524, 9.370, 2.774e6);
    // return new VernerRecombination(name, {HeII, HeIII, ELEC}, {+1, -1, -1}, vr_params);

    return new HydrogenicRecombination(name, {HeII, HeIII, ELEC}, {+1, -1, -1}, 2, "A");

  } else if (name == "HEIII_RECOMBINATION_CASE_B") {
    return new HydrogenicRecombination(name, {HeII, HeIII, ELEC}, {+1, -1, -1}, 2, "B");

  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in ReactionNetwork::GetReactionByName"
        << std::endl << "unknown reaction: '" << name <<"' ";
    ATHENA_ERROR(msg);
  }
}

Absorber* RadiationBand::GetAbsorberByName(std::string name, std::string band_name, ParameterInput *pin)
{
  ReactionNetwork *pnetwork = pmy_rad->pmy_block->pnetwork;

  std::string rxn_name = band_name + "_" + name;

  if (name == "HYDROGEN_IONIZATION") {
    HydrogenIonization *a = new HydrogenIonization(this, name, H, HII, pin);

    pnetwork->AddReaction(new Photoionization(rxn_name, a, ELEC));
    return a;
  }
  else if (name == "HELIUM_IONIZATION") {
    std::string xc_file = pin->GetString("radiation", "Helium_1(1)S_file");
    HeliumIonization *a = new HeliumIonization(this, name, He, HeII, pin, xc_file);

    pnetwork->AddReaction(new Photoionization(rxn_name, a, ELEC));
    return a;
  }
  else if (name == "HELIUM_TRIPLET_IONIZATION") {
    std::string xc_file = pin->GetString("radiation", "Helium_2(3)S_file");
    HeliumIonization *a = new HeliumIonization(this, name, He23S, HeII, pin, xc_file);

    pnetwork->AddReaction(new Photoionization(rxn_name, a, ELEC));
    return a; 
  }
  else if (name == "HELIUM_DOUBLE_IONIZATION") {
    // Real Eth, Real Emax, Real E0, Real sigma0, Real ya, Real P, Real yw, Real y0, Real y1
    VernerIonizationParams *vi_params = new VernerIonizationParams(54.42, 50000, 1.72, 13690, 32.88, 2.963, 0., 0., 0.);
    VernerIonization *a = new VernerIonization(this, name, HeII, HeIII, pin, vi_params);

    pnetwork->AddReaction(new Photoionization(rxn_name, a, ELEC));
    return a;
  }
  else {
    std::stringstream msg;
    msg << "### FATAL ERROR in RadiationBand::AddAbsorber"
        << std::endl << "unknown absorber: '" << name <<"' ";
    ATHENA_ERROR(msg);
  }
}

//----------------------------------------------------------------------------------------
// Additional Physics
//----------------------------------------------------------------------------------------


// void SourceTerms(MeshBlock *pmb, const Real time, const Real dt,
//   const AthenaArray<Real> &w, const AthenaArray<Real> &r,
//   const AthenaArray<Real> &bcc,
//   AthenaArray<Real> &du, AthenaArray<Real> &ds)
// {
//   int is = pmb->is, js = pmb->js, ks = pmb->ks;
//   int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

//   AthenaArray<Real> vol;
//   vol.NewAthenaArray(pmb->ncells1);

//   PassiveScalars *ps = pmb->pscalars;

//   for (int k = ks; k <= ke; ++k) {
//     for (int j = js; j <= je; ++j) {
//       pmb->pcoord->CellVolume(k,j,is,ie,vol);

//       for (int i = is; i <= ie; ++i) {
//       }
//     }
//   }
// }

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
          rs = semi_major_axis-x;
        }
        else {
          x = pmb->pcoord->x1v(i);
          y = pmb->pcoord->x2v(j);
          z = pmb->pcoord->x3v(k);

          // planet placed at (x,y,z) = (0,0,0)
          rsq = x*x + y*y + z*z;
          r = sqrt(rsq);

          // star placed at (x,y,z) = (a,0,0)
          rs_sq = (semi_major_axis-x)*(semi_major_axis-x) + y*y + z*z;
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
        gc = (G * Ms) / pow(semi_major_axis,3);

        // missing negative signs?
        g1(k, j, i) = gp * (-x) + (gs - gc) * (semi_major_axis-x);
        g2(k, j, i) = gp * (-y) + (gs - gc) * (-y);
        g3(k, j, i) = gp * (-z) + (gs - gc) * (-z);
      }
    }
  }
}

Real RadiationTime(RadiationBand *band, AthenaArray<Real> const &prim, Real time, int k, int j) {
  // tripathi
  Real time_factor = 5 * erf(time/8.e4 - 1.5)+5.1;
  // Real time_factor = (std::tanh(time/2000. - 5.)+1.)/2.;
  Real scaling = global_rad_scaling * band->band_scaling_factor;

  return scaling*time_factor;
}

//----------------------------------------------------------------------------------------
// Initial Conditions
//----------------------------------------------------------------------------------------
void SetInitialConditions(Real rad, Real &dens, Real &press, Real &v1, Real &v2, Real &v3) {
  v1 = v2 = v3 = 0;

    if (rad <= r_inner) {
      dens = rho_func(r_inner);
      press = press_func(dens);
    }
    else {
      dens  = rho_func(rad);
      press = press_func(dens);  
    }
}

void SetInitialAbundances(MeshBlock *pmb, PassiveScalars *ps) {
  int il, iu, jl, ju, kl, ku;
  Real rad;
  Real me  = ps->mass(ELEC);
  Real mh  = ps->mass(H);
  Real mhe = ps->mass(He);

  // validate H/He ratio
  bool ratio_error = false;
  if (H_He_ratio <= 0 && H_He_mass_ratio <= 0) {
    ratio_error = true;
  }
  else if (H_He_ratio > 0 && H_He_mass_ratio > 0){
    ratio_error = true;
  }
  if (ratio_error) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Problem File" << std::endl
        << "Must specify exactly one of H_He_ratio or H_He_mass_ratio" << std::endl;
    ATHENA_ERROR(msg);
  }
  if (H_He_ratio > 0) {
    // mass ratio undefined
    Real r = H_He_ratio;
    H_He_mass_ratio = (r * mh) / (r * mh + (1-r)*mhe);
  }
  else {
    // number ratio undefined
    Real q = H_He_mass_ratio;
    H_He_ratio = (q/mh) / (q / mh + (1-q)/mhe);
  }

  Real q = H_He_mass_ratio;
  Real norm;

  getMBBounds(pmb, il, iu, jl, ju, kl, ku);

  // set number densities, fill in mass later
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) { 
        rad = getRad(pmb->pcoord, i, j, k);

        if (rad <= r_e) {
          initial_abundances(H,k,j,i)   = q;
          initial_abundances(ELEC,k,j,i)  = 4*epsilon_concentration;
          initial_abundances(HII,k,j,i) = epsilon_concentration;
          initial_abundances(He,k,j,i) = (1.-q);
          initial_abundances(He23S,k,j,i) = epsilon_concentration;
          initial_abundances(HeII,k,j,i) = epsilon_concentration;
          initial_abundances(HeIII,k,j,i) = epsilon_concentration;
        }
        else {
          initial_abundances(H,k,j,i)   = epsilon_concentration;
          initial_abundances(ELEC,k,j,i)  = 1. + 2*epsilon_concentration;
          initial_abundances(HII,k,j,i) = q;
          initial_abundances(He,k,j,i) = epsilon_concentration;
          initial_abundances(He23S,k,j,i) = epsilon_concentration;
          initial_abundances(HeII,k,j,i) = 1.-q;
          initial_abundances(HeIII,k,j,i) = epsilon_concentration; 
        }

        // add up accross n
        norm = 0;
        for (int n = 0; n < NSCALARS; ++n)
        {
          // initial_abundances(n,k,j,i) *= ps->m(n);
          norm += initial_abundances(n,k,j,i);
        }

        for (int n = 0; n < NSCALARS; ++n)
        {
          initial_abundances(n,k,j,i) /= norm;
        }
      } // i
    } // j
  } // k

  mu = 0;
  for (int n = 0; n < NSCALARS; ++n)
  {
    mu += initial_abundances(n,0,0,0)/pmb->pscalars->mass(n);
  }
  mu = 1.0/mu;
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  initial_abundances.NewAthenaArray(NSCALARS, ncells3, ncells2, ncells1);
  SetInitialAbundances(this, pscalars);

  last_step_cons.NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1);
  last_step_cons_scalar.NewAthenaArray(NSCALARS, ncells3, ncells2, ncells1);

  // setup initial condition
  int il, iu, jl, ju, kl, ku;
  getMBBounds(this, il, iu, jl, ju, kl, ku);

  Real rad;
  Real dens, press, v1, v2, v3;

  v1 = v2 = v3 = 0;

  for (int i = il; i <= iu; ++i) {
    for (int j = jl; j <= ju; ++j) {
      for (int k = kl; k <= ku; ++k) {
        rad = getRad(pcoord, i, j, k);

        SetInitialConditions(rad, dens, press, v1, v2, v3);

        phydro->w(IPR,k,j,i) = press;
        phydro->w(IDN,k,j,i) = dens;
        phydro->w(IV1,k,j,i) = v1;
        phydro->w(IV2,k,j,i) = v2;
        phydro->w(IV3,k,j,i) = v3;

        // s0 -- neutral hydrogen
        for (int n = 0; n < NSCALARS; ++n)
        {
          pscalars->s(n,k,j,i) = initial_abundances(n,k,j,i) * dens;
        }
      }
    }
  }
  //-- end file loading

  for (int b = 0; b < prad->nbands; ++b)
  {
    prad->my_bands(b)->InitAbsorbers(pin);
  }

  pnetwork->Initialize();

  // set spectral properties
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      prad->CalculateRadiativeTransfer(phydro->w, pscalars->s, pmy_mesh->time, k, j);

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}

// resets things inside r_replenish
void MeshBlock::UserWorkInLoop() {
  // happens at last stage, at end

  Real gm1 = gas_gamma - 1.0;
  Real rad, dens, press, v1, v2, v3;

  bool nancheck = false;
  std::stringstream nan_msg;
  nan_msg << "### FATAL ERROR" << std::endl;
  nan_msg << "    nan value detected in (";

  bool converged = true;

  int il, iu, jl, ju, kl, ku;
  getMBBounds(this, il, iu, jl, ju, kl, ku);

  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        rad = getRad(pcoord, i, j, k);

        // nan checking
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

        // convergence check
        if (check_convergence && rad > r_replenish) {
          for (int nh = 0; nh < NHYDRO; ++nh)
          {
            Real current_value = phydro->u(nh, k, j, i);
            Real perc_difference = abs(current_value - last_step_cons(nh, k, j, i))/(current_value);
            if (perc_difference > convergence_level) {
              converged = false;
            }
            last_step_cons(nh, k, j, i) = current_value;
          }

          for (int ns = 0; ns < NSCALARS; ++ns)
          {
            Real current_value = pscalars->s(ns, k, j, i);
            Real perc_difference = abs(current_value - last_step_cons_scalar(ns, k, j, i))/(current_value);
            if (perc_difference > convergence_level) {
              converged = false;
            }
            last_step_cons_scalar(ns, k, j, i) = current_value;
          }
        }

        // reset conditions interior to r_replenish
        // probably best to do both cons and prim
        if (rad <= r_replenish) {
          SetInitialConditions(rad, dens, press, v1, v2, v3);

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
          for (int n = 0; n < NSCALARS; ++n)
          {
            pscalars->r(n,k,j,i) = initial_abundances(n,k,j,i);
            pscalars->s(n,k,j,i) = initial_abundances(n,k,j,i) * dens;
          }
        }
      } // i
    } // j
  } // k
  
  if (converged) {
    pmy_mesh->EndTimeIntegration();
  }
  return;
}

//----------------------------------------------------------------------------------------
// Utility
//----------------------------------------------------------------------------------------

bool ReactionNetwork::DoRunReactions(int k, int j, int i) {
  Real r = getRad(pmy_block->pcoord, i, j, k);
  if (r < r_replenish){
    return false;
  }
  return true;
}

Real rho_func(Real r) {
  Real Hinv = (G * Mp * mu) / (boltzmann * Teq * Rp * Rp);
  Real exp_r = (Rp - r);
  return rho_0 * std::exp(Hinv * exp_r);
}

Real press_func(Real rho) {
  return rho * Teq * boltzmann/mu;
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

void FixedBoundary(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
                   FaceField &bb, Real time, Real dt,
                   int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  return;
}