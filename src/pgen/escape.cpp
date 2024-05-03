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
#include "../reaction/reactions/collisional_ionization.hpp"
#include "../reaction/reactions/recombination.hpp"
#include "../reaction/reactions/helium_recombination.hpp"
#include "../reaction/reactions/collisions.hpp"

#include "../reemission/reemission.hpp"

void FixedBoundaryInnerX1(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
                   FaceField &bb, Real time, Real dt,
                   int il, int iu, int jl, int ju, int kl, int ku, int ngh);

namespace {
  // user set variables
  Real G, Mp, Ms, Rp, period, semi_major_axis;
  Real planck, speed_of_light;
  Real dfloor, pfloor, sfloor;
  Real gas_gamma, boltzmann;

  // radiation variables
  Real global_rad_scaling;
  bool radiation_time_ramp;

  std::string helium_alpha_file, helium_beta_file, helium_decay_file;
  std::string helium_collisions_file, helium_collision_reemission_branching_file;

  bool enable_reemission;

  // initial conditions
  std::string initial_condition_file;
  Real rho_0, Teq;
  Real r_replenish, r_reaction;
  AthenaArray<Real> initial_conditions;

  // species variables
  Real H_He_ratio, H_He_mass_ratio;
  enum speciesList {ELEC = 0,
                H = 1, HII = 2,
                He = 3, He23S = 4,
                HeII = 5, HeIII = 6};

  //convergence checking
  bool check_convergence;
  Real convergence_rel_tol;
  AthenaArray<Real> convergence_abs_tol;
  AthenaArray<Real> convergence_check_buffer;

} //namespace

  // Function declarations
  Real getRad(Coordinates *pcoord, int i, int j, int k);
  void getMBBounds(MeshBlock *pmb, int &il, int &iu, int &jl, int &ju, int &kl, int &ku);
  void gravity_func(MeshBlock *pmb, AthenaArray<Real> &g1, AthenaArray<Real> &g2, AthenaArray<Real> &g3);
  
  void NaNCheck(MeshBlock *pmb, const AthenaArray<Real> prim, int k, int j, int i);
  bool ConvergenceCheck(const AthenaArray<Real> prim, const AthenaArray<Real> prim_scalar, const AthenaArray<Real> scalar_mass, int k, int j, int i);
  void QualifyHHeRatio(PassiveScalars *ps);

  // void SourceTerms(MeshBlock *pmb, const Real time, const Real dt,
  //   const AthenaArray<Real> &w, const AthenaArray<Real> &r,
  //   const AthenaArray<Real> &bcc,
  //   AthenaArray<Real> &du, AthenaArray<Real> &ds);
  Real RadiationScaling(RadiationBand *band, AthenaArray<Real> const &prim, Real time, int k, int j);

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
  // EnrollUserExplicitSourceFunction(SourceTerms);

  dfloor = pin->GetOrAddReal("hydro", "dfloor", 0);
  pfloor = pin->GetOrAddReal("hydro", "pfloor", 0);
  sfloor = pin->GetOrAddReal("hydro", "sfloor", 0);

  check_convergence = pin->GetOrAddBoolean("problem", "check_convergence", false);
  convergence_rel_tol = pin->GetOrAddReal("problem", "convergence_rel_tol", 1.e-6);
  convergence_abs_tol.NewAthenaArray(NHYDRO+NSCALARS);

  convergence_abs_tol(IDN) = pin->GetOrAddReal("problem", "density_abs_tol", -1.);
  convergence_abs_tol(IVX) = pin->GetOrAddReal("problem", "velocity_abs_tol", -1.);
  convergence_abs_tol(IVY) = convergence_abs_tol(IVX);
  convergence_abs_tol(IVZ) = convergence_abs_tol(IVX);
  convergence_abs_tol(IPR) = pin->GetOrAddReal("problem", "pressure_abs_tol", -1.);

  char key[100];
  for (int n = 0; n < NSCALARS; ++n)
  {
    sprintf(key, "s%d.mass_concentration_abs_tol", n);
    convergence_abs_tol(NHYDRO+n) = pin->GetOrAddReal("scalar", key, -1.);
  }

  // Physical Constants
  G = pin->GetReal("problem","G");
  planck = pin->GetReal("problem","planck_constant");
  speed_of_light = pin->GetReal("problem", "speed_of_light");
  boltzmann = pin->GetReal("hydro", "boltzmann");

  // System Parameters
  Mp = pin->GetReal("problem","Mp");
  Ms = pin->GetReal("problem","Ms");
  Rp = pin->GetReal("problem","Rp");
  period = pin->GetReal("problem","period");
  // Convert period to semi-major-axis
  Real x = 4. * pow(M_PI,2.) / (G * Ms);
  semi_major_axis = pow( pow(period*86400.,2.)/x ,(1./3));

  // Initial Conditions
  gas_gamma = pin->GetReal("hydro","gamma");
  Teq = pin->GetReal("problem", "Teq");
  rho_0 = pin->GetReal("problem","rho_0");
  initial_condition_file = pin->GetString("problem","initial_conditions_file");

  H_He_ratio      = pin->GetOrAddReal("problem", "H_He_ratio", 0);
  H_He_mass_ratio = pin->GetOrAddReal("problem", "H_He_mass_ratio", 0);

  // Radiation parameters
  global_rad_scaling = pin->GetOrAddReal("radiation", "global_radiation_scaling", 1.);
  radiation_time_ramp = pin->GetOrAddBoolean("radiation", "time_ramp", true);

  EnrollUserRadiationScalingFunction(RadiationScaling);

  enable_reemission = pin->GetOrAddBoolean("reaction", "enable_reemission", false);

  helium_alpha_file = pin->GetOrAddString("reaction", "He_recombination_alpha_file", "");
  helium_beta_file = pin->GetOrAddString("reaction", "He_recombination_beta_file", "");
  helium_decay_file = pin->GetOrAddString("reaction", "He_recombination_decay_file", "");
  helium_collisions_file = pin->GetOrAddString("reaction", "He_e_collisions_file", "");
  helium_collision_reemission_branching_file = pin->GetOrAddString("reaction", "He_e_coll_reemission_branching_file", "");

  // Domain logic
  r_replenish = pin->GetOrAddReal("problem", "r_replenish_Rp", 0.0)*Rp;
  r_reaction = pin->GetOrAddReal("problem", "r_reaction_Rp", r_replenish/Rp)*Rp;

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
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, FixedBoundaryInnerX1);  
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
    ReactionTemplate *r = new HydrogenicRecombination(name, {H, HII, ELEC}, {+1, -1, -1}, 1, "1S");
    if (enable_reemission) {
      GroundStateReemission *h_reemission = new GroundStateReemission(pin, prad, pmy_block->pscalars->energy(HII));
      r->AssignReemission(h_reemission);
    }
    return r;

  } else if (name == "LYA_COOLING") {
    return new LyaCooling(name, H, ELEC);

  } else if (name == "HE_1S_RECOMBINATION") {
    ReactionTemplate *r = new HeliumRecombination(name, {He, HeII, ELEC}, {+1, -1, -1}, helium_alpha_file, helium_beta_file, 1);
    if (enable_reemission) {
      GroundStateReemission *he_reemission = new GroundStateReemission(pin, prad, pmy_block->pscalars->energy(HeII));
      r->AssignReemission(he_reemission);
    }
    return r;

  } else if (name == "HE_SINGLET_RECOMBINATION") {
    ReactionTemplate *r = new HeliumRecombination(name, {He, HeII, ELEC}, {+1, -1, -1}, helium_alpha_file, helium_beta_file, 2);
    if (enable_reemission) {
      HeliumReemisison *he_reemission = new HeliumReemisison(pin, prad, helium_decay_file, 1.,0.,0.);
      r->AssignReemission(he_reemission);
    }
    return r;

  } else if (name == "HE_TRIPLET_RECOMBINATION") {
    ReactionTemplate *r = new HeliumRecombination(name, {He, HeII, ELEC}, {+1, -1, -1}, helium_alpha_file, helium_beta_file, 3);
    return r;

  } else if (name == "HeElecCollision11") {
    ReactionTemplate *r = new HeElectronCollisions(name, {He, ELEC, He, ELEC}, {-1, -1, +1, +1}, helium_collisions_file, 1, 5);
    if (enable_reemission) {
      HeliumReemisison *he_reemission = new HeliumReemisison(pin, prad, helium_decay_file, 0,0,0);
      he_reemission->AssignBranchingFile(helium_collision_reemission_branching_file, 1, 2);
      r->AssignReemission(he_reemission);
    }

    return r;

  } else if (name == "HeElecCollision13") {
    return new HeElectronCollisions(name, {He, ELEC, He23S, ELEC}, {-1, -1, +1, +1}, helium_collisions_file, 2, 6);

  } else if (name == "HeElecCollision31") {
    ReactionTemplate *r = new HeElectronCollisions(name, {He23S, ELEC, He, ELEC}, {-1, -1, +1, +1}, helium_collisions_file, 3, 7);
    if (enable_reemission) {
      HeliumReemisison *he_reemission = new HeliumReemisison(pin, prad, helium_decay_file, 0,0,0);
      he_reemission->AssignBranchingFile(helium_collision_reemission_branching_file, 3, 4);
      r->AssignReemission(he_reemission);
    }

    return r;

  } else if (name == "HeElecCollision33") {
    return new HeElectronCollisions(name, {He23S, ELEC, He23S, ELEC}, {-1, -1, +1, +1}, helium_collisions_file, 4, 8);

  } else if (name == "HE23S_RADIO_DECAY") {
    ReactionTemplate *r = new HeliumTripletDecay(name, {He23S, He}, {-1, +1});
    if (enable_reemission) {
      GroundStateReemission *he_reemission = new GroundStateReemission(pin, prad, pmy_block->pscalars->energy(He23S));
      r->AssignReemission(he_reemission);
    }
    return r;

  } else if (name == "CHARGE_EXCHANGE_HE_H") {
    return new ChargeExchangeHeliumToHyd(name, {HeII, H, He, HII}, {-1, -1, +1, +1});

  } else if (name == "CHARGE_EXCHANGE_H_HE") {
    return new ChargeExchangeHydToHelium(name,  {He, HII, HeII, H}, {-1, -1, +1, +1});

  } else if (name == "TRIPLET_COLLISIONAL_RELAXATION") {
    return new TripletHydrogenCollision(name, {He23S, H, He, HII, ELEC}, {-1, -1, +1, +1, +1});

  } else if (name == "HEIII_RECOMBINATION") {
    // ReactionTemplate *r = new HydrogenicRecombination(name, {HeII, HeIII, ELEC}, {+1, -1, -1}, 2, "A");
    ReactionTemplate *r = new BadnellRecombination(name, {HeII, HeIII, ELEC}, {+1, -1, -1}, 1.8180e-10, 0.74920, 1.0170e+01, 2.7860e+06);
    if (enable_reemission) {
      Real energy_difference = pmy_block->pscalars->energy(HeIII) - pmy_block->pscalars->energy(HeII);
      GroundStateReemission *he_reemission = new GroundStateReemission(pin, prad, energy_difference);
      r->AssignReemission(he_reemission);
    }
    return r;

  } else if (name == "BREMSSTRAHLUNG") {
    return new Bremsstrahlung("Bremsstrahlung", ELEC, HII, HeII, HeIII);

  } else if (name == "HCI") {
    return new CollisionalIonization(name, {H, ELEC, HII, ELEC}, {-1, -1, +1, +2}, 0);

  } else if (name == "HeCI") {
    return new CollisionalIonization(name, {He, ELEC, HeII, ELEC}, {-1, -1, +1, +2}, 1);

  } else if (name == "HeIICI") {
    return new CollisionalIonization(name, {HeII, ELEC, HeIII, ELEC}, {-1, -1, +1, +2}, 2);

  } else if (name == "He23SCI") {
    return new CollisionalIonization(name, {He23S, ELEC, HeII, ELEC}, {-1, -1, +1, +2}, 3);

  } else if (name == "HeII_elec_collisions") {
    return new HeIIElecCollisions(name, HeII, ELEC);

  } else if (name == "HeI_elec_collisions") {
    return new HeIElecCollisions(name, He, ELEC);

  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in ReactionNetwork::GetReactionByName"
        << std::endl << "unknown reaction: '" << name <<"' ";
    ATHENA_ERROR(msg);
  }
}

Absorber* RadiationBand::GetAbsorberByName(std::string name, std::string band_name, ParameterInput *pin)
{
  ReactionNetwork *prad_network = pmy_rad->pmy_block->prad_network;

  std::string rxn_name = band_name + "_" + name;

  if (name == "HYDROGEN_IONIZATION") {
    HydrogenIonization *a = new HydrogenIonization(this, name, H, HII, pin);

    prad_network->AddReaction(new Photoionization(rxn_name, a, ELEC));
    return a;
  }
  else if (name == "HELIUM_IONIZATION") {
    std::string xc_file = pin->GetString("radiation", "Helium_1(1)S_file");
    HeliumIonization *a = new HeliumIonization(this, name, He, HeII, pin, xc_file);

    prad_network->AddReaction(new Photoionization(rxn_name, a, ELEC));
    return a;
  }
  else if (name == "HELIUM_TRIPLET_IONIZATION") {
    std::string xc_file = pin->GetString("radiation", "Helium_2(3)S_file");
    HeliumIonization *a = new HeliumIonization(this, name, He23S, HeII, pin, xc_file);

    prad_network->AddReaction(new Photoionization(rxn_name, a, ELEC));
    return a; 
  }
  else if (name == "HELIUM_DOUBLE_IONIZATION") {
    // Real Eth, Real Emax, Real E0, Real sigma0, Real ya, Real P, Real yw, Real y0, Real y1
    VernerIonizationParams *vi_params = new VernerIonizationParams(54.42, 50000, 1.72, 13690, 32.88, 2.963, 0., 0., 0.);
    VernerIonization *a = new VernerIonization(this, name, HeII, HeIII, pin, vi_params);

    prad_network->AddReaction(new Photoionization(rxn_name, a, ELEC));
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

Real RadiationScaling(RadiationBand *band, AthenaArray<Real> const &prim, Real time, int k, int j) {
  Real time_factor = 1;
  if (radiation_time_ramp)
    time_factor = (std::tanh((time-8.e4)/5e4)+1.)/2.;

  Real scaling = global_rad_scaling * band->band_scaling_factor;
  return scaling*time_factor;
}

//----------------------------------------------------------------------------------------
// Initial Conditions
//----------------------------------------------------------------------------------------

void QualifyHHeRatio(PassiveScalars *ps) {
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
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  initial_conditions.NewAthenaArray(NHYDRO+NSCALARS, ncells3, ncells2, ncells1);
  convergence_check_buffer.NewAthenaArray(NHYDRO+NSCALARS, ncells3, ncells2, ncells1);

  // setup initial condition
  int il, iu, jl, ju, kl, ku;
  getMBBounds(this, il, iu, jl, ju, kl, ku);

  // Read in data from file
  AthenaArray<Real> tmp_file_data;
  ReadDataTable(tmp_file_data, initial_condition_file);

  const int n_file = tmp_file_data.GetDim2();
  const int n_var = tmp_file_data.GetDim1();
  Real **file_data = new Real*[n_var];
  for(int i = 0; i < n_var; i++)
    file_data[i] = new Real[n_file];

  for (int i = 0; i < n_var; ++i)
  {
    for (int j = 0; j < n_file; ++j)
    {
      file_data[i][j] = tmp_file_data(j,i);
    }
  }
  tmp_file_data.DeleteAthenaArray();

  QualifyHHeRatio(pscalars);
  Real rad, rad_rp, var, number_density;

  for (int i = il; i <= iu; ++i) {
    for (int j = jl; j <= ju; ++j) {
      for (int k = kl; k <= ku; ++k) {
        rad = getRad(pcoord, i, j, k);
        rad_rp = rad/Rp;
        // interp file_data onto domain, scale appropriately
        // save to initial conditions

        number_density = 0;
        for (int n = 0; n < NHYDRO+NSCALARS; ++n)
        {
          //interp file_data onto domain
          var = interp1(rad_rp, file_data[n+1], file_data[0], n_file);

          //scale appropriately
          if (n == IDN) // density
            var *= rho_0;
          else if (n == IPR) // pressure
            var *= Teq;
          else if (n == NHYDRO + 1 || n == NHYDRO + 2)// Hydrogen Species
            var *= H_He_mass_ratio;
          else if (n > NHYDRO + 2) // Helium Species
            var *= (1-H_He_mass_ratio);

          // track number density of scalars for T->P conversion
          if (n >= NHYDRO)
            number_density += var / pscalars->mass(n-NHYDRO);

          //save to initial conditions
          initial_conditions(n,k,j,i) = var;
        } // n

        // convert temp to press and resave
        number_density = number_density * initial_conditions(IDN,k,j,i);
        initial_conditions(IPR,k,j,i) = initial_conditions(IPR,k,j,i) * number_density * boltzmann;

        // write to hydro variables
        for (int n = 0; n < NHYDRO; ++n)
          phydro->w(n,k,j,i) = initial_conditions(n,k,j,i);

        // write to scalar variables
        for (int n = 0; n < NSCALARS; ++n) {
          pscalars->r(n,k,j,i) = initial_conditions(n+NHYDRO,k,j,i);
          pscalars->s(n,k,j,i) = initial_conditions(n+NHYDRO,k,j,i) * initial_conditions(IDN,k,j,i);
        }
      } //k
    } //j
  } //i
  //-- end file loading
  for(int i = 0; i < n_var; i++)
      delete[] file_data[i];
  delete[] file_data;

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);

  for (int b = 0; b < prad->nbands; ++b)
  {
    prad->my_bands(b)->InitAbsorbers(pin);
  }

  prad_network->Initialize();
  pnetwork->Initialize();

  // set spectral properties
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      prad->CalculateRadiativeTransfer(phydro->w, pscalars->s, pmy_mesh->time, k, j);

  prad_network->ComputeReactionForcing(phydro->w, phydro->u, pscalars->s);
}

void NaNCheck(MeshBlock *pmb, const AthenaArray<Real> prim, int k, int j, int i) {
  bool tripped = false;
  std::string trip_var;

  // nan checking
  if (std::isnan(prim(IDN,k,j,i))) {
    tripped = true;
    trip_var = "dens";
  } else if (std::isnan(prim(IPR,k,j,i))) {
    tripped = true;
    trip_var = "press";
  } else if (std::isnan(prim(IV1,k,j,i))) {
    tripped = true;
    trip_var = "vel1";
  } else if (std::isnan(prim(IV2,k,j,i))) {
    tripped = true;
    trip_var = "vel2";
  } else if (std::isnan(prim(IV3,k,j,i))) {
    tripped = true;
    trip_var = "vel3";
  }

  if (tripped) {
    std::stringstream nan_msg;
    nan_msg << "### FATAL ERROR" << std::endl;
    nan_msg << "    nan value detected in (" << trip_var;
    nan_msg << ") at (k=" << k << ", j=" << j << ", i=" << i << ")." << std::endl;
    nan_msg << " at time: " << pmb->pmy_mesh->time << " cycle: " << pmb->pmy_mesh->ncycle << std::endl;
    std::cout << nan_msg.str() << std::endl;
    ATHENA_ERROR(nan_msg);
  }
}

bool ConvergenceCheck(const AthenaArray<Real> prim, const AthenaArray<Real> prim_scalar, const AthenaArray<Real> scalar_mass, int k, int j, int i) {
  bool converged_here = true;
  Real past_value, current_value;
  Real deviation, abs_tol, tol;

  for (int n = 0; n < NHYDRO+NSCALARS; ++n)
  {
    past_value = convergence_check_buffer(n, k, j, i);
    abs_tol = convergence_abs_tol(n);

    if (n < NHYDRO) { // Hydro Variables
      current_value = prim(n, k, j, i);
    }
    else { // Scalar Variables
      current_value = prim_scalar(n-NHYDRO, k, j, i);
    }

    deviation = abs(current_value - past_value);
    tol = std::max(convergence_rel_tol * abs(past_value), abs_tol);

    if (deviation > tol)
      converged_here = false;

    convergence_check_buffer(n,k,j,i) = current_value;
  }

  return converged_here;
}

// resets things inside r_replenish
void MeshBlock::UserWorkInLoop() {
  // happens at last stage, at end

  Real gm1 = gas_gamma - 1.0;
  Real rad, dens;

  bool converged = true;

  int il, iu, jl, ju, kl, ku;
  getMBBounds(this, il, iu, jl, ju, kl, ku);

  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        rad = getRad(pcoord, i, j, k);

        NaNCheck(this, phydro->w, k, j, i);

        // convergence check
        if (check_convergence && rad > r_replenish) {
          if (not ConvergenceCheck(phydro->w, pscalars->r, pscalars->mass, k, j, i))
            converged = false;
        }

        // reset conditions interior to r_replenish
        // probably best to do both cons and prim
        if (rad <= r_replenish) {
          // write to hydro variables
          for (int n = 0; n < NHYDRO; ++n) {
            phydro->w(n,k,j,i) = initial_conditions(n,k,j,i);
          }

          // write to scalar variables
          for (int n = 0; n < NSCALARS; ++n) {
            pscalars->r(n,k,j,i) = initial_conditions(n+NHYDRO,k,j,i);
            pscalars->s(n,k,j,i) = initial_conditions(n+NHYDRO,k,j,i) * initial_conditions(IDN,k,j,i);
          }

          //cons
          dens = phydro->w(IDN,k,j,i);
          phydro->u(IEN,k,j,i) = phydro->w(IPR,k,j,i)/gm1;
          phydro->u(IDN,k,j,i) = dens;
          phydro->u(IM1,k,j,i) = phydro->w(IV1,k,j,i)*dens;
          phydro->u(IM2,k,j,i) = phydro->w(IV2,k,j,i)*dens;
          phydro->u(IM3,k,j,i) = phydro->w(IV3,k,j,i)*dens;
        }
      } // i
    } // j
  } // k

  if (check_convergence && converged) {
    pmy_mesh->EndTimeIntegration();
  }
  return;
}

//----------------------------------------------------------------------------------------
// Utility
//----------------------------------------------------------------------------------------

bool ReactionNetwork::DoRunReactions(int k, int j, int i) {
  Real r = getRad(pmy_block->pcoord, i, j, k);
  if (r < r_reaction){
    return false;
  }
  return true;
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


void FixedBoundaryInnerX1(MeshBlock *pmb, Coordinates *pcoord, AthenaArray<Real> &prim,
                   FaceField &bb, Real time, Real dt,
                   int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for (int n = 0; n < NHYDRO; ++n) {
          prim(n,k,j,il-i) = initial_conditions(n,k,j,il-i);
        }
      }
    }
  }
}
