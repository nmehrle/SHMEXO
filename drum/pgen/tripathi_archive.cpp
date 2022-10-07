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
#include "../radiation/freedman_mean.hpp"
#include "../radiation/hydrogen_absorber.hpp"
#include "../radiation/user_absorber.hpp"
#include "../radiation/null_absorber.hpp"
#include "../scalars/scalars.hpp"

// global parameters
Real T_iso, rho_ambient, p_ambient;
//gravity
Real G, Mp, Ms, Rp, period, a;
Real wave_to_meters_conversion;

//x is ion fraction (1-rho_n/rho_tot)
// stores initial conditions
AthenaArray<Real> rho0,press0,x0;
Real gas_gamma;

Real dist_, ref_dist_;

AthenaArray<Real> rad_flux_ion;
AthenaArray<Real> rad_eng;
AthenaArray<Real> ionization_rate, recombination_rate, recombinative_cooling, lya_cooling, output_du;

//constants
Real A0=6.30431812E-22; // (m^2)
Real nu_0=3.2898419603E15; // (1/s) Rydberg frequency
Real c=2.998E8; // m/s
Real mh=1.674E-27; // kg
Real nm_0=91.126705058; // Hydrogen ionization wavelength in nm
Real Ry=2.1798723611E-18;// J (joules)

//----------------------------------------------------------------------------------------
// Output Variables
//----------------------------------------------------------------------------------------
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(10);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "g1");
  SetUserOutputVariableName(2, "rad_eng");
  SetUserOutputVariableName(3, "ionization_rate");
  SetUserOutputVariableName(4, "recombination_rate");
  SetUserOutputVariableName(5, "recombinative_cooling");
  SetUserOutputVariableName(6, "lya_cooling");
  SetUserOutputVariableName(7, "du");
  SetUserOutputVariableName(8, "vol");
  SetUserOutputVariableName(9, "farea");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  Real dq[1+NVAPOR], rh;

  AthenaArray<Real> vol, farea;
  vol.NewAthenaArray(ncells1);
  farea.NewAthenaArray(ncells1+1);

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      pcoord->Face1Area(k, j, is, ie+1, farea);
      pcoord->CellVolume(k,j, is, ie,   vol);

      for (int i = is-NGHOST; i <= ie+NGHOST; ++i) {
        Real R = pthermo->GetRd();
        Real T = phydro->w(IPR,k,j,i)/(R * phydro->w(IDN,k,j,i));
        Real ion_f = pscalars->r(1,k,j,i);
        T = T * (1.-ion_f/2.);

        user_out_var(0,k,j,i) = T;
        user_out_var(1,k,j,i) = phydro->hsrc.g1(k,j,i);

        user_out_var(2,k,j,i) = rad_eng(k,j,i);
        user_out_var(3,k,j,i) = ionization_rate(k,j,i);
        user_out_var(4,k,j,i) = recombination_rate(k,j,i);
        user_out_var(5,k,j,i) = recombinative_cooling(k,j,i);
        user_out_var(6,k,j,i) = lya_cooling(k,j,i);
        user_out_var(7,k,j,i) = output_du(k,j,i);

        user_out_var(8,k,j,i) = vol(i);
        user_out_var(9,k,j,i) = farea(i);
      }
    }
  }
}
//----------------------------------------------------------------------------------------
// Absorber Info
//----------------------------------------------------------------------------------------
Real AbsorptionCoefficient(Absorber const *pabs, Real wave, Real const prim[], int k, int j, int i)
{
  // convert microns to meters
  Real freq = c / (wave * wave_to_meters_conversion);

  PassiveScalars *ps = pabs->pmy_band->pmy_rad->pmy_block->pscalars;
  Real n_neutral = ps->s(0,k,j,i)/mh;

  // cover low energy and edge case
  if (freq < nu_0)
    return 0.;
  else if (freq == nu_0)
    return A0*n_neutral;

  return A0 * n_neutral; // 1/m
}

Real EnergyAbsorption(Absorber *pabs, Real wave, Real flux, int k, int j, int i)
{
  Real wave_nm = (wave * wave_to_meters_conversion) * 1e9;
  if (wave_nm > nm_0) {
    // energy not absorbed
    // should also be reflected in tau_lambda = 0
    return flux;
  }
  else {
    // fraction of energy turned into heat
    // removes ionizatoin energy cost
    Real energy_fraction = 0.15;

    rad_flux_ion(k,j,i) += (1-energy_fraction) * flux;

    return flux*energy_fraction;
  }
}

void RadiationBand::AddAbsorber(std::string name, std::string file, ParameterInput *pin)
{ 
  std::stringstream msg;
  if (name == "FREEDMAN") {
    pabs->AddAbsorber(FreedmanMean(this));
  } else if (name == "NULL") {
    pabs->AddAbsorber(NullAbsorber(this));
  } else if (name == "HYDROGEN") {
    pabs->AddAbsorber(HydrogenAbsorber(this, pin));
  } else if (name == "HYDROGEN_IONIZATION") {
    UserDefinedAbsorber hydrogen_ionization(this);

    hydrogen_ionization.EnrollUserAbsorptionCoefficientFunc(AbsorptionCoefficient);
    hydrogen_ionization.EnrollUserEnergyAbsorptionFunc(EnergyAbsorption);
    pabs->AddAbsorber(hydrogen_ionization);
  }else {
    msg << "### FATAL ERROR in RadiationBand::AddAbsorber"
        << std::endl << "unknown absorber: '" << name <<"' ";
    ATHENA_ERROR(msg);
  }
}

//----------------------------------------------------------------------------------------
// Additional Physics
//----------------------------------------------------------------------------------------
void MeshBlock::UserWorkInLoop() {
  // happens at last stage, at end

  Real gm1 = gas_gamma - 1.0;
  Real x,y,z, rad, rad_sq, dens, press, ion_f;

  // calculate energy that goes into ionizing
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        x = pcoord->x1v(i);
        y = pcoord->x2v(j);
        z = pcoord->x3v(k);

        rad_sq = x*x + y*y + z*z;
        rad = sqrt(rad_sq);

        // reset conditions interior to 0.75 Rp
        // probably best to do both cons and prim
        if (rad <= 0.75 * Rp) {
          dens  = rho0(k,j,i);
          press = press0(k,j,i);
          ion_f = x0(k,j,i);

          // prim
          phydro->w(IPR,k,j,i) = press;
          phydro->w(IDN,k,j,i) = dens;
          phydro->w(IV1,k,j,i) = 0;
          phydro->w(IV2,k,j,i) = 0;
          phydro->w(IV3,k,j,i) = 0;

          //cons
          phydro->u(IEN,k,j,i) = press/gm1;
          phydro->u(IDN,k,j,i) = dens;
          phydro->u(IM1,k,j,i) = 0;
          phydro->u(IM2,k,j,i) = 0;
          phydro->u(IM3,k,j,i) = 0;

          // scalars
          // s0 -- neutral hydrogen
          pscalars->s(0,k,j,i) = (1.0-ion_f) * dens;
          pscalars->r(0,k,j,i) = (1.0-ion_f);
          // s1 -- ionized hydrogen
          pscalars->s(1,k,j,i) = ion_f * dens;
          pscalars->r(1,k,j,i) = ion_f;
        }
      } // i
    } // j
  } // k
  return;
}

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
        // just for output land
        rad_eng(k,j,i) = pmb->prad->du(k,j,i); // J m-3 s-1

        // number densities
        Real n_ion = ps->s(1,k,j,i)/mh;
        Real n_neu = ps->s(0,k,j,i)/mh;

        //ionization
        // negative sign bc downward energy transfer
        // Real energy_ion = -dt * rad_flux_ion(k,j,i);
        // Real n_ion_gain = energy_ion/Ry/vol(i);

        // Real n = ps->s(0,k,j,i)/mh;
        Real Fx = pmb->prad->pband->bflxdn(k,j,i)/2.563e-18;// Divide by 16eV per photon to get photon number flux
        Real I = Fx * (n_neu * A0);
        Real n_ion_gain = dt * I;

        if (n_ion_gain > n_neu) {
          std::stringstream msg;
          msg << "### FATAL ERROR in SourceTerms" << std::endl
              << "    n_ion_gain ("<< n_ion_gain <<") is greater than n_neutral ("<< n_neu <<")." << std::endl
              << "    Re-run with lower timestep." << std::endl;
        }

        // recombination
        Real R = pmb->pthermo->GetRd();
        Real T = w(IPR,k,j,i)/(R * w(IDN,k,j,i));
        Real ion_f = ps->r(1,k,j,i);
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

        du(IEN,k,j,i) -= (recomb_cooling_rate + lya_cooling_rate) * dt; // J m-3
        ds(0,k,j,i) += ( n_recomb - n_ion_gain) * mh;
        ds(1,k,j,i) += (-n_recomb + n_ion_gain) * mh;

        ionization_rate(k,j,i) = n_ion_gain/dt;
        recombination_rate(k,j,i) = n_recomb/dt;

        lya_cooling(k,j,i) = lya_cooling_rate;
        recombinative_cooling(k,j,i) = recomb_cooling_rate;
        output_du(k,j,i) = du(IEN,k,j,i);
      }
    }
  }

  // clear eng deposited for next iteration
  rad_flux_ion.ZeroClear();
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

        // planet
        gc1 = - (G * Mp) / (rad_sq);
        // star
        gc2 = (G * Ms) / pow((a-rad),2);
        // centrifugal
        gc3 = - (G * Ms) * (a-rad) / pow(a,3);

        g1(k, j, i) = gc1+gc2+gc3;
      }
    }
  }
}

Real RadiationTime(AthenaArray<Real> const &prim, Real time, int k, int j, int il, int iu) {
  Real rad_scaling = (ref_dist_*ref_dist_)/(dist_*dist_);

  // tripathi
  Real time_factor = 5 * erf(time/8.e4 - 1.5)+5.1;
  // Real time_factor = (std::tanh(time/2000. - 5.)+1.)/2.;

  return rad_scaling*time_factor;
}

//----------------------------------------------------------------------------------------
// User Setup
//----------------------------------------------------------------------------------------
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserExplicitGravityFunction(Gravity);
  EnrollUserExplicitSourceFunction(SourceTerms);


  // Radiation parameters
  dist_ = pin->GetOrAddReal("radiation", "distance", 1.);
  ref_dist_ = pin->GetOrAddReal("radiation", "reference_distance", 1.);

  wave_to_meters_conversion = pin->GetOrAddReal("radiation","wave_to_meters",1.e-9);
  EnrollUserRadiationTimeFunc(RadiationTime);

  // Gravity/System Parameters
  G = pin->GetReal("problem","G");
  Mp = pin->GetReal("problem","Mp");
  Ms = pin->GetReal("problem","Ms");
  Rp = pin->GetReal("problem","Rp");
  period = pin->GetReal("problem","period");

  Real x = 4. * pow(M_PI,2.) / (G * Ms);
  a = pow( pow(period*86400.,2.)/x ,(1./3));

  gas_gamma = pin->GetReal("hydro","gamma");

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

//----------------------------------------------------------------------------------------
// Initial Conditions
//----------------------------------------------------------------------------------------
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  rad_flux_ion.NewAthenaArray(ncells3, ncells2, ncells1);
  rad_flux_ion.ZeroClear();

  rad_eng.NewAthenaArray(ncells3,ncells2,ncells1);
  ionization_rate.NewAthenaArray(ncells3,ncells2,ncells1);
  recombination_rate.NewAthenaArray(ncells3,ncells2,ncells1);
  recombinative_cooling.NewAthenaArray(ncells3,ncells2,ncells1);
  lya_cooling.NewAthenaArray(ncells3,ncells2,ncells1);
  output_du.NewAthenaArray(ncells3,ncells2,ncells1);

  rho0.NewAthenaArray(ncells3,ncells2,ncells1);
  press0.NewAthenaArray(ncells3,ncells2,ncells1);
  x0.NewAthenaArray(ncells3,ncells2,ncells1);

  int nx1 = pmy_mesh->mesh_size.nx1;


  //-- file_loading
  // checks for stuff I haven't included yet
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


  if (NSCALARS != 2) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Problem Generator" << std::endl
        << "    NSCALARS ("<< NSCALARS <<") must be exactly 2." << std::endl;
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

  Real *NEUF = new Real [numRows-1];
  Real *IONF = new Real [numRows-1];

  int i = 0;
  while (std::getline(inp, line)) {
    Real dn_i, pr_i, v1_i, v2_i, v3_i, neu, ion;
    sscanf(line.c_str(), "%le %le %le %le %le %le %le\n", &dn_i, &pr_i, &v1_i, &v2_i, &v3_i, &neu, &ion);
    V1[i] = v1_i;
    V2[i] = v2_i;
    V3[i] = v3_i;

    DN[i] = dn_i;
    PR[i] = pr_i;

    NEUF[i] = neu;
    IONF[i] = ion;

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

        Real press = PR[ii] + p_ambient;
        Real dens  = DN[ii] + rho_ambient;

        phydro->w(IPR,k,j,i) = press;
        phydro->w(IDN,k,j,i) = dens;
        phydro->w(IV1,k,j,i) = V1[ii];
        phydro->w(IV2,k,j,i) = V2[ii];
        phydro->w(IV3,k,j,i) = V3[ii];

        // s0 -- neutral hydrogen
        pscalars->s(0,k,j,i) = NEUF[ii] * dens;

        // s1 -- ionized hydrogen
        pscalars->s(1,k,j,i) = IONF[ii] * dens;


        rho0(k,j,i)   = dens;
        press0(k,j,i) = press;
        x0(k,j,i)     = IONF[ii];
      }
    }
  }
  //-- end file loading

  // set spectral properties
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      prad->CalculateRadiativeTransfer(phydro->w, pmy_mesh->time, k, j, is, ie+1);

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);

  delete[] V1;
  delete[] V2;
  delete[] V3;
  delete[] PR;
  delete[] DN;
}