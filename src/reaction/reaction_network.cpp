// C/C++ headers
#include <sstream>

// Athena++ header
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../utils/utils.hpp"
#include "../radiation/radiation.hpp"
#include "../scalars/scalars.hpp"
#include "reaction_network.hpp"
#include "reaction.hpp"

ReactionNetwork::ReactionNetwork(MeshBlock *pmb, ParameterInput *pin){
  pmy_block = pmb;
  pscalars = pmb->pscalars;
  prad = pmb->prad;

  pfirst = nullptr;
  plast = nullptr;
  num_reactions = 0;

  boltzmann = pin->GetReal("hydro", "boltzmann");

  eV_conversion = pin->GetReal("problem", "eV_conversion");
  ry_conversion = pin->GetReal("problem", "ry_conversion");

  implicit_flag = pin->GetOrAddBoolean("reaction", "implicit_reactions", "False");
  consumption_tolerance = pin->GetOrAddReal("reaction", "consumption_tolerance", 0.0);
  temperature_.NewAthenaArray(pmb->ncells3, pmb->ncells2, pmb->ncells1);
  
  speed_of_light = pin->GetReal("problem", "speed_of_light");
  planck_constant = pin->GetReal("problem", "planck_constant");
}

ReactionNetwork::~ReactionNetwork() {
  Reaction *prxn = pfirst;
  Reaction *pnext;
  while (prxn != nullptr) {
    pnext = prxn->next;
    delete prxn;
    prxn = pnext;
  }
}

void ReactionNetwork::ReadReactionsFromInput(ParameterInput *pin) {
  std::string rxn_name = "";
  std::string null_name = "NO_MORE_REACTIONS";
  int rxn_num = 1;
  char key[100];

  // collects rxn1, rxn2, rxn3 etc. 
  // must be in order with no skips
  // this limitation can be fixed later
  sprintf(key, "rxn%d", rxn_num);
  rxn_name = pin->GetOrAddString("reaction", key, null_name);
  while (rxn_name != null_name) {
    Reaction* prxn = GetReactionByName(rxn_name, pin);
    AddReaction(prxn);
  
    rxn_num +=1;
    sprintf(key, "rxn%d", rxn_num);
    rxn_name = pin->GetOrAddString("reaction", key, null_name);
  }
}

void ReactionNetwork::SetImplicitReactions(bool implicit_flag_) {
  implicit_flag = implicit_flag_;
}

void ReactionNetwork::AddReaction(Reaction *prxn) {
  prxn->pmy_network = this;
  Reaction *p = plast;

  if (p == nullptr) {
    pfirst = prxn;
    prxn->my_rxn_num = 0;
    num_reactions = 1;
  }
  else {
    p->next = prxn;
    prxn->prev = p;
    prxn->my_rxn_num = num_reactions;
    ++num_reactions;
  }
  plast = prxn;
}

void ReactionNetwork::Initialize() {
  de_rate.NewAthenaArray(num_reactions, pmy_block->ncells3, pmy_block->ncells2, pmy_block->ncells1);
  dn_rate.NewAthenaArray(num_reactions, NSCALARS, pmy_block->ncells3, pmy_block->ncells2, pmy_block->ncells1);
  jacobian.NewAthenaArray(NSCALARS, NSCALARS, pmy_block->ncells3, pmy_block->ncells2, pmy_block->ncells1);

  my_reactions.NewAthenaArray(num_reactions);

  Reaction *p = pfirst;
  for (int r = 0; r < num_reactions; ++r)
  {
    my_reactions(r) = p;
    p=p->next;
  }
}

bool __attribute__((weak)) ReactionNetwork::DoRunReactions(int k, int j, int i) { 
  return true;
}

Reaction* __attribute__((weak)) ReactionNetwork::GetReactionByName(std::string name, ParameterInput *pin) {
  return nullptr;
}


void ReactionNetwork::ResetRates() {
  dn_rate.ZeroClear();
  de_rate.ZeroClear();
  jacobian.ZeroClear();
  
  for (int i = 0; i < prad->nbands; ++i)
  {
    prad->my_bands(i)->emission_coefficient.ZeroClear();
  }
}

// computes dn, de
void ReactionNetwork::ComputeReactionForcing(const AthenaArray<Real> prim, const AthenaArray<Real> cons, const AthenaArray<Real> cons_scalar) {
  MeshBlock *pmb = pmy_block;
  PassiveScalars *pscalars = pmb->pscalars;
  ResetRates();

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  // loop over space
  // compute rate values at each point,
  // add to total value
  Real temp;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {

        if (DoRunReactions(k,j,i)) {
          pmb->peos->Temperature(prim, cons_scalar, pscalars->mass, temp, k, j, i);
          temperature_(k,j,i) = temp;

          // populates de_rate(k,j,i), dn_rate(k,j,i)
          // by summing contributions from reactions
          for (int r = 0; r < num_reactions; ++r)
          {
            Reaction *p = my_reactions(r);
            p->react(k, j, i);
          }
        }
      } // i
    } // j
  } // k
}

// Applies dn, de
void ReactionNetwork::IntegrateReactions(const Real dt, const AthenaArray<Real> cons_scalar, AthenaArray<Real> &du, AthenaArray<Real> &ds)
{
  MeshBlock *pmb = pmy_block;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  Real drho[NSCALARS];
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {

        if (DoRunReactions(k,j,i)) {
          for (int r = 0; r < num_reactions; ++r) {
            du(IEN, k, j, i) += de_rate(r, k,j,i) * dt;
          }

          if (implicit_flag)
            IntegrateScalarsImplicit(dt, drho, k, j, i);
          else
            IntegrateScalarsExplicit(dt, drho, k, j, i);

          for (int n = 0; n < NSCALARS; ++n)
          {
            CheckScalarConflict(n,k,j,i, drho, cons_scalar);

            // add to scalar change array
            ds(n,k,j,i) += drho[n];
          }
        }
      } // i
    } // j
  } // k
}

void ReactionNetwork::IntegrateScalarsImplicit(const Real dt, Real drho[NSCALARS], int k, int j, int i) {
  Real dn_rate_total;

  // populate jacobian matrix
  // and RHS vector
  MatrixNSR jacobi_matrix;
  VectorNSR dn_RHS;
  for (int l = 0; l < NSCALARS; ++l)
  {
    for (int m = 0; m < NSCALARS; ++m)
    {
      jacobi_matrix(l,m) = jacobian(l,m,k,j,i);
    }

    dn_rate_total = 0;
    for (int r = 0; r < num_reactions; ++r) {
      dn_rate_total += dn_rate(r,l,k,j,i);
    }
    dn_RHS(l) = dn_rate_total * dt;
  }

  MatrixNSR LHS = MatrixNSR::Identity() - dt * jacobi_matrix;
  VectorNSR dn_sol = LHS.inverse() * dn_RHS;

  for (int n = 0; n < NSCALARS; ++n)
  {
    // convert density rate to density
    drho[n] = dn_sol(n) * pscalars->mass(n);
  }
}

void ReactionNetwork::IntegrateScalarsExplicit(const Real dt, Real drho[NSCALARS], int k, int j, int i) {
  Real dn_rate_total;
  for (int n = 0; n < NSCALARS; ++n)
  {
    dn_rate_total = 0;
    for (int r = 0; r < num_reactions; ++r) {
      dn_rate_total += dn_rate(r,n,k,j,i);
    }
    // convert density rate to density
    drho[n] = dn_rate_total * dt * pscalars->mass(n);
  }
}

void ReactionNetwork::CheckScalarConflict(int n, int k, int j, int i, Real drho[NSCALARS], const AthenaArray<Real> cons_scalar)
{
  Real s_n = cons_scalar(n,k,j,i);
  Real drho_n = drho[n];
  if(drho_n < 0 && -drho_n > s_n) {
    Real deviation = (-drho_n - s_n)/s_n;

    if (std::abs(deviation) < consumption_tolerance) {
      std::stringstream msg;
      msg << " consumption_tolerance activating " << std::endl
          <<"Error at position (k,j,i) = " << k << ", " << j << ", " << i << "." << std::endl
          << "drho ("<<drho_n<<") for scalar (" << n << ") exceeds scalar concentration"
          << "(" << cons_scalar(n,k,j,i)<< ") causing negative density." << std::endl
          << std::endl;
      std::cout << msg.str();
      drho[n] = -s_n;
    } else {
      std::stringstream msg;
      msg << "##### FATAL ERROR in ReactionNetwork::ComputeReactionForcing" << std::endl
          << "Error at position (k,j,i) = " << k << ", " << j << ", " << i << "." << std::endl
          << "drho ("<<drho_n<<") for scalar (" << n << ") exceeds scalar concentration"
          << "(" << cons_scalar(n,k,j,i)<< ") causing negative density." << std::endl
          << "fraction deviation " << deviation << " exceeds consumption tolerance "<<consumption_tolerance
          << std::endl;
      ATHENA_ERROR(msg);
    }
  }
}

Real ReactionNetwork::NewBlockTimeStep(Real dt_ratio) {
  Real real_max = std::numeric_limits<Real>::max();
  Real min_dt_rxn = real_max;
  Real dn_dt, min_dt_n, s_n;

  MeshBlock *pmb = pmy_block;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  for (int n = 0; n < NSCALARS; ++n) {
    min_dt_n = real_max;

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          s_n = pscalars->s(n,k,j,i);
          dn_dt = 0;
          for (int r = 0; r < num_reactions; ++r) {
            dn_dt += dn_rate(r,n,k,j,i);
          }

          if (dn_dt < 0) {
            min_dt_n = std::min(min_dt_n, s_n / (-dn_dt));
          }
        } // i
      } // j
    } // k
    min_dt_n = min_dt_n / pscalars->mass(n);
    min_dt_rxn = std::min(min_dt_rxn, min_dt_n);
  } // n

  min_dt_rxn *= dt_ratio;
  return min_dt_rxn;
}