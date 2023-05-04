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

ReactionNetwork::ReactionNetwork(MeshBlock *pmb, ParameterInput *pin){
  pmy_block = pmb;
  pscalars = pmb->pscalars;

  pfirst = nullptr;
  plast = nullptr;
  num_reactions = 0;

  boltzmann = pin->GetReal("hydro", "boltzmann");

  eV_conversion = pin->GetReal("problem", "eV_conversion");
  ry_conversion = pin->GetReal("problem", "ry_conversion");

  implicit_reactions = pin->GetOrAddBoolean("reaction", "implicit_reactions", "False");
  consumption_tolerance = pin->GetOrAddReal("reaction", "consumption_tolerance", 0.0);
  sfloor = pin->GetOrAddReal("hydro", "sfloor", 1.e-12);
  temperature_.NewAthenaArray(pmb->ncells3, pmb->ncells2, pmb->ncells1);

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
  // InitUserReactions(pin);
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

bool __attribute__((weak)) DoRunReactions(int k, int j, int i) { 
  return true;
}

void ReactionNetwork::ResetRates() {
  dn_rate.ZeroClear();
  de_rate.ZeroClear();
  jacobian.ZeroClear();
}

void ReactionNetwork::ComputeReactionForcing(const Real dt, const AthenaArray<Real> prim, const AthenaArray<Real> cons, const AthenaArray<Real> cons_scalar, AthenaArray<Real> &du, AthenaArray<Real> &ds) {
  MeshBlock *pmb = pmy_block;
  PassiveScalars *pscalars = pmb->pscalars;
  Hydro *phydro = pmb->phydro;
  ResetRates();

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  // loop over space
  // compute rate values at each point,
  // add to total value
  Real temp;
  Real drho[NSCALARS];
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {

        if (DoRunReactions(k,j,i)) {
          pmb->peos->Temperature(prim, cons_scalar, pscalars->mass, temp, k, j, i);
          temperature_(k,j,i) = temp;

          // populates de_rate(k,j,i), dn_rate(k,j,i)
          // by summing contributions from reacations
          for (int r = 0; r < num_reactions; ++r)
          {
            Reaction *p = my_reactions(r);
            p->react(k, j, i);

            du(IEN, k, j, i) += de_rate(r, k,j,i) * dt;
          }

          ComputeScalarDensityChange(dt, drho, k, j, i);
          for (int n = 0; n < NSCALARS; ++n)
          { 

            // check this drho isnt in error
            if (drho[n] < 0 && -drho[n] > cons_scalar(n,k,j,i)) {
              Real s_rho = cons_scalar(n,k,j,i);
              Real deviation = (-drho[n] - s_rho)/s_rho;

              // unnecessary abs
              if (std::abs(deviation) < consumption_tolerance) {
                drho[n] = s_rho * (1.0 - sfloor);
              } else {
                std::stringstream msg;
                msg << "##### FATAL ERROR in ReactionNetwork::ComputeReactionForcing" << std::endl
                    << "Error at position (k,j,i) = " << k << ", " << j << ", " << i << "." << std::endl
                    << "drho ("<<drho[n]<<") for scalar (" << n << ") exceeds scalar concentration"
                    << "(" << cons_scalar(n,k,j,i)<< ") causing negative density." << std::endl
                    << "fraction deviation " << deviation << " exceeds consumption tolerance "<<consumption_tolerance
                    << std::endl;
                ATHENA_ERROR(msg);
              }
            }

            // add to scalar change array
            ds(n,k,j,i) += drho[n];
          }
        }
      } // i
    } // j
  } // k
}

//check for implicit bug here
void ReactionNetwork::ComputeScalarDensityChange(const Real dt, Real drho[NSCALARS], int k, int j, int i) {
  Real dn_rate_total;

  if (implicit_reactions) {
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
  } // if (implicit reactions)
  else {
    for (int n = 0; n < NSCALARS; ++n)
    {
      dn_rate_total = 0;
      for (int r = 0; r < num_reactions; ++r) {
        dn_rate_total += dn_rate(r,n,k,j,i);
      }
      // convert density rate to density
      drho[n] = dn_rate_total * dt * pscalars->mass(n);
    }
  } // else (non-implicit reaction)
}