// C/C++ headers
#include <sstream>

// Athena++ header
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
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

  InitUserReactions(pin);
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
  dn_rate.NewAthenaArray(NSCALARS, pmy_block->ncells3, pmy_block->ncells2, pmy_block->ncells1);

  my_reactions.NewAthenaArray(num_reactions);

  Reaction *p = pfirst;
  for (int r = 0; r < num_reactions; ++r)
  {
    my_reactions(r) = p;
    p=p->next;
  }
}

void ReactionNetwork::ComputeReactionForcing(const Real dt, AthenaArray<Real> &du, AthenaArray<Real> &ds) {
  MeshBlock *pmb = pmy_block;
  PassiveScalars *pscalars = pmb->pscalars;
  dn_rate.ZeroClear();
  de_rate.ZeroClear();

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  // loop over space
  // compute rate values at each point,
  // add to total value
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {

        // populates de_rate(k,j,i), dn_rate(k,j,i)
        // by summing contributions from reacations
        for (int r = 0; r < num_reactions; ++r)
        {
          Reaction *p = my_reactions(r);
          p->react(dn_rate, de_rate, k, j, i);
        }

        du(IEN, k, j, i) += de_rate(k,j,i) * dt;

        for (int n = 0; n < NSCALARS; ++n)
        {
          ds(n,k,j,i) += dn_rate(n,k,j,i) * dt;
        }

      } // i
    } // j
  } // k
}






