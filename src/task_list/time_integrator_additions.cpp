// C/C++ headers
#include <iostream>
#include <sstream>

// Athena++ headers
#include "../athena.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../chemistry/chemistry.hpp"
#include "../radiation/radiation.hpp"
#include "../reaction/reaction_network.hpp"
#include "../scalars/scalars.hpp"
#include "task_list.hpp"

//----------------------------------------------------------------------------------------
// Functions for implicit correction
TaskStatus TimeIntegratorTaskList::UpdateHydro(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  Real dt = pmb->pmy_mesh->dt;

  if (stage <= nstages) {
    if (ph->implicit_flag == 1)
      ph->ImplicitCorrectionFull(ph->du, ph->w, stage_wghts[stage-1].beta*dt);
    else if (ph->implicit_flag == 2)
      ph->ImplicitCorrectionReduced(ph->du, ph->w, stage_wghts[stage-1].beta*dt);
    Real wghts[3];
    wghts[0] = 1.;
    wghts[1] = 1.;
    wghts[2] = 0.;
    pmb->WeightedAve(ph->u, ph->du, ph->u2, wghts);
  }

  return TaskStatus::next;
}

TaskStatus TimeIntegratorTaskList::UpdateScalars(MeshBlock *pmb, int stage) {
  PassiveScalars *ps = pmb->pscalars;
  Hydro *ph = pmb->phydro;

  if (stage <= nstages) {
    if (ph->implicit_flag == 1)
      ps->ImplicitCorrectionFull(ps->ds, ph->implicit_correction);
    else if (ph->implicit_flag == 2)
      ps->ImplicitCorrectionReduced(ps->ds, ph->implicit_correction);
    Real wghts[3];
    wghts[0] = 1.;
    wghts[1] = 1.;
    wghts[2] = 0.;
    pmb->WeightedAve(ps->s, ps->ds, ps->s1, wghts);
  }

  return TaskStatus::next;
}

//----------------------------------------------------------------------------------------
// Functions to integrate chemistry
TaskStatus TimeIntegratorTaskList::IntegrateChemistry(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;

  if (stage != nstages) return TaskStatus::next;

  // frictional heating
  // pmb->pchem->AddFrictionalHeating(ph->u, ph->w, pmb->pmy_mesh->dt);

  // do slow chemistry
  // pmb->pchem->EvolveOneStep(ph->u, pmb->pmy_mesh->time, pmb->pmy_mesh->dt);

  // do fast chemistry
  // pmb->pthermo->SaturationAdjustment(ph->u);

  return TaskStatus::next;
}
//----------------------------------------------------------------------------------------

// Functions to calculate radiation flux
TaskStatus TimeIntegratorTaskList::CalculateRadiationFlux(MeshBlock *pmb, int stage)
{
  #if RADIATION_ENABLED
    if (stage != nstages) return TaskStatus::next;

    Radiation *prad = pmb->prad;
    PassiveScalars *pscalars = pmb->pscalars;
    Hydro *phydro=pmb->phydro;

    int js = pmb->js; int ks = pmb->ks;
    int je = pmb->je; int ke = pmb->ke;

    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        Real t_start_stage = pmb->pmy_mesh->time;
        prad->CalculateRadiativeTransfer(phydro->w, pscalars->s, t_start_stage, k, j);
      }
    }
  #endif
  return TaskStatus::next;
}

TaskStatus TimeIntegratorTaskList::IntegrateRadiation(MeshBlock *pmb, int stage)
{
  #if RADIATION_ENABLED
    if (stage != nstages) return TaskStatus::next;

    ReactionNetwork *prad_network = pmb->prad_network;
    PassiveScalars *pscalars = pmb->pscalars;
    Hydro *phydro = pmb->phydro;

    Real dt = pmb->pmy_mesh->dt;
    prad_network->ComputeReactionForcing(phydro->w, phydro->u, pscalars->s);
    prad_network->IntegrateReactions(dt, pscalars->s, phydro->u, pscalars->s);
  #endif
  return TaskStatus::next;
}

TaskStatus TimeIntegratorTaskList::IntegrateReactions(MeshBlock *pmb, int stage)
{
  #if REACTION_ENABLED
    if (stage != nstages) return TaskStatus::next;

    ReactionNetwork *pnetwork = pmb->pnetwork;
    PassiveScalars *pscalars = pmb->pscalars;
    Hydro *phydro = pmb->phydro;

    Real dt = pmb->pmy_mesh->dt;
    pnetwork->ComputeReactionForcing(phydro->w, phydro->u, pscalars->s);
    pnetwork->IntegrateReactions(dt, pscalars->s, phydro->u, pscalars->s);
  #endif
  return TaskStatus::next;
}