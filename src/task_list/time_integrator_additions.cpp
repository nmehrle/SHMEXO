// C/C++ headers
#include <iostream>

// Athena++ headers
#include "../athena.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../chemistry/chemistry.hpp"
#include "../radiation/radiation.hpp"
#include "../scalars/scalars.hpp"
#include "task_list.hpp"

//----------------------------------------------------------------------------------------
// Functions for implicit correction
enum TaskStatus TimeIntegratorTaskList::UpdateHydro(MeshBlock *pmb, int stage) {
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

enum TaskStatus TimeIntegratorTaskList::UpdateScalars(MeshBlock *pmb, int stage) {
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
  pmb->pchem->AddFrictionalHeating(ph->u, ph->w, pmb->pmy_mesh->dt);

  // do slow chemistry
  pmb->pchem->EvolveOneStep(ph->u, pmb->pmy_mesh->time, pmb->pmy_mesh->dt);

  // do fast chemistry
  pmb->pthermo->SaturationAdjustment(ph->u);

  return TaskStatus::next;
}

//----------------------------------------------------------------------------------------
// Functions to calculate radiation flux
TaskStatus TimeIntegratorTaskList::CalculateRadiationFlux(MeshBlock *pmb, int stage) {
  // only do radiation at first rk step  -- quickening assumption for now
  if (stage != 1) return TaskStatus::next;

  Radiation *prad = pmb->prad;
  Hydro *phydro=pmb->phydro;

  if (prad->current > 0.) {
    prad->current -= pmb->pmy_mesh->dt;  // radiation is in cool-down
  } else {
    prad->current = prad->cooldown;

    prad->ClearRadFlux();

    int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
    int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
    int jl, ju, kl, ku;

    jl = js, ju = je, kl = ks, ku = ke;

    if (pmb->block_size.nx2 > 1) {
      if (pmb->block_size.nx3 == 1) // 2D
        jl = js-1, ju = je+1, kl = ks, ku = ke;
      else // 3D
        jl = js-1, ju = je+1, kl = ks-1, ku = ke+1;
    }

    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        prad->CalculateFluxes(phydro->w, pmb->pmy_mesh->time, k, j, is, ie+1);
  }
  return TaskStatus::next;
}

// after integrate hydro
// before update hydro
TaskStatus TimeIntegratorTaskList::AddSourceTermsRadiation(MeshBlock *pmb, int stage) {
  Radiation *prad = pmb->prad;
  Hydro *ph       = pmb->phydro;

  if (stage <= nstages) {
    Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);
    prad->AddRadiationSourceTerm(dt, ph->du);
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::next;
}
