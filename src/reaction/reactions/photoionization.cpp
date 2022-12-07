// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "../../mesh/mesh.hpp"
#include "../../scalars/scalars.hpp"
#include "../../radiation/radiation.hpp"
#include "../../radiation/absorber/absorber.hpp"
#include "../reaction_network.hpp"
#include "photoionization.hpp"

Photoionization::Photoionization(ReactionNetwork *pnetwork, std::string name, Absorber *pabs,
  int ion_num, int elec_num, Real ion_eng):
    Reaction(pnetwork, name)
{
  pmy_abs = pabs;
  scalar_num = pabs->scalar_num;
  ion_scalar_num = ion_num;
  ionization_energy = ion_eng;
  electron_scalar_num = elec_num;
  my_name = name;
}

void Photoionization::react(AthenaArray<Real> &dn_rate, AthenaArray<Real> &de_rate, int k, int j, int i) {
  Real reaction_energy, contribution;

  Real dn = 0;
  Real d_energy   = 0;
  Spectrum *spec = pmy_abs->pmy_band->spec;
  for (int n = 0; n < pmy_abs->pmy_band->nspec; ++n)
  {
    reaction_energy = pmy_abs->energyAbsorbed(n, k, j, i);
    contribution = reaction_energy * pmy_abs->h(n) / ionization_energy;

    dn += contribution * spec[n].wgt;
    d_energy   += reaction_energy * pmy_abs->q(n) * spec[n].wgt;
  }

  dn_rate(scalar_num, k, j, i) -= dn;
  dn_rate(ion_scalar_num, k, j, i) += dn;
  dn_rate(electron_scalar_num, k, j, i) += dn;
  de_rate(my_rxn_num, k, j, i) += d_energy;
}