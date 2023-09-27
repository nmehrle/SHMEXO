// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "../../mesh/mesh.hpp"
#include "../../scalars/scalars.hpp"
#include "../../radiation/radiation.hpp"
#include "../../radiation/absorber/absorber.hpp"
#include "../../radiation/absorber/ionizing_absorber.hpp"
#include "../reaction_network.hpp"
#include "../reaction.hpp"
#include "photoionization.hpp"

Photoionization::Photoionization(std::string name, IonizingAbsorber *pabs,
  int elec_num):
    Reaction(name)
{
  pmy_abs = pabs;
  scalar_num_ = pabs->scalar_num;
  ion_num_ = pabs->ion_num;
  electron_num_ = elec_num;

  ionization_energy = pmy_abs->ionization_energy;
}

void Photoionization::react(int k, int j, int i) {
  PassiveScalars *ps = pmy_network->pscalars;
  Real reaction_energy, contribution;

  Real dn = 0;
  Real d_energy = 0;
  Spectrum *spec = pmy_abs->pmy_band->spec;
  for (int n = 0; n < pmy_abs->pmy_band->nspec; ++n)
  {
    reaction_energy = pmy_abs->energyAbsorbed(n, k, j, i);
    contribution = reaction_energy * pmy_abs->h(n)/ionization_energy;

    dn += contribution;
    d_energy += reaction_energy * pmy_abs->q(n);
  }

  pmy_network->dn_rate(my_rxn_num, scalar_num_, k, j, i) -= dn;
  pmy_network->dn_rate(my_rxn_num, ion_num_, k, j, i) += dn;
  pmy_network->dn_rate(my_rxn_num, electron_num_, k, j, i) += dn;

  pmy_network->de_rate(my_rxn_num, k, j, i) += d_energy;
}