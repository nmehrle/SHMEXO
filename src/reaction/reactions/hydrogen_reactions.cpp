// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "../../mesh/mesh.hpp"
#include "../../scalars/scalars.hpp"
#include "../../radiation/radiation.hpp"
#include "../../radiation/absorber/absorber.hpp"
#include "../reaction_network.hpp"
#include "hydrogen_reactions.hpp"

H_recombination::H_recombination(ReactionNetwork *pnetwork, std::string name,
  int neu_num, int ion_num, int elec_num):
    Reaction(pnetwork, name)
{
  scalar_num = neu_num;
  ion_scalar_num = ion_num;
  electron_scalar_num = elec_num;
  my_name = name;
}

Lya_cooling::Lya_cooling(ReactionNetwork *pnetwork, std::string name,
  int neu_num, int ion_num, int elec_num):
    Reaction(pnetwork, name)
{
  scalar_num = neu_num;
  ion_scalar_num = ion_num;
  electron_scalar_num = elec_num;
  my_name = name;
}

void H_recombination::react(int k, int j, int i) {

  PassiveScalars *ps = pmy_network->pscalars;

  Real T = pmy_network->temperature_(k,j,i);

  Real n_ion = ps->s(ion_scalar_num, k, j, i) / ps->mass(ion_scalar_num);
  Real n_elec = ps->s(electron_scalar_num, k, j, i)/ ps->mass(electron_scalar_num);

  Real alpha_B  = 2.59E-13 * pow(T/1.E4,-0.7); // cm3 s-1
  Real n_recomb = alpha_B * n_ion * n_elec;

  Real beta = 6.11E-10 * pow(T,-0.89); // cm3 s-1
  Real e_recomb = beta * (pmy_network->boltzmann * T) * (n_ion * n_elec); //J s-1 m-3

  int species[3] = {scalar_num, ion_scalar_num, electron_scalar_num};
  int sign[3]    = {+1, -1, -1};
  for (int l = 0; l < 3; ++l)
  {
    pmy_network->dn_rate(species[l], k, j, i) += sign[l] * n_recomb;
    pmy_network->jacobian(species[l], ion_scalar_num, k, j, i)      += sign[l] * alpha_B * n_elec;
    pmy_network->jacobian(species[l], electron_scalar_num, k, j, i) += sign[l] * alpha_B * n_ion;
  }

  pmy_network->de_rate(my_rxn_num, k, j, i) -= e_recomb;
}

void Lya_cooling::react(int k, int j, int i) {

  PassiveScalars *ps = pmy_network->pscalars;

  Real T = pmy_network->temperature_(k,j,i);

  Real n_ion = ps->s(ion_scalar_num, k, j, i) / ps->mass(ion_scalar_num);
  Real n_neu = ps->s(scalar_num, k, j, i) / ps->mass(scalar_num);

  Real beta = -7.5E-19 * exp(-118348./T); // erg cm3 s-1
  Real de_lya = beta * n_ion * n_neu; // J m-3 s-1

  pmy_network->de_rate(my_rxn_num, k, j, i) += de_lya;
}