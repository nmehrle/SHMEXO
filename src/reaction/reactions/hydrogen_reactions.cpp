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

void H_recombination::react(AthenaArray<Real> &dn_rate, AthenaArray<Real> &de_rate, int k, int j, int i) {

  PassiveScalars *ps = pmy_network->pscalars;

  Real T = pmy_network->temperature_(k,j,i);

  Real n_ion = ps->s(ion_scalar_num, k, j, i) / ps->m(ion_scalar_num);
  Real n_elec = ps->s(electron_scalar_num, k, j, i)/ ps->m(electron_scalar_num);

  Real alpha_B  = 2.59E-19 * pow(T/1.E4,-0.7); // m3 s-1
  Real n_recomb = alpha_B * n_ion * n_elec;

  Real recomb_cooling_const = -6.11E-16 * pow(T,-0.89); // m3 s-1
  Real recomb_cooling_rate = recomb_cooling_const * (pmy_network->boltzmann * T) * (n_ion * n_elec); //J s-1 m-3

  dn_rate(scalar_num, k, j, i) += n_recomb;
  dn_rate(ion_scalar_num, k, j, i) -= n_recomb;
  dn_rate(electron_scalar_num, k, j, i) -= n_recomb;

  de_rate(my_rxn_num, k, j, i) += recomb_cooling_rate;
}

void Lya_cooling::react(AthenaArray<Real> &dn_rate, AthenaArray<Real> &de_rate, int k, int j, int i) {

  PassiveScalars *ps = pmy_network->pscalars;

  Real T = pmy_network->temperature_(k,j,i);

  Real n_ion = ps->s(ion_scalar_num, k, j, i) / ps->m(ion_scalar_num);
  Real n_neu = ps->s(scalar_num, k, j, i) / ps->m(scalar_num);

  Real lya_cooling_const = -7.5E-32 * exp(-118348./T); // J m3 s-1
  Real lya_cooling_rate = lya_cooling_const * n_ion * n_neu; // J m-3 s-1

  de_rate(my_rxn_num, k, j, i) += lya_cooling_rate;
}