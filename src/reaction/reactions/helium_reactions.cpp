// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "../../mesh/mesh.hpp"
#include "../../scalars/scalars.hpp"
#include "../../radiation/radiation.hpp"
#include "../../radiation/absorber/absorber.hpp"
#include "../reaction_network.hpp"
#include "helium_reactions.hpp"

He_recombination::He_recombination(ReactionNetwork *pnetwork, std::string name,
  int neu_num, int ion_num, int elec_num):
    Reaction(pnetwork, name)
{
  scalar_num = neu_num;
  ion_scalar_num = ion_num;
  electron_scalar_num = elec_num;
  my_name = name;
}

void He_recombination::react(int k, int j, int i) {
  PassiveScalars *ps = pmy_network->pscalars;
  Real T = pmy_network->temperature_(k,j,i);

  Real n_ion = ps->s(ion_scalar_num, k, j, i) / ps->mass(ion_scalar_num);
  Real n_elec = ps->s(electron_scalar_num, k, j, i)/ ps->mass(electron_scalar_num);

  Real alpha_1 = 1e-11 / (0.00217556 * pow(T, 1.12036526) + 0.32053997 * pow(T, 0.61526097));
  Real beta_1 = 1e-11 / (0.00218954 * pow(T, 1.18645797) + 0.352778 * pow(T, 0.62645323));

  Real n_recomb = n_elec * n_ion * alpha_1;
  Real e_recomb = n_elec * n_ion * (pmy_network->boltzmann * T) * beta_1;

  int species[3] = {scalar_num, ion_scalar_num, electron_scalar_num};
  int sign[3]    = {+1, -1, -1};
  for (int l = 0; l < 3; ++l)
  {
    pmy_network->dn_rate(species[l], k, j, i) += sign[l] * n_recomb;
    pmy_network->jacobian(species[l], ion_scalar_num, k, j, i)  += sign[l] * alpha_1 * n_elec;
    pmy_network->jacobian(species[l], electron_scalar_num, k, j, i) += sign[l] * alpha_1 * n_ion;
  }

  pmy_network->de_rate(my_rxn_num, k, j, i) -= e_recomb;
}

He_23S_recombination::He_23S_recombination(ReactionNetwork *pnetwork, std::string name,
  int neu_num, int ion_num, int elec_num):
    Reaction(pnetwork, name)
{
  scalar_num = neu_num;
  ion_scalar_num = ion_num;
  electron_scalar_num = elec_num;
  my_name = name;
}

void He_23S_recombination::react(int k, int j, int i) {
  PassiveScalars *ps = pmy_network->pscalars;
  Real T = pmy_network->temperature_(k,j,i);

  Real n_ion = ps->s(ion_scalar_num, k, j, i) / ps->mass(ion_scalar_num);
  Real n_elec = ps->s(electron_scalar_num, k, j, i)/ ps->mass(electron_scalar_num);

  Real alpha_3 = 1e-11 / (0.00258173 * pow(T, 0.9848205) + 0.10883234 * pow(T, 0.58864659));
  Real beta_3 = 1e-11 / (0.00427277 * pow(T, 0.99204123) + 0.12332369 * pow(T, 0.57871911));

  Real n_recomb = n_elec * n_ion * alpha_3;
  Real e_recomb = n_elec * n_ion * (pmy_network->boltzmann * T) * beta_3;

  int species[3] = {scalar_num, ion_scalar_num, electron_scalar_num};
  int sign[3]    = {+1, -1, -1};
  for (int l = 0; l < 3; ++l)
  {
    pmy_network->dn_rate(species[l], k, j, i) += sign[l] * n_recomb;
    pmy_network->jacobian(species[l], ion_scalar_num, k, j, i)  += sign[l] * alpha_3 * n_elec;
    pmy_network->jacobian(species[l], electron_scalar_num, k, j, i) += sign[l] * alpha_3 * n_ion;
  }

  pmy_network->de_rate(my_rxn_num, k, j, i) -= e_recomb;
}