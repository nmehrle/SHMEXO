// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "../../mesh/mesh.hpp"
#include "../../scalars/scalars.hpp"
#include "../../radiation/radiation.hpp"
#include "../../radiation/absorber/absorber.hpp"
#include "../reaction_network.hpp"
#include "verner_recombination.hpp"

VernerRecombination::VernerRecombination(std::string name,
  int neu_num, int ion_num, int elec_num, Real a_in, Real b_in, Real T0_in, Real T1_in, Real E_ion_in, Real E_neutral_in):
    Reaction(name)
{
  scalar_num = neu_num;
  ion_scalar_num = ion_num;
  electron_scalar_num = elec_num;

  a = a_in;
  b = b_in;
  T0 = T0_in;
  T1 = T1_in;

  E_ion = E_ion_in;
  E_neutral = E_neutral_in;
}

void VernerRecombination::react(int k, int j, int i) {
  PassiveScalars *ps = pmy_network->pscalars;

  Real T = pmy_network->temperature_(k,j,i);

  Real n_ion = ps->s(ion_scalar_num, k, j, i) / ps->mass(ion_scalar_num);
  Real n_elec = ps->s(electron_scalar_num, k, j, i)/ ps->mass(electron_scalar_num);

  Real alpha_1 = sqrt(T/T0);
  Real alpha_2 = pow(1+sqrt(T/T0), 1-b);
  Real alpha_3 = pow(1+sqrt(T/T1), 1+b);
  Real alpha = a * pow(alpha_1 * alpha_2 * alpha_3, -1.0); // cm3 s-1

  Real n_recomb = alpha * n_ion * n_elec; // s-1 cm-3

  Real e_recomb = (E_ion - E_neutral) * n_recomb; // erg s-1 cm-3

  int species[3] = {scalar_num, ion_scalar_num, electron_scalar_num};
  int sign[3]    = {+1, -1, -1};
  for (int l = 0; l < 3; ++l)
  {
    pmy_network->dn_rate(my_rxn_num, species[l], k, j, i) += sign[l] * n_recomb;
    pmy_network->jacobian(species[l], ion_scalar_num, k, j, i)      += sign[l] * alpha * n_elec;
    pmy_network->jacobian(species[l], electron_scalar_num, k, j, i) += sign[l] * alpha * n_ion;
  }

  pmy_network->de_rate(my_rxn_num, k, j, i) += e_recomb;
}