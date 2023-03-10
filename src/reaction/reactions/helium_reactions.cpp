// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "../../math/core.h"
#include "../../math/interpolation.h"

#include "../../mesh/mesh.hpp"
#include "../../utils/utils.hpp"
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
    pmy_network->dn_rate(my_rxn_num, species[l], k, j, i) += sign[l] * n_recomb;
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
    pmy_network->dn_rate(my_rxn_num, species[l], k, j, i) += sign[l] * n_recomb;
    pmy_network->jacobian(species[l], ion_scalar_num, k, j, i)  += sign[l] * alpha_3 * n_elec;
    pmy_network->jacobian(species[l], electron_scalar_num, k, j, i) += sign[l] * alpha_3 * n_ion;
  }

  pmy_network->de_rate(my_rxn_num, k, j, i) -= e_recomb;
}

He_e_collisions::He_e_collisions(ReactionNetwork *pnetwork, std::string data_file, std::string name, int singlet_num, int triplet_num, int elec_num):
    Reaction(pnetwork, name)
{
  ReadDataTable(file_data, data_file);

  singlet_scalar_num = singlet_num;
  triplet_scalar_num = triplet_num;
  electron_scalar_num = elec_num;

  n_file = file_data.GetDim1();
  file_T = new Real[n_file];
  temp_q = new Real[n_col_type];
  temp_L = new Real[n_col_type];

  file_q = new Real*[n_col_type];
  file_L = new Real*[n_col_type];

  for (int i = 0; i < n_col_type; ++i)
  {
    file_q[i] = new Real[n_file];
    file_L[i] = new Real[n_file];
  }

  // File Columns
  // 0 -- T
  // 1-4  -- q11, q13, q31, q33
  // 5-8  -- L11, L13, L31, L33
  for (int i = 0; i < n_file; ++i)
  {
    file_T[i] = file_data(0,i);
    for (int j = 0; j < n_col_type; ++j)
    {
      file_q[j][i] = file_data(1+j, i);
      file_L[j][i] = file_data(5+j, i);
    }
  }
}

He_e_collisions::~He_e_collisions() {
  delete[] file_T;

  for (int i = 0; i < n_col_type; ++i)
  {
    delete[] file_q[i];
    delete[] file_L[i];
  }

  delete[] file_q;
  delete[] file_L;

  delete[] temp_q;
  delete[] temp_L; 
}

void He_e_collisions::react(int k, int j, int i) {
  PassiveScalars *ps = pmy_network->pscalars;
  Real T = pmy_network->temperature_(k,j,i);

  Real n_singlet = ps->s(singlet_scalar_num, k, j, i) / ps->mass(singlet_scalar_num);
  Real n_triplet = ps->s(triplet_scalar_num, k, j, i) / ps->mass(triplet_scalar_num);
  Real n_elec    = ps->s(electron_scalar_num, k, j, i)/ ps->mass(electron_scalar_num);

  Real k_T = (pmy_network->boltzmann * T);
  Real T_dependance = std::sqrt((pmy_network->ry_conversion * pmy_network->eV_conversion)/k_T);
  T_dependance = T_dependance * exp(-1.0/k_T);

  for (int ii = 0; ii < n_col_type; ++ii)
  {
    temp_q[ii] = interp1(T, file_q[ii], file_T, n_file) * T_dependance; // cm3/s
    temp_L[ii] = interp1(T, file_L[ii], file_T, n_file) * T_dependance; // erg cm3 /s
  }
  
 Real singlet_reaction_factor = (n_singlet * n_elec);
 Real triplet_reaction_factor = (n_triplet * n_elec);
 
 // ----- singlet -> singlet collisions
 // pmy_network->dn_rate(my_rxn_num, singlet_scalar_num, k, j, i) -= temp_q[x11] * singlet_reaction_factor;
 // pmy_network->dn_rate(my_rxn_num, singlet_scalar_num, k, j, i) += temp_q[x11] * singlet_reaction_factor;
 pmy_network->de_rate(my_rxn_num, k, j, i) -= temp_L[x11] * singlet_reaction_factor;

 // ----- singlet -> triplet collisions
 pmy_network->dn_rate(my_rxn_num, singlet_scalar_num, k, j, i) -= temp_q[x13] * singlet_reaction_factor;
 pmy_network->dn_rate(my_rxn_num, triplet_scalar_num, k, j, i) += temp_q[x13] * singlet_reaction_factor;
 pmy_network->de_rate(my_rxn_num, k, j, i) -= temp_L[x13] * singlet_reaction_factor;

 // ----- triplet -> singlet collisions
 pmy_network->dn_rate(my_rxn_num, triplet_scalar_num, k, j, i) -= temp_q[x31] * triplet_reaction_factor;
 pmy_network->dn_rate(my_rxn_num, singlet_scalar_num, k, j, i) += temp_q[x31] * triplet_reaction_factor;
 pmy_network->de_rate(my_rxn_num, k, j, i) -= temp_L[x31] * triplet_reaction_factor;

 // ----- triplet -> triplet collisions
 // pmy_network->dn_rate(my_rxn_num, triplet_scalar_num, k, j, i) -= temp_q[x33] * triplet_reaction_factor;
 // pmy_network->dn_rate(my_rxn_num, triplet_scalar_num, k, j, i) += temp_q[x33] * triplet_reaction_factor;
 pmy_network->de_rate(my_rxn_num, k, j, i) -= temp_L[x33] * triplet_reaction_factor;
}