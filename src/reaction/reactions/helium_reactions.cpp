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

He_recombination::He_recombination(std::string name,
  int neu_num, int ion_num, int elec_num):
    Reaction(name)
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

He_23S_recombination::He_23S_recombination(std::string name,
  int neu_num, int ion_num, int elec_num):
    Reaction(name)
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

He_e_collisions::He_e_collisions(std::string data_file, std::string name, int singlet_num, int triplet_num, int elec_num):
    Reaction(name)
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
 pmy_network->de_rate(my_rxn_num, k, j, i) -= temp_L[x11] * singlet_reaction_factor;

 // ----- singlet -> triplet collisions
 int s2t_species[2] = {singlet_scalar_num, triplet_scalar_num};
 int s2t_sign[2]   = {-1, +1};

 for (int l = 0; l < 2; ++l) {
   pmy_network->dn_rate(my_rxn_num, s2t_species[l], k, j, i) += s2t_sign[l] * temp_q[x13] * singlet_reaction_factor;
   pmy_network->jacobian(s2t_species[l], singlet_scalar_num, k, j ,i) += s2t_sign[l] * temp_q[x13] * n_elec;
   pmy_network->jacobian(s2t_species[l], electron_scalar_num, k, j ,i) += s2t_sign[l] * temp_q[x13] * n_singlet;
 }

 // pmy_network->dn_rate(my_rxn_num, singlet_scalar_num, k, j, i) -= temp_q[x13] * singlet_reaction_factor;
 // pmy_network->dn_rate(my_rxn_num, triplet_scalar_num, k, j, i) += temp_q[x13] * singlet_reaction_factor;

 pmy_network->de_rate(my_rxn_num, k, j, i) -= temp_L[x13] * singlet_reaction_factor;

 // ----- triplet -> singlet collisions
 int t2s_species[2] = {triplet_scalar_num, triplet_scalar_num};
 int t2s_sign[2]   = {-1, +1};

 for (int l = 0; l < 2; ++l) {
   pmy_network->dn_rate(my_rxn_num, t2s_species[l], k, j, i) += t2s_sign[l] * temp_q[x31] * triplet_reaction_factor;
   pmy_network->jacobian(t2s_species[l], triplet_scalar_num, k, j ,i) += t2s_sign[l] * temp_q[x31] * n_elec;
   pmy_network->jacobian(t2s_species[l], electron_scalar_num, k, j ,i) += t2s_sign[l] * temp_q[x31] * n_triplet;
 }

 // pmy_network->dn_rate(my_rxn_num, triplet_scalar_num, k, j, i) -= temp_q[x31] * triplet_reaction_factor;
 // pmy_network->dn_rate(my_rxn_num, singlet_scalar_num, k, j, i) += temp_q[x31] * triplet_reaction_factor;

 pmy_network->de_rate(my_rxn_num, k, j, i) -= temp_L[x31] * triplet_reaction_factor;

 // ----- triplet -> triplet collisions
 pmy_network->de_rate(my_rxn_num, k, j, i) -= temp_L[x33] * triplet_reaction_factor;
}


He_triplet_decay::He_triplet_decay(std::string name, int He_trip_num, int He_singlet_num, Real E_triplet):
    Reaction(name)
{
  He_trip_num_    = He_trip_num;
  He_singlet_num_ = He_singlet_num;
  E_triplet_ = E_triplet;
}

void He_triplet_decay::react(int k, int j, int i) {
  PassiveScalars *ps = pmy_network->pscalars;

  Real n_triplet = ps->s(He_trip_num_, k, j, i) / ps->mass(He_trip_num_);
  Real alpha = 1.272E-4; // s^-1 Lampon+ 2020, Table 2
  Real n_rxn = alpha * n_triplet;
  Real e_rxn = n_rxn * E_triplet_;

  pmy_network->dn_rate(my_rxn_num, He_trip_num_, k, j, i) -= n_rxn;
  pmy_network->dn_rate(my_rxn_num, He_singlet_num_, k, j, i) += n_rxn;

  pmy_network->jacobian(He_trip_num_, He_trip_num_, k, j, i) -= alpha;
  pmy_network->jacobian(He_singlet_num_, He_trip_num_, k, j, i) += alpha;

  pmy_network->de_rate(my_rxn_num, k, j, i) += e_rxn;
}


//HETRIP, HYD, HE, HPLUS, ELEC
CollisionalRelaxation_HeH::CollisionalRelaxation_HeH(std::string name, int r1_num, int r2_num, int p1_num, int p2_num, int p3_num, Real E_reactant, Real E_product):
  Reaction(name)
{
  r1_num_ = r1_num;
  r2_num_ = r2_num;
  p1_num_ = p1_num;
  p2_num_ = p2_num;
  p3_num_ = p3_num;

  E_reactant_ = E_reactant;
  E_product_  = E_product;
}

void CollisionalRelaxation_HeH::react(int k, int j, int i) {
  PassiveScalars *ps = pmy_network->pscalars;

  Real n_r1 = ps->s(r1_num_, k, j, i)/ps->mass(r1_num_);
  Real n_r2 = ps->s(r2_num_, k, j, i)/ps->mass(r2_num_);

  Real alpha = 5.0E-10;

  Real n_rxn = alpha * n_r1 * n_r2;
  Real e_rxn = (E_reactant_ - E_product_) * n_rxn;

  int species[5] = {r1_num_, r2_num_, p1_num_, p2_num_, p3_num_};
  int sign[5] = {-1, -1, +1, +1, +1};

  for (int l = 0; l < 5; ++l) {
    pmy_network->dn_rate(my_rxn_num, species[l], k, j, i) += sign[l] * n_rxn;
    pmy_network->jacobian(species[l], r1_num_, k, j, i) += sign[l] * alpha * n_r2;
    pmy_network->jacobian(species[l], r2_num_, k, j, i) += sign[l] * alpha * n_r1;
  }

  pmy_network->de_rate(my_rxn_num, k, j, i) += e_rxn;
}