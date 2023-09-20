// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "../../math/core.h"
#include "../../math/interpolation.h"

#include "../../mesh/mesh.hpp"
#include "../../utils/utils.hpp"
#include "../../scalars/scalars.hpp"
#include "../reaction_network.hpp"
#include "../reaction.hpp"
#include "collisions.hpp"

Real ChargeExchangeHeliumToHyd::alpha(Real T, int k, int j, int i) {
  return 1.25E-15 * pow(300.0/T, -0.25);
}

Real ChargeExchangeHydToHelium::alpha(Real T, int k, int j, int i) {
  Real t1 = 1.75E-11 * pow(300.0/T, 0.75);
  Real t2 = exp(-128000.0/T);
  return t1 * t2;
}

// Real ChargeExchangeHeliumToHyd::beta(Real T, int k, int j, int i) {
//   PassiveScalars *ps = pmy_network->pscalars;

//   Real E_reactant = ps->energy(r1_num_) + ps->energy(r2_num_);
//   Real E_product = ps->energy(p1_num_) + ps->energy(p2_num_);

//   return this->alpha(T, k, j, i) * (E_reactant - E_product);
// }

// Real ChargeExchangeHydToHelium::beta(Real T, int k, int j, int i) {
//   PassiveScalars *ps = pmy_network->pscalars;

//   Real E_reactant = ps->energy(r1_num_) + ps->energy(r2_num_);
//   Real E_product = ps->energy(p1_num_) + ps->energy(p2_num_);

//   return this->alpha(T, k, j, i) * (E_reactant - E_product);
// }

Real TripletHydrogenCollision::alpha(Real T, int k, int j, int i) {
  return 5.0E-10;
}

// todo
// Real TripletHydrogenCollision::beta(Real T, int k, int j, int i) {
//   PassiveScalars *ps = pmy_network->pscalars;

//   Real E_reactant = ps->energy(r1_num_) + ps->energy(r2_num_);
//   Real E_product = ps->energy(p1_num_) + ps->energy(p2_num_) + ps->energy(p3_num_);

//   return this->alpha(T, k, j, i) * (E_reactant - E_product);
// }

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

  Real k_T_ev = (pmy_network->boltzmann * T) / pmy_network->eV_conversion;
  Real T_dependance = std::sqrt((pmy_network->ry_conversion / k_T_ev));
  T_dependance = T_dependance * exp(-1.0/k_T_ev);

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