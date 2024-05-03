// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "../../scalars/scalars.hpp"
#include "../reaction_network.hpp"
#include "../reaction.hpp"
#include "additional_reactions.hpp"

LyaCooling::LyaCooling(std::string name, int scalar_num, int elec_num):
    ReactionTemplate(name, {scalar_num, elec_num}, {-1, -1})
{
  ;
}

Real LyaCooling::beta(Real T, int k, int j, int i) {
  Real cen_factor = std::pow(1+std::sqrt(T/1E5),-1);
  return -7.5E-19 * exp(-118348./T) * cen_factor;
}

Real HeliumTripletDecay::alpha(Real T, int k, int j, int i) {
  return 1.272E-4; // s^-1 Lampon+ 2020, Table 2
}

// is it zero?
// Energy released as 19eV photon
// What does it do lol
Real HeliumTripletDecay::beta(Real T, int k, int j, int i) {
  // PassiveScalars *ps = pmy_network->pscalars;

  // Real E_reactant = ps->energy(r1_num_);
  // Real E_product = ps->energy(p1_num_);

  // return this->alpha(T, k, j, i) * (E_reactant - E_product);
  return 0;
}

Bremsstrahlung::Bremsstrahlung(std::string name, int elec_num,
  int HII_num, int HeII_num, int HeIII_num):
    Reaction(name)
{
  elec_num_ = elec_num;
  HII_num_ = HII_num;
  HeII_num_ = HeII_num;
  HeIII_num_ = HeIII_num;
}

void Bremsstrahlung::react(int k, int j, int i) {
  PassiveScalars *ps = pmy_network->pscalars;
  Real T = pmy_network->temperature_(k,j,i);

  Real n_e = ps->s(elec_num_, k, j, i)/ps->mass(elec_num_);
  Real n_HII = ps->s(HII_num_, k, j, i)/ps->mass(HII_num_);
  Real n_HeII = ps->s(HeII_num_, k, j, i)/ps->mass(HeII_num_);
  Real n_HeIII = ps->s(HeIII_num_, k, j, i)/ps->mass(HeIII_num_);

  // Black 1981
  // https://articles.adsabs.harvard.edu/pdf/1981MNRAS.197..553B
  Real e_rxn = -1.42E-27 * pow(T,1/2.) * (n_HII + n_HeII + 4*n_HeIII)*n_e;

  pmy_network->de_rate(my_rxn_num, k, j, i) += e_rxn;

  // if (this->my_reemission != nullptr)
  //   this->ProducePhotons(n_rxn, T, k, j, i);
}

HeIIElecCollisions::HeIIElecCollisions(std::string name, int scalar_num, int elec_num):
    ReactionTemplate(name, {scalar_num, elec_num}, {-1, -1})
{
  ;
}

Real HeIIElecCollisions::beta(Real T, int k, int j, int i) {
  Real cen_factor = std::pow(1+std::sqrt(T/1E5),-1);
  return -5.54e-17 * pow(T,-0.397) * exp(-473638./T) * cen_factor;
}

HeIElecCollisions::HeIElecCollisions(std::string name, int scalar_num, int elec_num):
    ReactionTemplate(name, {scalar_num, elec_num}, {-1, -2})
{
  ;
}

Real HeIElecCollisions::alpha(Real T, int k, int j, int i) {
  return 0;
}

Real HeIElecCollisions::beta(Real T, int k, int j, int i) {
  Real cen_factor = std::pow(1+std::sqrt(T/1E5),-1);
  return -9.1e-27 * pow(T,-0.1687) * exp(-13179./T) * cen_factor;
}