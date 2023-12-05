// C/C++ headers
#include <sstream>
#include <string>

// Athena++ header
#include "reaction.hpp"
#include "reaction_network.hpp"
#include "../reemission/reemission.hpp"
#include "../scalars/scalars.hpp"
#include "../radiation/radiation.hpp"
#include "../mesh/mesh.hpp"

Reaction::Reaction(std::string name):
  my_name(name)
{
  pmy_network = nullptr;
  next = nullptr;
  prev = nullptr;
}

NullReaction::NullReaction():
  Reaction("Null Reaction")
{
}

ReactionTemplate::ReactionTemplate(std::string name, std::vector<int> species, std::vector<Real> stoichiometry):
  Reaction(name), num_species_(species.size()), number_density_(species.size(), 0.0),
  species_(species), stoichiometry_(stoichiometry)
{
  if (species.size() != stoichiometry.size()) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ReactionTemplate::Constructor" << std::endl
        << "Length of input species array (" << species.size()
        << ") does not match length of input stoichiometry array (" << stoichiometry.size()
        << ") " << std::endl << "For reaction " << name << "." << std::endl;
    ATHENA_ERROR(msg);
  }
}

ReactionTemplate::~ReactionTemplate() {
  if (my_reemission != nullptr)
    delete my_reemission;
}

void ReactionTemplate::AssignReemission(Reemission *reemission_) {
  my_reemission = reemission_;
}

void ReactionTemplate::react(int k, int j, int i) {
  PassiveScalars *ps = pmy_network->pscalars;
  Real T = pmy_network->temperature_(k,j,i);

  Real reactant_number_density = 1;
  // Compute total reactant density for rate calculation
  for (int r = 0; r < num_species_; ++r)
  {
    // check if this species is reactant (i.e. negative stoichiometry)
    if (stoichiometry_[r] >= 0)
      continue;
    number_density_[r] = ps->s(species_[r], k, j, i)/ps->mass(species_[r]);
    reactant_number_density *= pow(number_density_[r], abs(stoichiometry_[r]));
  }

  Real n_rxn = reactant_number_density * this->alpha(T, k, j, i);
  Real e_rxn = reactant_number_density * this->beta(T, k, j, i);

  // Record dn(s)/dt for each species
  for (int s = 0; s < num_species_; ++s)
  {
    pmy_network->dn_rate(my_rxn_num, species_[s], k, j, i) += stoichiometry_[s]*n_rxn;

    // Record jacobian elements d2 n(s) / (dt dn(r))
    for (int r = 0; r < num_species_; ++r)
    {
      if (stoichiometry_[r] >= 0)
        continue;
      pmy_network->jacobian(species_[s], species_[r], k, j, i) += stoichiometry_[s] * n_rxn * (abs(stoichiometry_[r])/number_density_[r]);
    }
  }

  pmy_network->de_rate(my_rxn_num, k, j, i) += e_rxn;

  if (this->my_reemission != nullptr)
    this->ProducePhotons(n_rxn, T, k, j, i);
}

void ReactionTemplate::ProducePhotons(Real n_rxn, Real T, int k, int j, int i) {
  Radiation *prad = pmy_network->prad;
  RadiationBand *pband;
  Real wave, energy, frac;

  for (int b = 0; b < prad->nbands; ++b)
  {
    pband = prad->my_bands(b);
    for (int n = 0; n < pband->nspec; ++n)
    {
      wave = pband->spec[n].wave;
      frac = my_reemission->ReemissionFunction(b, n, T, k, j, i);

      // photon energy
      energy = (pmy_network->planck_constant * pmy_network->speed_of_light)/wave;

      // energy converted to photons
      energy = energy * n_rxn * frac;

      pband->emission_coefficient(n, k, j, i) += energy;
    }
  }
}

Real ReactionTemplate::beta(Real T, int k, int j, int i) {
  PassiveScalars *ps = pmy_network->pscalars;

  Real E_reactant = 0;
  Real E_product  = 0;

  for (int s = 0; s < num_species_; ++s)
  {
    // Check if reactant (negative stoichiometry)
    if (stoichiometry_[s] < 0)
      E_reactant += ps->energy(species_[s]);
    // Check if product (positive stoichiometry)
    else if (stoichiometry_[s] > 0)
      E_product += ps->energy(species_[s]);
  }

  return this->alpha(T, k, j, i) * (E_reactant - E_product);
}