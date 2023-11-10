// C/C++ headers
#include <sstream>
#include <string>

#include "../../math/core.h"
#include "../../math/interpolation.h"

// Athena++ header
#include "../../utils/utils.hpp"
#include "../../scalars/scalars.hpp"
#include "../../mesh/mesh.hpp"
#include "../../radiation/radiation.hpp"
#include "../reaction_network.hpp"
#include "../reaction.hpp"

#include "helium_recombination.hpp"

HeliumRecombination::HeliumRecombination(std::string name, std::vector<int> species, std::vector<Real> stoichiometry, std::string alpha_file_, std::string beta_file_, int file_column_index_):
  ReactionTemplate(name, species, stoichiometry)
{
  ReadAlphaBetaFiles(alpha_file_, beta_file_, file_column_index_);
}

HeliumRecombination::~HeliumRecombination() {

}

void HeliumRecombination::ReadAlphaBetaFiles(std::string alpha_file, std::string beta_file, int column_index)
{
  // Read alpha/beta files
  ReadDataTableForInterp(alpha_file, alpha_T, alpha_data, n_alpha_file, 0, column_index, true);
  ReadDataTableForInterp(beta_file, beta_T, beta_data, n_beta_file, 0, column_index, true);

  // Adjust values for units
  // file written as log10(T), (T^(1/2) * [alpha,beta])
  for (int i = 0; i < n_alpha_file; ++i)
  {
    alpha_T[i]    = pow(10, alpha_T[i]);
    alpha_data[i] = alpha_data[i]/pow(alpha_T[i],0.5);
  }

  for (int i = 0; i < n_beta_file; ++i)
  {
    beta_T[i]    = pow(10, beta_T[i]);
    beta_data[i] = alpha_data[i]/pow(beta_T[i],0.5);
  }
}

Real HeliumRecombination::alpha(Real T, int k, int j, int i) {
  if (T < alpha_T[0]) {
    ;
  }
  else if (T > alpha_T[n_alpha_file -1])
  {
    ;
  }

  Real alpha_val = interp1(T, alpha_data.data(), alpha_T.data(), n_alpha_file);
  return alpha_val;
}

Real HeliumRecombination::beta(Real T, int k, int j, int i) {
  if (T < beta_T[0]) {
    ;
  }
  else if (T > beta_T[n_beta_file -1])
  {
    ;
  }
  Real beta_val = interp1(T, beta_data.data(), beta_T.data(), n_beta_file);
  return (-pmy_network->boltzmann * T) * beta_val;
}