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

Real TripletHydrogenCollision::alpha(Real T, int k, int j, int i) {
  return 5.0E-10;
}

HeElectronCollisions::HeElectronCollisions(std::string name, std::vector<int> species, std::vector<Real> stoichiometry, std::string data_file_, int alpha_column_, int beta_column_):
  ReactionTemplate(name, species, stoichiometry)
{
  // Read alpha/beta files
  ReadDataTableForInterp(data_file_, file_T, alpha_data, n_file, 0, alpha_column_, true);
  ReadDataTableForInterp(data_file_, file_T, beta_data, n_file, 0, beta_column_, true);
}

Real HeElectronCollisions::TempScalingFactor(Real T) {
  Real k_T_ev = (pmy_network->boltzmann * T) / pmy_network->eV_conversion;
  Real T_dependance = std::sqrt((pmy_network->ry_conversion / k_T_ev));
  return T_dependance;
}

Real HeElectronCollisions::alpha(Real T, int k, int j, int i) {
  Real alpha_val;
  if (T < file_T[0]) {
    return 0;
  }
  else if (T > file_T[n_file -1])
  {
    alpha_val = alpha_data[n_file-1];
  }
  else {
    alpha_val = interp1(T, alpha_data.data(), file_T.data(), n_file);
  }

  return alpha_val * TempScalingFactor(T);
}

Real HeElectronCollisions::beta(Real T, int k, int j, int i) {
  Real beta_val;
  if (T < file_T[0]) {
    return 0;
  }
  else if (T > file_T[n_file -1])
  {
    beta_val = beta_data[n_file-1];
  }
  else {
    beta_val = interp1(T, beta_data.data(), file_T.data(), n_file);
  }

  return -beta_val * TempScalingFactor(T);
}
