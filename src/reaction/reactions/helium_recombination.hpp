// C/C++ headers
#include <string>

// Athena++ headers
#include "../../athena.hpp"
#include "../reaction_network.hpp"
#include "../reaction.hpp"

class MeshBlock;
struct float_triplet;

class HeliumRecombination: public ReactionTemplate {
public:
  HeliumRecombination(std::string name, std::vector<int> species, std::vector<Real> stoichiometry, std::string alpha_file_, std::string beta_file_, int file_column_index_);
  ~HeliumRecombination();

  Real alpha(Real T, int k, int j, int i);
  Real beta(Real T, int k, int j, int i);

protected:
  int n_alpha_file, n_beta_file;
  std::vector<Real> alpha_T, alpha_data;
  std::vector<Real> beta_T, beta_data;

  void ReadAlphaBetaFiles(std::string alpha_file, std::string beta_file, int column_index);
};
