// C/C++ headers
#include <string>

// Athena++ headers
#include "../../athena.hpp"
#include "../../scalars/scalars.hpp"
#include "../reaction_network.hpp"
#include "../reaction.hpp"


class CollisionalIonization: public ReactionTemplate {
public:
  CollisionalIonization(std::string name, std::vector<int> species, std::vector<Real> stoichiometry, int mode);
  Real alpha(Real T, int k, int j, int i);
  Real beta(Real T, int k, int j, int i);

protected:
  Real alpha_coeff, beta_coeff, exponent;
  Real Tmin;
};