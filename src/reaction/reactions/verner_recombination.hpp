// C/C++ headers
#include <string>

// Athena++ headers
#include "../../athena.hpp"
#include "../reaction_network.hpp"

class MeshBlock;
class ReactionNetwork;
class Reaction;
class Absorber;
class PassiveScalars;

class VernerRecombination: public Reaction {
public:
  VernerRecombination(std::string name, int neu_num, int ion_num, int elec_num, Real a_in, Real b_in, Real T0_in, Real T1_in, Real E_ion_in, Real E_neutral_in);
  void react(int k, int j, int i);

protected:
  int scalar_num, ion_scalar_num, electron_scalar_num;

  Real a,b,T1,T0;
  Real E_neutral, E_ion;

  Absorber *pmy_abs;
};