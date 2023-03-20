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

class H_recombination: public Reaction {
public:
  H_recombination(ReactionNetwork *pnetwork, std::string name, int neu_num, int ion_num, int elec_num);
  void react(int k, int j, int i);

protected:
  int scalar_num, ion_scalar_num, electron_scalar_num;
  Absorber *pmy_abs;
};

class Lya_cooling: public Reaction {
public:
  Lya_cooling(ReactionNetwork *pnetwork, std::string name, int neu_num, int ion_num, int elec_num);
  void react(int k, int j, int i);

protected:
  int scalar_num, ion_scalar_num, electron_scalar_num;
  Absorber *pmy_abs;
};