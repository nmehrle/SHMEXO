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

class Photoionization: public Reaction {
public:
  Photoionization(ReactionNetwork *pnetwork, std::string name, Absorber *pabs, 
    int ion_num, int elec_num, Real ion_eng);
  void react(AthenaArray<Real> &dn_rate, AthenaArray<Real> &de_rate, int k, int j, int i);

  std::string my_name;
protected:
  int scalar_num, ion_scalar_num, electron_scalar_num;
  Real ionization_energy;
  Absorber *pmy_abs;
};