// C/C++ headers
#include <string>

// Athena++ headers
#include "../../athena.hpp"

class MeshBlock;
class ReactionNetwork;
class Reaction;
class Absorber;
class PassiveScalars;

class Photoionization: public Reaction {
public:
  Photoionization(ReactionNetwork *pnetwork, std::string name, Absorber *pabs,
    int num, int ion_num);
  void react(AthenaArray<Real> &dn_rate, AthenaArray<Real> &de_rate, int k, int j, int i);
protected:
  int scalar_num, ion_scalar_num;
  Real ionization_energy;
  Absorber *pmy_abs;
};