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

class He_recombination: public Reaction {
public:
  He_recombination(ReactionNetwork *pnetwork, std::string name, int neu_num, int ion_num, int elec_num);
  void react(AthenaArray<Real> &dn_rate, AthenaArray<Real> &de_rate, int k, int j, int i);

  std::string my_name;
protected:
  int scalar_num, ion_scalar_num, electron_scalar_num;
  Absorber *pmy_abs;
};

class He_23S_recombination: public Reaction {
public:
  He_23S_recombination(ReactionNetwork *pnetwork, std::string name, int neu_num, int ion_num, int elec_num);
  void react(AthenaArray<Real> &dn_rate, AthenaArray<Real> &de_rate, int k, int j, int i);

  std::string my_name;
protected:
  int scalar_num, ion_scalar_num, electron_scalar_num;
  Absorber *pmy_abs;
};

// class He_e_collisions: public Reaction {
// public:
//   He_e_collisions(ReactionNetwork *pnetwork, std::string name, int neu_num, int ion_num, int elec_num);
//   void react(AthenaArray<Real> &dn_rate, AthenaArray<Real> &de_rate, int k, int j, int i);

//   std::string my_name;
// protected:
//   int scalar_num, ion_scalar_num, electron_scalar_num;
//   Absorber *pmy_abs;
// }