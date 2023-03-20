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
  void react(int k, int j, int i);
protected:
  int scalar_num, ion_scalar_num, electron_scalar_num;
  Absorber *pmy_abs;
};

class He_23S_recombination: public Reaction {
public:
  He_23S_recombination(ReactionNetwork *pnetwork, std::string name, int neu_num, int ion_num, int elec_num);
  void react(int k, int j, int i);
protected:
  int scalar_num, ion_scalar_num, electron_scalar_num;
  Absorber *pmy_abs;
};

class He_e_collisions: public Reaction {
public:
  He_e_collisions(ReactionNetwork *pnetwork, std::string data_file, std::string name, int singlet_num, int triplet_num, int elec_num);
  ~He_e_collisions();
  void react(int k, int j, int i);
protected:
  int n_col_type = 4;

  int singlet_scalar_num, triplet_scalar_num, electron_scalar_num;
  AthenaArray<Real> file_data;

  int n_file;

  Real *file_T;
  Real **file_q;
  Real **file_L;
  enum {x11 = 0, x13 = 1, x31 = 2, x33 = 3};

  Real *temp_q;
  Real *temp_L;
};