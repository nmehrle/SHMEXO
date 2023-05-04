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
  He_recombination(std::string name, int neu_num, int ion_num, int elec_num);
  void react(int k, int j, int i);
protected:
  int scalar_num, ion_scalar_num, electron_scalar_num;
  Absorber *pmy_abs;
};

class He_23S_recombination: public Reaction {
public:
  He_23S_recombination(std::string name, int neu_num, int ion_num, int elec_num);
  void react(int k, int j, int i);
protected:
  int scalar_num, ion_scalar_num, electron_scalar_num;
  Absorber *pmy_abs;
};

class He_e_collisions: public Reaction {
public:
  He_e_collisions(std::string data_file, std::string name, int singlet_num, int triplet_num, int elec_num);
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

class He_triplet_decay: public Reaction {
public:
  He_triplet_decay(std::string name, int He_trip_num, int He_singlet_num, Real E_triplet);

  void react(int k, int j, int i);
protected:
  int He_trip_num_, He_singlet_num_;
  Real E_triplet_;
};


class CollisionalRelaxation_HeH: public Reaction {
public:
  // CollisionalRelaxation_HeH(ReactionNetwork *pnetwork, std::string name, int HEtrip_num, int H_num, int HE_num, Real E_reactant, Real E_product);

  CollisionalRelaxation_HeH(std::string name, int r1_num, int r2_num, int p1_num, int p2_num, int p3_num, Real E_reactant, Real E_product);
  void react(int k, int j, int i);
protected:
  int r1_num_;
  int r2_num_;
  int p1_num_;
  int p2_num_;
  int p3_num_;

  Real E_reactant_;
  Real E_product_;
};