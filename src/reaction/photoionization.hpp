// C/C++ headers
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "reaction.hpp"

class Reaction;
class Absorber;
class IonizingAbsorber;

class Photoionization: public Reaction {
public:
  Photoionization(std::string name, IonizingAbsorber *pabs, int elec_num);
  void react(int k, int j, int i);

protected:
  int scalar_num_;
  int ion_num_;
  int electron_num_;

  Real ionization_energy;
  IonizingAbsorber *pmy_abs;
};