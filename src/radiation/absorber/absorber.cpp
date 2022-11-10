// C/C++ headers
#include <sstream>

// Athena++ headers
#include "../../athena.hpp"
#include "absorber.hpp"
#include "../radiation.hpp"
#include "../../mesh/mesh.hpp"
#include "../../scalars/scalars.hpp"

Absorber::Absorber(RadiationBand *pband) {
  pmy_band = pband;
  ps = pmy_band->pmy_rad->pmy_block->pscalars;
}

Absorber::~Absorber() {
}

void Absorber::SetScalar(int n) {
  scalar_num = n;
  // check n is valid value
  if (n == -1) {
    // absorber configured without mass value
    // set negative mass to cause error if mass is requested
    mass = -1.;
  }
  else if (n < NSCALARS) {
    // Check scalar has valid mass
    mass = ps->m(n);
    if (mass <= 0) {
      std::stringstream msg;
      msg << "##### Fatal Error in Absorber Constructor." << std::endl
          << "Scalar mass for scalar ("<<n<<") must be a positive value." << std::endl
          << "Currently ("<<mass<<")";
      ATHENA_ERROR(msg);  
    }
  }
  else {
    std::stringstream msg;
    msg << "##### Fatal Error in Absorber Constructor." << std::endl
        << "Scalar number value ("<<n<<") must be either -1 (to use complete density)"
        << std::endl
        << " or between 0 and NSCALARS ("<<NSCALARS<<").";
    ATHENA_ERROR(msg);
  }
}