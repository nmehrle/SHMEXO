// C/C++ headers
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <cassert>  // assert

// Athena++ headers
#include "absorber.hpp"
#include "null_absorber.hpp"

Real NullAbsorber::AbsorptionCoefficient(Real wave, Real const prim[]) const
{
  return 0;
}
