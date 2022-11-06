// C/C++ headers
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <cassert>  // assert

// Athena++ headers
#include "../athena_arrays.hpp"
#include "../globals.hpp"
#include "absorber.hpp"
#include "radiation.hpp"
#include "../mesh/mesh.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "null_absorber.hpp"

Real NullAbsorber::AbsorptionCoefficient(Real wave, Real const prim[]) const
{
  return 0;
}
