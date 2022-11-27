// C/C++ headers
#include <sstream>

// Athena++ header
#include "rtsolver.hpp"
#include "simple_rtsolver.hpp"

#include "../../parameter_input.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../radiation.hpp"
#include "../absorber/absorber.hpp"


RTSolver::RTSolver(RadiationBand *pband, ParameterInput *pin)
{
  pmy_band = pband;
}