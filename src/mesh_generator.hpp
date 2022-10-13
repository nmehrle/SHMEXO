// C/C++ header
#include <ctime>
#include <cmath>

// Athena++ headers
#include "athena.hpp"
#include "athena_arrays.hpp"
#include "parameter_input.hpp"
#include "coordinates/coordinates.hpp"
#include "mesh/mesh.hpp"
#include "globals.hpp"
#include "utils/utils.hpp"


class MeshGenerator {
public:
  MeshGenerator(){};
  MeshGenerator(Real input_r_min, Real input_r_max, int input_n, ParameterInput *pin);

  Real MeshSpacing(Real x);

protected:
  // mesh_size params
  Real r_min, r_max;
  int n;

  // input params
  Real outer_cell_ratio, r_turnover;

  // calculated params
  Real A2, B2, x_turnover;

  // Functions for calculating parameters
  Real ImplicitTurnoverPoint(Real xt);
  Real ImplicitTurnoverPointDerivative(Real xt);
  Real CalculateTurnoverPoint(Real x0, int maxit, Real convergence);
  Real CalculateA2();
  Real CalculateB2();

  // Functions for calculating mesh spacing
  Real LowXFunction(Real x);
  Real HighXFunction(Real x);
};