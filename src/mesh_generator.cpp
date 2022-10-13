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
#include "mesh_generator.hpp"

MeshGenerator::MeshGenerator(Real input_r_min, Real input_r_max, int input_n, ParameterInput *pin):
  r_min(input_r_min), r_max(input_r_max), n(input_n)
{
  // Values from parameter file
  outer_cell_ratio = pin->GetOrAddReal("problem","outer_cell_ratio",1.0);
  int turnover_maxit = pin->GetOrAddReal("problem","turnover_maxit",100);
  Real turnover_convergence = pin->GetOrAddReal("problem","turnover_convergence",1.e-6);
  Real turnover_initial_guess = pin->GetOrAddReal("problem","turnover_initial_guess",0.5);

  Real r_turnover_Rp = pin->GetReal("problem","r_turnover_Rp");
  Real Rp = pin->GetReal("problem","Rp");
  r_turnover = r_turnover_Rp * Rp;

  // calculate mesh spacing parameters
  // with edge cases
  if (r_turnover == r_min) {
    x_turnover = 0.0;
  }
  else if (r_turnover == r_max) {
    x_turnover = 1.0;
  }
  else {
    x_turnover = CalculateTurnoverPoint(turnover_initial_guess, turnover_maxit, turnover_convergence);
  }
  A2 = CalculateA2();
  B2 = CalculateB2();
}

// Functions for calculating parameters
Real MeshGenerator::ImplicitTurnoverPoint(Real xt) {
  Real c = (r_turnover - r_min)/(n * log(outer_cell_ratio));
  Real term1 = pow(outer_cell_ratio, ((1.0-xt)*n) );
  Real term2 = (r_max-r_turnover)/c * xt;
  return term1 - 1.0 - term2;
}

Real MeshGenerator::ImplicitTurnoverPointDerivative(Real xt) {
  Real c = (r_turnover - r_min)/(n * log(outer_cell_ratio));
  Real term1_a = -n * log(outer_cell_ratio);
  Real term1_b = pow(outer_cell_ratio, ((1.0-xt)*n) );
  Real term2 = (r_max-r_turnover)/c;
  return term1_a * term1_b - term2;
}

Real MeshGenerator::CalculateTurnoverPoint(Real x0, int maxit, Real convergence) {
  // update iterations
  maxit-=1;

  Real x1 = x0 - ImplicitTurnoverPoint(x0)/ImplicitTurnoverPointDerivative(x0);

  if (std::abs(x1 - x0) < convergence) {
    return x1;
  }

  // reached maxit and not converged
  if (maxit == 0) {
    std::stringstream msg;
    msg << "##### FATAL ERROR in MeshGenerator::CalculateTurnoverPoint" << std::endl
        << "      Reached max iterations without convergence." << std::endl
        << "      maxit       -- " << maxit << std::endl
        << "      abs(x1-x0)  -- " << abs(x1 - x0) << std::endl
        << "      convergence -- " << convergence << std::endl;
    ATHENA_ERROR(msg);
  }

  // more iterations to go, have not converged
  return CalculateTurnoverPoint(x1, maxit, convergence);
}

Real MeshGenerator::CalculateA2() {
  if (x_turnover == 0.0) {
    return (r_max-r_min)/(pow(outer_cell_ratio,n) -1.0);
  }
  Real c = (r_turnover - r_min)/(n * log(outer_cell_ratio));
  return c/x_turnover;
}

Real MeshGenerator::CalculateB2() {
  return r_turnover - A2;
}

// Functions for calculating mesh spacing
Real MeshGenerator::LowXFunction(Real x) {
  Real a = (r_turnover - r_min)/x_turnover;
  return a * x + r_min;
}

Real MeshGenerator::HighXFunction(Real x) {
  Real exp = n*(x-x_turnover);
  return A2 * pow(outer_cell_ratio, exp) + B2;
}

Real MeshGenerator::MeshSpacing(Real x) {
  if (x < x_turnover) {
    return LowXFunction(x);
  }
  else if (x == x_turnover) {
    return r_turnover;
  }
  else {
    return HighXFunction(x);
  }
}