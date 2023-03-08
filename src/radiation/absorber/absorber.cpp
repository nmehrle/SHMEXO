// C/C++ headers
#include <sstream>

// Athena++ headers
#include "../../athena.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../../scalars/scalars.hpp"
#include "../../reaction/reaction_network.hpp"
#include "../radiation.hpp"
#include "absorber.hpp"

Absorber::Absorber(RadiationBand *pband, int my_scalar_number, ParameterInput *pin) {
  pmy_band = pband;
  MeshBlock *pmb = pmy_band->pmy_rad->pmy_block;
  ps = pmb->pscalars;

  // associate a scalar to get number densities
  SetScalar(my_scalar_number);

  // initialize absorber data
  crossSection.NewAthenaArray(pband->nspec);

  // assign default h, q values
  h.NewAthenaArray(pband->nspec);
  q.NewAthenaArray(pband->nspec);

  for (int n = 0; n < pband->nspec; ++n)
  {
    crossSection(n) = 0.0;
    h(n) = 0.0;
    q(n) = 1.0;
  }

  speed_of_light = pin->GetReal("problem", "speed_of_light");
  planck_constant = pin->GetReal("problem", "planck_constant");

  absorptionCoefficient.NewAthenaArray(pband->nspec, pmb->ncells3, pmb->ncells2, pmb->ncells1);
  energyAbsorbed.NewAthenaArray(pband->nspec, pmb->ncells3, pmb->ncells2, pmb->ncells1);
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
    mass = ps->mass(n);
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

void Absorber::CalculateAbsorptionCoefficient(AthenaArray<Real> const& prim, AthenaArray<Real> const& cons_scalar, int n, int k, int j, int i) {
  Spectrum *spec = pmy_band->spec;
  Real number_density;

  if (scalar_num == -1) {
    if (mass <= 0.0) {
      std::stringstream msg;
      msg << "##### Fatal Error in Absorber::CalculateAbsorptionCoefficient." << std::endl
          << "Default absorption coefficient calculation requires number density, which requires a particle mass." << std::endl
          << " No valid mass has been assigned. Consider using scalars to assign a mass or using a non-default CalculateAbsorptionCoefficient function.";
      ATHENA_ERROR(msg);
    }
    else {
      number_density = prim(IDN,k,j,i)/mass;
    }
  }
  else {
    number_density = cons_scalar(scalar_num,k,j,i)/mass;
  }

  for (int n = 0; n < pmy_band->nspec; ++n)
  {
    absorptionCoefficient(n,k,j,i) = number_density * crossSection(n);
  }
}