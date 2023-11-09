// C/C++ headers
#include <sstream>

// Athena++ headers
#include "../athena.hpp"
#include "../parameter_input.hpp"
#include "../radiation/radiation.hpp"
#include "reemission.hpp"

Reemission::Reemission(ParameterInput *pin, Radiation *prad_)
{
  prad = prad_;
}

Reemission::~Reemission() {
}

Real Reemission::FixedWaveReemission(int b, int n, Real wave) {
  RadiationBand *pband = prad->my_bands(b);
  Real left = pband->spec[n].wave - pband->spec_bin_width/2.0;
  Real right = pband->spec[n].wave + pband->spec_bin_width/2.0;

  if (left > wave)
    return 0.;
  if (right < wave)
    return 0.;
  return 1.;
}

HydrogenReemission::HydrogenReemission(ParameterInput *pin, Radiation *prad_, Real ionization_energy_):
  Reemission(pin, prad_)
{
  ionization_energy = ionization_energy_;
  planck_constant = pin->GetReal("problem","planck_constant");
  speed_of_light = pin->GetReal("problem", "speed_of_light");
}

Real HydrogenReemission::ReemissionFunction(int b, int n, Real T, int k, int j, int i)
{
  Real critical_wave = (planck_constant * speed_of_light)/ionization_energy;

  return FixedWaveReemission(b, n, critical_wave);
}