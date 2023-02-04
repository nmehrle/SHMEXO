#ifndef HELIUM_IONIZATION_HPP
#define HELIUM_IONIZATION_HPP

// Athena++ header
#include "absorber.hpp"
#include "ionizing_absorber.hpp"

class HeliumIonization: public IonizingAbsorber {
public:
  HeliumIonization(RadiationBand *pband, int my_scalar_number, std::string name, ParameterInput *pin);
  ~HeliumIonization() {};

protected:
  // Helium Cross Section Constants from Verner1996
  // https://articles.adsabs.harvard.edu/pdf/1996ApJ...465..487V

  Real Eth  = 2.459E1; // eV
  Real Emax = 5.000E4; // eV
  Real E0   = 1.361E1; // eV
  Real sig0 = 9.492E2; // Mb

  Real eV_conversion = 1.602E-19; // J/eV
  Real mb_conversion = 1.E-18; // m^2/Mb

  // following are dimensionless
  Real ya = 1.469;
  Real P  = 3.188;
  Real yw = 2.039;
  Real y0 = 4.434E-1;
  Real y1 = 2.136;

  void CalculateCrossSections(Spectrum const *spec, int nspec);
};

#endif