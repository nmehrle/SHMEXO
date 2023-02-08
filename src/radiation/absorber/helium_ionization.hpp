#ifndef HELIUM_IONIZATION_HPP
#define HELIUM_IONIZATION_HPP

// Athena++ header
#include "absorber.hpp"
#include "ionizing_absorber.hpp"

class HeliumIonization: public IonizingAbsorber {
public:
  HeliumIonization(RadiationBand *pband, int my_scalar_number, std::string name, ParameterInput *pin, std::string my_xc_file);
  ~HeliumIonization() {};

protected:
  // Helium Cross Section Constants from Verner1996
  // https://articles.adsabs.harvard.edu/pdf/1996ApJ...465..487V
  std::string xc_file;

  Real eV_conversion = 1.602E-12; // erg/eV
  Real mb_conversion = 1.E-18; // cm^2/Mb
  Real ry_conversion = 13.605693123; // eV/Ry

  void CalculateCrossSections(Spectrum const *spec, int nspec);
};

#endif