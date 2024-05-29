#ifndef FILE_READ_IONIZATION_HPP
#define FILE_READ_IONIZATION_HPP

// Athena++ header
#include "absorber.hpp"
#include "ionizing_absorber.hpp"

class FileReadIonization: public IonizingAbsorber {
public:
  FileReadIonization(RadiationBand *pband, std::string name, int my_scalar_number, int my_ion_number, ParameterInput *pin, std::string my_xc_file);
  ~FileReadIonization() {};

protected:
  std::string xc_file;

  void CalculateCrossSections(Spectrum const *spec, int nspec);
};

#endif