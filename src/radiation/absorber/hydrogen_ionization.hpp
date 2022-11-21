#ifndef HYDROGEN_IONIZATION_HPP
#define HYDROGEN_IONIZATION_HPP

// Athena++ header
#include "absorber.hpp"

class HydrogenIonization: public Absorber {
public:
  HydrogenIonization(RadiationBand *pband, int my_scalar_number);
  ~HydrogenIonization() {};

  void CalculateAsorptionCoefficient(AthenaArray<Real> const& prim, int n, int k, int j, int i);

  Real ionization_energy=2.1798723611E-18; // Joule
protected:
  // cross section constant
  Real A0=6.30431812E-22; // m2
  // rydberg frequency
  Real nu_0=3.2898419603E15; // s-1
  Real c=2.998E8; // m/s

  void CalculateCrossSections(Spectrum const& spec, int nspec);
  void CalculateEnergyFunctions(Spectrum const& spec, int nspec);
};

#endif