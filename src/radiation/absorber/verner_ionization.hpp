#ifndef VERNER_IONIZATION_HPP
#define VERNER_IONIZATION_HPP

// Photoionization cross sections from Verner+ 1996b
// Atomic Data for Astrophysics. II. New Analytic FITS for Photoionization Cross Sections of Atoms and Ions
// https://ui.adsabs.harvard.edu/abs/1996ApJ...465..487V/abstract
// Athena++ header
#include "absorber.hpp"
#include "ionizing_absorber.hpp"

class VernerIonization: public IonizingAbsorber {
public:
  VernerIonization(RadiationBand *pband, int my_scalar_number, std::string name, ParameterInput *pin, Real Eth_in, Real Emax_in, Real E0_in, Real sigma0_in, Real ya_in, Real P_in, Real yw_in, Real y0_in, Real y1_in);
  ~VernerIonization() {};

protected:
  Real Eth, Emax, E0, sigma0, ya, P, yw, y0, y1;


  void CalculateCrossSections(Spectrum const *spec, int nspec);
};

#endif