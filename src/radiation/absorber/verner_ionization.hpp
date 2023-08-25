#ifndef VERNER_IONIZATION_HPP
#define VERNER_IONIZATION_HPP

// Photoionization cross sections from Verner+ 1996b
// Atomic Data for Astrophysics. II. New Analytic FITS for Photoionization Cross Sections of Atoms and Ions
// https://ui.adsabs.harvard.edu/abs/1996ApJ...465..487V/abstract
// Athena++ header
#include "absorber.hpp"
#include "ionizing_absorber.hpp"

struct VernerIonizationParams {
  Real Eth, Emax, E0, sigma0, ya, P, yw, y0, y1;

  VernerIonizationParams(Real Eth_, Real Emax_, Real E0_, Real sigma0_, Real ya_, Real P_, Real yw_, Real y0_, Real y1_):
    Eth(Eth_), Emax(Emax_), E0(E0_), sigma0(sigma0_), ya(ya_), P(P_), yw(yw_), y0(y0_), y1(y1_)
  {
  }
};

class VernerIonization: public IonizingAbsorber {
public:
  VernerIonization(RadiationBand *pband, std::string name, int my_scalar_number, int my_ion_number, ParameterInput *pin, VernerIonizationParams *params);
  ~VernerIonization() {};

protected:
  Real Eth, Emax, E0, sigma0, ya, P, yw, y0, y1;


  void CalculateCrossSections(Spectrum const *spec, int nspec);
};

#endif