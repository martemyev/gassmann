#ifndef SEISMIC_HPP
#define SEISMIC_HPP

#include <cmath>

double bulk_modulus(double rho, double vp, double vs);
double VoigtReussHill_averaging(double F1, double K1, double F2, double K2);
double shear_modulus(double rho, double vs);
double bulk_modulus_fluid_mix(double S_water, double K_water, double K_oil);
double density_fluid_mix(double S_water, double rho_water, double rho_oil);
double bulk_modulus_porous_rock_frame(double Ksat, double Ko, double Kfl, double phi);
double bulk_modulus_saturated(double Kstar, double Ko, double Kfl, double phi);
double bulk_density(double rho_grain, double rho_fluid, double phi);
double compres_velocity(double K, double G, double rho);
double shear_velocity(double G, double rho);


class SeismicProperties {
public:
  double rho; ///< Density
  double vp;  ///< P-wave velocity
  double vs;  ///< S-wave velocity
  double K;   ///< Bulk modulus
  double G;   ///< Shear modulus

  SeismicProperties(double rho_, double vp_, double vs_)
    : rho(rho_), vp(vp_), vs(vs_)
    , K(bulk_modulus(rho, vp, vs))
    , G(shear_modulus(rho, vs))
  {}
};


#endif // SEISMIC_HPP
