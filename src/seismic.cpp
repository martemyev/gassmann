#include "seismic.hpp"

double bulk_modulus(double rho, double vp, double vs) {
  return rho * (vp*vp - 4./3.*vs*vs);
}

double VoigtReussHill_averaging(double F1, double K1, double F2, double K2) {
  const double K_Reuss = 1.0 / (F1/K1 + F2/K2);
  const double K_Voigt = F1*K1 + F2*K2;
  return (0.5 * (K_Voigt + K_Reuss));
}

double shear_modulus(double rho, double vs) {
  return rho*vs*vs;
}

double bulk_modulus_fluid_mix(double S_water, double K_water, double K_oil) {
  const double S_oil = 1.0 - S_water;
  return 1.0 / (S_water/K_water + S_oil/K_oil);
}

double density_fluid_mix(double S_water, double rho_water, double rho_oil) {
  const double S_oil = 1.0 - S_water;
  return (S_water*rho_water + S_oil*rho_oil);
}

double bulk_modulus_porous_rock_frame(double Ksat, double Ko, double Kfl, double phi) {
  double num = Ksat * (phi*Ko/Kfl + 1.0 - phi) - Ko;
  double den = phi*Ko/Kfl + Ksat/Ko - 1.0 - phi;
  return num / den;
}

double bulk_modulus_saturated(double Kstar, double Ko, double Kfl, double phi) {
  double num = (1.0 - Kstar/Ko);
  num *= num;
  double den = phi/Kfl + (1.0-phi)/Ko - Kstar/Ko/Ko;
  return Kstar + num/den;
}

double bulk_density(double rho_grain, double rho_fluid, double phi) {
  return (rho_grain*(1.0-phi) + rho_fluid*phi);
}

double compres_velocity(double K, double G, double rho) {
  return sqrt((K + 4./3.*G) / rho);
}

double shear_velocity(double G, double rho) {
  return sqrt(G / rho);
}
