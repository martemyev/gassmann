#include "gassmann.hpp"


void Gassmann::compute_substituted(const std::vector<double>& rho_array,
                                   const std::vector<double>& vp_array,
                                   const std::vector<double>& vs_array,
                                   const std::vector<double>& sat0_array,
                                   const std::vector<double>& sat1_array,
                                   const std::vector<double>& phi_array,
                                   std::vector<double>& rho_new_array,
                                   std::vector<double>& vp_new_array,
                                   std::vector<double>& vs_new_array) const
{
  const int n_values = rho_array.size();

  for (int i = 0; i < n_values; ++i)
  {
    const double phi = phi_array[i];
    const SeismicProperties old(rho_array[i], vp_array[i], vs_array[i]);

    const double old_saturation = sat0_array[i];
    double K_fl_mix   = SeismicProperties::bulk_modulus_fluid_mix(old_saturation, water.K, oil.K);
    double rho_fl_mix = SeismicProperties::density_fluid_mix(old_saturation, water.rho, oil.rho);

    const double Kstar = SeismicProperties::bulk_modulus_porous_rock_frame(old.K, mineral_matrix.K, K_fl_mix, phi);

    const double new_saturation = sat1_array[i];
    K_fl_mix   = SeismicProperties::bulk_modulus_fluid_mix(new_saturation, water.K, oil.K);
    rho_fl_mix = SeismicProperties::density_fluid_mix(new_saturation, water.rho, oil.rho);

    const double Ksat = SeismicProperties::bulk_modulus_saturated(Kstar, mineral_matrix.K, K_fl_mix, phi);

    rho_new_array[i] = SeismicProperties::bulk_density(grain.rho, rho_fl_mix, phi);
    vp_new_array[i]  = SeismicProperties::compres_velocity(Ksat, old.G, rho_new_array[i]);
    vs_new_array[i]  = SeismicProperties::shear_velocity(old.G, rho_new_array[i]);
  }
}


