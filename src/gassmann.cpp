#include "gassmann.hpp"
#include "utilities.hpp"

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

  double phi, K_fl_mix, rho_fl_mix, old_saturation, new_saturation;
  std::vector<double> Kstar(n_values);
  std::vector<double> Ksat(n_values);

  for (int i = 0; i < n_values; ++i)
  {
    phi = phi_array[i];
    const SeismicProperties old(rho_array[i], vp_array[i], vs_array[i]);

    old_saturation = sat0_array[i];
    K_fl_mix   = bulk_modulus_fluid_mix(old_saturation, water.K, oil.K);
    rho_fl_mix = density_fluid_mix(old_saturation, water.rho, oil.rho);

    Kstar[i] = mineral_matrix_K; //bulk_modulus_porous_rock_frame(old.K, mineral_matrix_K, K_fl_mix, phi);

    new_saturation = sat1_array[i];
    K_fl_mix   = bulk_modulus_fluid_mix(new_saturation, water.K, oil.K);
    rho_fl_mix = density_fluid_mix(new_saturation, water.rho, oil.rho);

    Ksat[i] = bulk_modulus_saturated(Kstar[i], mineral_matrix_K, K_fl_mix, phi);

    rho_new_array[i] = bulk_density(grain_density, rho_fl_mix, phi);
    vp_new_array[i]  = compres_velocity(Ksat[i], old.G, rho_new_array[i]);
    vs_new_array[i]  = shear_velocity(old.G, rho_new_array[i]);
  }

  write_binary("Kstar.bin", n_values, &Kstar[0]);
  write_binary("Ksat.bin",  n_values, &Ksat[0]);
}


