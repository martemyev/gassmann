#ifndef GASSMANN_HPP
#define GASSMANN_HPP

#include "seismic.hpp"

#include <vector>

class Gassmann {
public:
  Gassmann(SeismicProperties water_, SeismicProperties oil_,
           double mineral_matrix_K_, double grain_density_)
    : water(water_), oil(oil_), mineral_matrix_K(mineral_matrix_K_),
      grain_density(grain_density_)
  {}

  void compute_substituted(const std::vector<double>& rho_array,
                           const std::vector<double>& vp_array,
                           const std::vector<double>& vs_array,
                           const std::vector<double>& sat0_array,
                           const std::vector<double>& sat1_array,
                           const std::vector<double>& phi_array,
                           std::vector<double>& rho_new_array,
                           std::vector<double>& vp_new_array,
                           std::vector<double>& vs_new_array) const;

private:
  SeismicProperties water;
  SeismicProperties oil;
  double mineral_matrix_K;
  double grain_density;
};


#endif // GASSMANN_HPP
