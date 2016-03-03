#ifndef GASSMANN_HPP
#define GASSMANN_HPP

#include "seismic.hpp"

#include <vector>

class Gassmann {
public:
  Gassmann(SeismicProperties water_, SeismicProperties oil_,
           SeismicProperties mineral_matrix_, SeismicProperties grain_)
    : water(water_), oil(oil_), mineral_matrix(mineral_matrix_), grain(grain_)
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
  SeismicProperties mineral_matrix;
  SeismicProperties grain;
};


#endif // GASSMANN_HPP
