#include "gassmann.hpp"
#include "seismic.hpp"
#include "utilities.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>

const double RHO_WATER = 1000.0;
const double VP_WATER  = 1500.0;
const double VS_WATER  = 0.0;

const double RHO_OIL = 600.0;
const double VP_OIL  = 1200.0;
const double VS_OIL  = 0.0;

const double K_QUARTZ = 37e+9;
const double K_CLAY   = 15e+9; // wet clay
const double F_CLAY   = 0.15;  // fraction of clay
const double F_QUARTZ = 1.0 - F_CLAY;
const double RHO_QUARTZ = 2650.0;



void fill_up(int argc, char **argv, const std::string &option,
             std::vector<double> &array, std::string &filename) {
  const std::string file_option = option + "-file";
  const int n_cells = array.size();
  int pos;

  if ((pos = argcheck(argc, argv, file_option.c_str()))) {
    filename = argv[pos+1];
    read_binary(filename, n_cells, &array[0]);
    return;
  }

  if ((pos = argcheck(argc, argv, option.c_str()))) {
    const double value = atof(argv[pos+1]);
    array.clear();
    array.resize(n_cells, value);
    return;
  }

  throw std::runtime_error("No options: " + option + " nor " + file_option);
}



int main(int argc, char **argv) {
  if (argc == 1 || argcheck(argc, argv, "-h") || argcheck(argc, argv, "--help"))
  {
    std::cout << "\nUsage:\n" << argv[0] <<
                 " (-rho value or -rho-file file.bin)"
                 " (-vp value or -vp-file file.bin)"
                 " (-vs value or -vs-file file.bin)"
                 " (-phi value or -phi-file file.bin)"
                 " (-s0 value or -s0-file file.bin)"
                 " (-s1 value or -s1-file file.bin)"
                 " -nx n_cells_x"
                 " -ny n_cells_y"
                 " [-nz n_cells_z"
                 " -rho-out rho_out_file.bin"
                 " -vp-out vp_out_file.bin"
                 " -vs-out vs_out_file.bin]"
                 "\nwhere\n"
                 "rho_file     binary file with density in kg/m^3\n"
                 "vp_file      binary file with P-wave velocity in m/s\n"
                 "vs_file      binary file with S-wave velocity in m/s\n"
                 "s0_file      binary file with previous water saturation\n"
                 "s1_file      binary file with new water saturation\n"
                 "phi_file     binary file with porosity\n"
                 "phi_value    constant porosity\n"
                 "n_cells_x    number of cells along x-direction\n"
                 "n_cells_y    number of cells along y-direction\n"
                 "n_cells_z    number of cells along z-direction\n"
                 "rho_out_file binary file with new density\n"
                 "vp_out_file  binary file with new P-wave velocity\n"
                 "vs_out_file  binary file with new S-wave velocity\n"
                 "\n";
    std::cout << std::endl;
    return 1;
  }

  try {

    SeismicProperties water(RHO_WATER, VP_WATER, VS_WATER);
    SeismicProperties oil(RHO_OIL, VP_OIL, VS_OIL);
    double mineral_matrix_K = VoigtReussHill_averaging(F_CLAY, K_CLAY,
                                                       F_QUARTZ, K_QUARTZ);
    double grain_density = RHO_QUARTZ;

    std::cout << "Constants (assumptions):\n"
                 "Bulk modulus of the mineral matrix:\n"
                 "  K of quartz   " << K_QUARTZ*1e-9 << " GPa\n"
                 "  K of clay     " << K_CLAY*1e-9 << " GPa\n"
                 "  F of quartz   " << F_QUARTZ*100 << " %\n"
                 "  F of clay     " << F_CLAY*100 << " %\n"
                 "  K of matrix   " << mineral_matrix_K*1e-9 << " GPa\n"
                 "Grain density: " << grain_density << " kg/m^3\n";
    std::cout << std::endl;

    std::string rho_file, vp_file, vs_file, phi_file, sat0_file, sat1_file;
    std::string rho_out_file, vp_out_file, vs_out_file;
    int dim = 2; // default dimension
    int nx = 0, ny = 0, nz = 0;
    int pos;

    if ((pos = argcheck(argc, argv, "-nx")))
      nx = atoi(argv[pos+1]);
    else throw std::runtime_error("Number of cells in x-direction is not provided");

    if ((pos = argcheck(argc, argv, "-ny")))
      ny = atoi(argv[pos+1]);
    else throw std::runtime_error("Number of cells in y-direction is not provided");

    if ((pos = argcheck(argc, argv, "-nz"))) {
      nz = atoi(argv[pos+1]);
      dim = 3;
    }

    int n_cells = nx*ny;
    if (dim == 3) n_cells *= nz;

    std::vector<double> rho_array(n_cells);
    std::vector<double> vp_array(n_cells);
    std::vector<double> vs_array(n_cells);
    std::vector<double> phi_array(n_cells);
    std::vector<double> sat0_array(n_cells);
    std::vector<double> sat1_array(n_cells);

    fill_up(argc, argv, "-rho", rho_array,  rho_file);
    fill_up(argc, argv, "-vp",  vp_array,   vp_file);
    fill_up(argc, argv, "-vs",  vs_array,   vs_file);
    fill_up(argc, argv, "-phi", phi_array,  phi_file);
    fill_up(argc, argv, "-s0",  sat0_array, sat0_file);
    fill_up(argc, argv, "-s1",  sat1_array, sat1_file);

    if ((pos = argcheck(argc, argv, "-rho-out")))
      rho_out_file = argv[pos+1];
    else
      rho_out_file = file_stem(rho_file) + "_new.rho";

    if ((pos = argcheck(argc, argv, "-vp-out")))
      vp_out_file = argv[pos+1];
    else
      vp_out_file = file_stem(vp_file) + "_new.vp";

    if ((pos = argcheck(argc, argv, "-vs-out")))
      vs_out_file = argv[pos+1];
    else
      vs_out_file = file_stem(vs_file) + "_new.vs";

    std::vector<double> rho_new_array(n_cells);
    std::vector<double> vp_new_array(n_cells);
    std::vector<double> vs_new_array(n_cells);

    Gassmann gassmann(water, oil, mineral_matrix_K, grain_density);
    gassmann.compute_substituted(rho_array, vp_array, vs_array, sat0_array,
                                 sat1_array, phi_array, rho_new_array,
                                 vp_new_array, vs_new_array);

    write_binary(rho_out_file, n_cells, &rho_new_array[0]);
    write_binary(vp_out_file,  n_cells, &vp_new_array[0]);
    write_binary(vs_out_file,  n_cells, &vs_new_array[0]);

  } catch (const std::exception& e) {
    std::cout << "\n" << e.what() << "\n" << std::endl;
  } catch (...) {
    std::cout << "\nUnknown exception!\n" << std::endl;
  }

  return 0;
}

