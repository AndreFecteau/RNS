#ifndef GNUPLOT_PRIMITIVE_VARIABLES_REDUCED_H
#define GNUPLOT_PRIMITIVE_VARIABLES_REDUCED_H

#include <fstream>

#include "../Usefull_Headers/Variable_Vector_Isolator.h"

template <typename grid_type, typename flow_type>
void plot_reduced(std::string file_name, grid_type grid, flow_type flow , int number_of_cell_in_output){

  std::ofstream gnu_input_file;
  gnu_input_file.open (file_name + ".dat");
  gnu_input_file.precision(40);

  gnu_input_file << "#rho u p T Y" << std::endl;
  for (size_t i = 0; i < grid.number_of_cells(); i++) {
    if(i%(grid.number_of_cells()/number_of_cell_in_output) == 0 || i == grid.number_of_cells()-1){
    auto temp = Variable_Vector_Isolator<grid_type>(grid.global_solution_vector[i], 1.4);
    // double mach = temp.u()/sqrt(1.4*temp.p()/temp.rho())
    gnu_input_file << (i+0.5) * grid.dx() << " " << temp.rho() << " " << temp.u()  << " " << temp.p() << " " << temp.T() << " " << temp.Y() << " " << temp.u()/sqrt(1.4*temp.p()/temp.rho())
                   << " " << flow.lambda * temp.rho() * temp.Y() * exp(-flow.theta()/temp.T())  << std::endl;
    }
  }
  gnu_input_file.close();
}

#endif //#ifndef GNUPLOT_PRIMITIVE_VARIABLES_REDUCED_H
