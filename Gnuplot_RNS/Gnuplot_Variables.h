#ifndef GNUPLOT_VARIABLES_H
#define GNUPLOT_VARIABLES_H

#include <iostream>
#include <fstream>
#include <cstdarg>
#include <vector>
#include "../Usefull_Headers/Variable_Vector_Isolator.h"

template <typename global_solution_vector_type>
void plot_v(std::string file_name, global_solution_vector_type solution_vector, double dx){
  std::ofstream gnu_input_file;
  gnu_input_file.open (file_name + ".dat");
  gnu_input_file << "#rho u E Y" << std::endl;
  for (size_t i = 0; i < solution_vector.size(); i++) {
    gnu_input_file << i * dx << " " << solution_vector[i][0] << " " << solution_vector[i][1]  << " " << solution_vector[i][2]  << " " << solution_vector[i][3] << std::endl;
  }
  gnu_input_file.close();
}

#endif //#ifndef #define GNUPLOT_VARIABLES_H
