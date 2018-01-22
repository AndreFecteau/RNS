#ifndef GNUPLOT_PRIMITIVE_VARIABLES_H
#define GNUPLOT_PRIMITIVE_VARIABLES_H

#include <iostream>
#include <fstream>
#include <cstdarg>
#include <vector>
#include "Variable_Vector_Isolator.h"

template <typename global_solution_vector_type>
void plot(std::string file_name, global_solution_vector_type solution_vector, double dx){
  using solution_vector_type = typename global_solution_vector_type::value_type;
  std::ofstream gnu_input_file;
  gnu_input_file.open (file_name + ".dat");
  gnu_input_file << "#rho u p T Y" << std::endl;
  for (size_t i = 0; i < solution_vector.size(); i++) {
    Variable_Vector_Isolator<solution_vector_type> temp = Variable_Vector_Isolator<solution_vector_type>(solution_vector[i], 1.4);
    gnu_input_file << i * dx << " " << temp.rho() << " " << temp.u()  << " " << temp.p() << " " << temp.T() << " " << temp.Y() << std::endl;
  }
  gnu_input_file.close();
}

#endif //#ifndef GNUPLOT_RNS_PRIMITIVE_VARIABLES_H
