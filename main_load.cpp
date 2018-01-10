#include "Solver/Solver.h"
#include "Low_Mach_Solver/RK4_Low_Mach_Solver.h"
#include "Usefull_Headers/Variable_Vector_Isolator.h"
#include "Usefull_Headers/Gnuplot_Primitive_Variables.h"
#include "Eigen/Core"
#include "Eigen/Dense"
#include <cmath>

typedef Eigen::Matrix<double, 5, 1> Vector5d;
using global_solution_vector_type = std::vector<Vector5d>;
using solution_vector_type = typename global_solution_vector_type::value_type;
#include "Usefull_Headers/Initial_Conditions.h"

int main(){
  double final_time = 100;
  int number_of_cells = 6000;
  std::string filename = "Movie/Refinement1_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";

  auto solver = Solver<global_solution_vector_type>();
  unserialize_to_file(solver, filename);
  solver.solve();
};
