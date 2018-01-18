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
  // double Pr = 0.75;
  // double Le = 1.0;
  // double Q = 4;
  // double theta = 30;

  double Pr = 0.75;
  double Le = 0.3;
  double Q = 9.0;
  double theta =500.0/9.0;
  double gamma = 1.4;
  double mf = 0.005;
  int    number_of_cells = 6000;
  double final_time = 100;
  int    frames = 100;

  double T_ignition = 1.0;
  double lambda = 0.0;
  double CFL = 0.5;
  double x_min = 0.0;
  double x_max;
  double lambda_max;
  double lambda_min;
  double lambda_run;
  // std::string filename = "../Movie/_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";
  // std::string filename = "../Movie/Test_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";
  // std::string filename = "../Movie/Case_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";
  global_solution_vector_type initial_solution;
  initial_solution.resize(number_of_cells);

  std::cout << "//////////////////////" << std::endl;
  std::cout << "Initial Conditions" << std::endl;
  std::cout << "//////////////////////" << std::endl;


  lambda_max = 100000;
  lambda_min = 70000;
  lambda_run = 95287;
  std::string filename = "Movie/Refinement9_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";
  RK4_low_mach(lambda, number_of_cells, initial_solution, Le, Q, theta, T_ignition, gamma, x_max, mf);
  auto solver = Solver<global_solution_vector_type>(initial_solution,Pr, Le, Q/(mf*mf*(gamma-1)), theta/(gamma*mf*mf), mf, lambda_run, lambda_min, lambda_max, gamma, number_of_cells, CFL, x_max, x_min, final_time, frames, filename);
  solver.solve();

};
