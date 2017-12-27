#include "Solver/Solver.h"
#include "Low_Mach_Solver/RK4_Low_Mach_Solver.h"
#include "Usefull_Headers/Variable_Vector_Isolator.h"
#include "Usefull_Headers/Gnuplot_Primitive_Variables.h"
#include "Eigen/Core"
#include "Eigen/Dense"
#include <cmath>
#include <chrono>

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
  int    number_of_cells = 600;
  double final_time = 0.01;
  int    frames = 2;

  double T_ignition = 1.0;
  double lambda = 0.0;
  double CFL = 0.5;
  double x_min = 0.0;
  double x_max;
  double lambda_max;
  double lambda_min;
  double lambda_run;
  int init_position;
  int position;
  // std::string filename = "../Movie/_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";
  // std::string filename = "../Movie/Test_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";
  // std::string filename = "../Movie/Case_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";
  global_solution_vector_type initial_solution;
  initial_solution.resize(number_of_cells);

  std::cout << "//////////////////////" << std::endl;
  std::cout << "Initial Conditions" << std::endl;
  std::cout << "//////////////////////" << std::endl;

  lambda_max = 96000;
  lambda_min = 93250;
  lambda_run = 90000;
  for (int i = 0; i < 200; ++i) {
  std::string filename = "Movie/Explicit_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";
  RK4_low_mach(lambda, number_of_cells, initial_solution, Le, Q, theta, T_ignition, gamma, x_max, mf);
  init_position = flame_position_algorithm(initial_solution, gamma);
  auto start = std::chrono::high_resolution_clock::now();
  std::cout << "//////////////////////" << std::endl;
  std::cout << "Solver, Lambda =" << lambda_run << std::endl;
  std::cout << "//////////////////////" << std::endl;
  auto solver = Solver<global_solution_vector_type>(initial_solution,Pr, Le, Q/(mf*mf*(gamma-1)), theta/(gamma*mf*mf), mf, lambda_run, gamma, number_of_cells, CFL, x_max, x_min, final_time, frames, filename);
  position = solver.solve();
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count() << " s\n";
  bisection_lambda(lambda_min, lambda_max, lambda_run, init_position, position);
}
};
