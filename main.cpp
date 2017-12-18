#include "Solver/Solver.h"
#include "Low_Mach_Solver/RK4_Low_Mach_Solver.h"
#include "Usefull_Headers/Variable_Vector_Isolator.h"
#include "Usefull_Headers/Gnuplot_Primitive_Variables.h"
#include "Eigen/Core"
#include "Eigen/Dense"
#include <cmath>
#include <chrono>
#comment to check
typedef Eigen::Matrix<double, 5, 1> Vector5d;
using global_solution_vector_type = std::vector<Vector5d>;
using solution_vector_type = typename global_solution_vector_type::value_type;
#include "Usefull_Headers/Initial_Conditions.h"

template <typename T>
std::string tostring(T name) {
  return std::to_string(static_cast<int>(name));
}

void bisection_lambda(double& lambda_min, double& lambda_max, double& lambda_run, double& T_sum_initial, double& T_sum) {
  if (T_sum < T_sum_initial){
    lambda_min = lambda_run;
    lambda_run = (lambda_min + lambda_max)*0.5;
  } else {
    lambda_max = lambda_run;
    lambda_run = (lambda_min + lambda_max)*0.5;
  }
}

double flame_position_algorithm(global_solution_vector_type global_solution_vector, double gamma) {
  double sum = 0.0;
  for (size_t i = 0; i < global_solution_vector.size(); ++i){
    Variable_Vector_Isolator<solution_vector_type> var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i], gamma);
    sum += var_vec.T();
  }
  return sum;
}


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
  double final_time = 2.0;
  int    frames = 2;

  double T_ignition = 1.0;
  double lambda = 0.0;
  double CFL = 0.5;
  double x_min = 0.0;
  double x_max;
  double lambda_max;
  double lambda_min;
  double lambda_run;
  double T_sum_initial;
  double T_sum;
  // std::string filename = "../Movie/_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";
  // std::string filename = "../Movie/Test_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";
  // std::string filename = "../Movie/Case_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";
  global_solution_vector_type initial_solution;
  initial_solution.resize(number_of_cells);

  std::cout << "//////////////////////" << std::endl;
  std::cout << "Initial Conditions" << std::endl;
  std::cout << "//////////////////////" << std::endl;

  T_sum_initial = flame_position_algorithm(initial_solution, gamma);
  lambda_max = 94350;
  lambda_min = 92500;
  lambda_run = 93250;
  for (int i = 0; i < 200; ++i) {
  std::string filename = "Movie/Refinement_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";
  RK4_low_mach(lambda, number_of_cells, initial_solution, Le, Q, theta, T_ignition, gamma, x_max, mf);
  auto start = std::chrono::high_resolution_clock::now();
  std::cout << "//////////////////////" << std::endl;
  std::cout << "Solver, Lambda =" << lambda_run << std::endl;
  std::cout << "//////////////////////" << std::endl;
  auto solver = Solver<global_solution_vector_type>(initial_solution,Pr, Le, Q/(mf*mf*(gamma-1)), theta/(gamma*mf*mf), mf, lambda_run, gamma, number_of_cells, CFL, x_max, x_min, final_time, frames, filename);
  T_sum = solver.solve();
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count() << " s\n";
  bisection_lambda(lambda_min, lambda_max, lambda_run, T_sum_initial, T_sum);
}
};
