#include "Solver.h"
#include "RK4_Solver_Incompressible_Reactive.h"
#include "Variable_Vector_Isolator.h"
#include "Gnuplot_Primitive_Variables.h"
#include "Eigen/Core"
#include "Eigen/Dense"
#include <cmath>
#include <chrono>

typedef Eigen::Matrix<double, 5, 1> Vector5d;
using global_solution_vector_type = std::vector<Vector5d>;
using solution_vector_type = typename global_solution_vector_type::value_type;
#include "Initial_Conditions.h"

template <typename T>
std::string tostring(T name) {
  return std::to_string(static_cast<int>(name));
}

///////////////////////////////////////////////////////////////////////////////
// Bisection for flame speed
///////////////////////////////////////////////////////////////////////////////
// inline void bisection_lambda(double& lambda_min, double& lambda_max, double& lambda_curr, double check) {
//   if (check == 0){
//     lambda_min = lambda_curr;
//     lambda_curr = lambda_curr + (lambda_max-lambda_curr) / 2;
//   } else if (check == 1){
//     lambda_max = lambda_curr;
//     lambda_curr = lambda_curr - (lambda_curr - lambda_min)/2;
//   }
// }

double bisection_lambda(double& lambda_min, double& lambda_max, double& lambda_high, double& lambda_low, double& p_high, double& p_low, bool& check) {
  double lambda_run;
  if (fabs(p_high-28571.4) > fabs(p_low-28571.4)){
    lambda_max = lambda_high;
    lambda_high = lambda_low;
    p_high = p_low;
    lambda_run = (lambda_low + lambda_min) * 0.5;
    lambda_low = lambda_run;
    check = 1;
  } else {
    lambda_min = lambda_low;
    lambda_low = lambda_high;
    p_low = p_high;
    lambda_run = (lambda_high + lambda_max) * 0.5;
    lambda_high = lambda_run;
    check = 0;
  }
  return lambda_run;
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
  int    number_of_cells = 1500;
  double final_time = 2000.0;
  int    frames = 100;

  double T_ignition = 1.0;
  double lambda = 0.0;
  double CFL = 0.5;
  double x_min = 0.0;
  double x_max;
  double lambda_max;
  double lambda_min;
  double lambda_high = 0;
  double lambda_low;
  double lambda_run;
  double p_low = 0;
  double p_high = 0;
  bool check = 0;
  // std::string filename = "../Movie/Test_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";
  // std::string filename = "../Movie/Refinement_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";
  // std::string filename = "../Movie/Case_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";
  global_solution_vector_type initial_solution;
  initial_solution.resize(number_of_cells);

  std::cout << "//////////////////////" << std::endl;
  std::cout << "Initial Conditions" << std::endl;
  std::cout << "//////////////////////" << std::endl;

  // case_4(final_time, number_of_cells, initial_solution, gamma, x_max, x_min);
  // RK4_Incompressible(lambda, number_of_cells, initial_solution, Le, Q, theta, T_ignition, gamma, x_max, mf);
  // lambda_max = 1.1 * lambda;
  // lambda_min = 0.9 * lambda;
  // lambda_run = 1.05 * lambda;
  // lambda_high = lambda_run;
  // lambda_low = 1.0 * lambda;
  std::string filename = "../Movie/itt_" + tostring(final_time) + "_" + tostring(number_of_cells) + "_";
  // auto solver = Solver<global_solution_vector_type>(initial_solution,Pr, Le, Q/(mf*mf*(gamma-1)), theta/(gamma*mf*mf), mf, lambda_low, gamma, number_of_cells, CFL, x_max, x_min, final_time, frames, filename);
  // p_low = solver.solve();
  // for (int i = 0; i < 200; ++i) {
  RK4_Incompressible(lambda, number_of_cells, initial_solution, Le, Q, theta, T_ignition, gamma, x_max, mf);
  lambda = 91000;
  auto start = std::chrono::high_resolution_clock::now();
  std::cout << "//////////////////////" << std::endl;
  std::cout << "Solver, Lambda =" << lambda_run << std::endl;
  std::cout << "//////////////////////" << std::endl;
  auto solver = Solver<global_solution_vector_type>(initial_solution,Pr, Le, Q/(mf*mf*(gamma-1)), theta/(gamma*mf*mf), mf, lambda, gamma, number_of_cells, CFL, x_max, x_min, final_time, frames, filename);
  // if (check == 1) {
    p_low = solver.solve();
  // } else {
  //   p_high = solver.solve();
  // }
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count() << " s\n";
  // std::cout << lambda_max << " : " << lambda_min << " : " << lambda_high << " : " << lambda_low << " : " << p_high << " : " << p_low << " : " << lambda_run << std::endl;
  // lambda_run = bisection_lambda(lambda_min, lambda_max, lambda_high, lambda_low, p_high, p_low, check);
  // std::cout << lambda_max << " : " << lambda_min << " : " << lambda_high << " : " << lambda_low << " : " << p_high << " : " << p_low << " : " << lambda_run << std::endl;
// }
};
