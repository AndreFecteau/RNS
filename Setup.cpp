  #include "Solver/Solver.h"
#include "Low_Mach_Solver/RK4_Low_Mach_Solver.h"
#include "Usefull_Headers/Variable_Vector_Isolator.h"
#include "Gnuplot_RNS/Gnuplot_Primitive_Variables.h"
#include "Solver/Implicit_Marching_Generic.h"
#include "Solver/Explicit_Marching.h"
#include "Eigen/Core"
#include "Eigen/Dense"
#include <cmath>

typedef Eigen::Matrix<double, 4, 1> Vector_type;
using matrix_type = Eigen::Matrix<double, 4,4>;
using global_solution_vector_type = std::vector<Vector_type>;
using solution_vector_type = typename global_solution_vector_type::value_type;
#include "Usefull_Headers/Initial_Conditions.h"

using implicit_marching_type = Implicit_Marching<global_solution_vector_type, matrix_type>;
using explicit_marching_type = Explicit_Marching<global_solution_vector_type, matrix_type>;



void bisection_lambda(double& lambda_min, double& lambda_max, double& lambda_run, bool check) {
  if (check == 1){
    lambda_min = lambda_run;
    lambda_run = (lambda_min + lambda_max)*0.5;
  } else {
    lambda_max = lambda_run;
    lambda_run = (lambda_min + lambda_max)*0.5;
  }
}

int main(){
  // double Pr = 0.75;
  // double Le = 1.0;
  // double Q = 4;
  // double theta = 30;

  double Pr = 0.75;
  double Le = 0.3;
  double mf = 0.005;
  double gamma = 1.4;
  double Q_low_mach = 9.0;
  double T_ignition = 1.0;
  double theta_low_mach =500.0/9.0;
  double Q = Q_low_mach/(mf*mf*(gamma-1));
  double theta =theta_low_mach/(gamma*mf*mf);

  int    number_of_cells = 500;
  double frame_time = 1;

  double lambda = 0.0;
  double x_min = 0.0;
  double x_max;
  double lambda_max;
  double lambda_min;
  double lambda_run;
  double target_residual = 1e-16;

  double Theta = 1.0;
  double zeta = 0.5;
  double CFL =  2;

  // std::string filename = "Movie/Plot_Euler_" + tostring(frame_time) + "_" + tostring(number_of_cells) + "_";
  std::string filename = "Movie/Test2_" + tostring(number_of_cells) + "_";

  global_solution_vector_type initial_solution;
  initial_solution.resize(number_of_cells);

  std::cout << "//////////////////////" << std::endl;
  std::cout << "Initial Conditions" << std::endl;
  std::cout << "//////////////////////" << std::endl;

  lambda_max = 100000;
  lambda_min = 90000;
  lambda_run = 95290;

  // while(1>0) {
  // case_4(frame_time, number_of_cells, initial_solution, gamma, x_max, x_min);
  RK4_low_mach_initial_conditions(lambda, number_of_cells, initial_solution, Le, Q_low_mach,
               theta_low_mach, T_ignition, gamma, x_max, mf);
  auto explicit_march = explicit_marching_type(Pr, Le, Q, theta, mf, gamma,
                        number_of_cells, CFL, (x_max - x_min)/number_of_cells);
  auto implicit_march = implicit_marching_type(Pr, Le, Q, theta, mf, gamma,
                        number_of_cells, CFL, (x_max - x_min)/number_of_cells, Theta, zeta);

  plot<global_solution_vector_type>(filename + std::to_string(static_cast<int>(lambda_run)) + "_0",
                                    initial_solution, (x_max - x_min)/number_of_cells);

  auto solver = Solver<global_solution_vector_type, matrix_type>(initial_solution, lambda_run,
                                                                 filename);
  // bool check = solver.solve<implicit_marching_type>(implicit_march, target_residual,
  //                                                           frame_time, gamma);
  solver.solve<implicit_marching_type>(implicit_march, target_residual, frame_time, gamma);
  // solver.solve<explicit_marching_type>(explicit_march, target_residual, frame_time, gamma);
  // solver.solve<explicit_euler_marching_type>(explicit_euler_march, target_residual, frame_time, gamma);

  // bisection_lambda(lambda_min, lambda_max, lambda_run, check);
// }

};
