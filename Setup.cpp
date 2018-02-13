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
  double theta_low_mach =500.0/9.0;
  double Q = Q_low_mach/(mf*mf*(gamma-1));
  double theta =theta_low_mach/(gamma*mf*mf);

  int    number_of_cells = 2000;
  double frame_time = 1e1;
  double T_ignition = 1.0;
  double theta_low_mach =5.0*(1.0+Q_low_mach)*(1.0+Q_low_mach) / Q_low_mach;
  double Q = Q_low_mach/(mf*mf*(gamma-1));
  double theta =theta_low_mach/(gamma*mf*mf);

  int    number_of_cells =100000;
  double frame_time = 1e3;

  double lambda = 0.0;
  double x_min = 0.0;
  double x_max;
  double T_ignition = 1.0;
  double lambda_max;
  double lambda_min;
  double lambda_run;
  double target_residual = 1e-14;

  double Theta = 1.0;
  double zeta = 0.0;
  double CFL =  1e6;
  std::ofstream gnu_input_file;
  gnu_input_file.open("Convergence_Plot.dat", std::ios_base::app);
  gnu_input_file << "#number_of_cells residual time" << std::endl;

  global_solution_vector_type initial_solution;
  initial_solution.resize(number_of_cells);
  RK4_low_mach_initial_conditions(lambda, number_of_cells, initial_solution, Le, Q_low_mach,
    theta_low_mach, T_ignition, gamma, x_max, mf);
  while(CFL > 100) {
    // double CFL =  number_of_cells/1000;
  // std::string filename = "Movie/Plot_Euler_" + tostring(frame_time) + "_" + tostring(number_of_cells) + "_";
  // std::string filename = "Movie/Test_Implicit_Residual_" + tostring(number_of_cells) + "_";
  std::string filename = "Movie/Test_Explicit_Residual_" + tostring(number_of_cells) + "_" + tostring(CFL) + "_";
  // std::string filename = "Movie/Exact_" + tostring(number_of_cells) + "_";


  std::cout << "//////////////////////" << std::endl;
  std::cout << "Initial Conditions"  << std::endl;
  std::cout << "//////////////////////" << std::endl;

  lambda_max = 95700;
  lambda_min = 94000;
  lambda_run = 95400;

  // straight_line(number_of_cells, initial_solution, x_max, x_min, mf, gamma);
  // manufactured_solution(number_of_cells, initial_solution, x_max, x_min);
  // case_4(frame_time, number_of_cells, initial_solution, gamma, x_max, x_min);
  // RK4_low_mach_initial_conditions(lambda, number_of_cells, initial_solution, Le, Q_low_mach,
  //              theta_low_mach, T_ignition, gamma, x_max, mf);
  auto explicit_march = explicit_marching_type(Pr, Le, Q, theta, mf, gamma,
                        number_of_cells, CFL, (x_max - x_min)/number_of_cells);
  auto implicit_march = implicit_marching_type(Pr, Le, Q, theta, mf, gamma,
                        number_of_cells, CFL, (x_max - x_min)/number_of_cells, Theta, zeta);

  plot<global_solution_vector_type>(filename+"0",
                                    initial_solution, (x_max - x_min)/number_of_cells);

  auto solver = Solver<global_solution_vector_type, matrix_type>(initial_solution, lambda,
                                                                 filename);

  bool check = solver.solve<implicit_marching_type>(implicit_march, target_residual, frame_time, gamma);
  // bool check = solver.solve<explicit_marching_type>(explicit_march, target_residual, frame_time, gamma);
  // solver.solve<explicit_euler_marching_type>(explicit_euler_march, target_residual, frame_time, gamma);
  CFL /= 10;
  // bisection_lambda(lambda_min, lambda_max, lambda_run, check);
}

};
