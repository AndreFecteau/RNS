#define HYPERBOLIC
#define VISCOUS
#define SOURCE
// #define EXPLICIT
// #define IMPLICIT

#include "Solver/Solver.h"
#include "Low_Mach_Solver/RK4_Low_Mach_Solver.h"
#include "Usefull_Headers/Variable_Vector_Isolator.h"
#include "Gnuplot_RNS/Gnuplot_Primitive_Variables.h"
#include "Solver/Implicit_Marching_Generic.h"
#include "Solver/Explicit_Marching.h"
#include "Eigen/Core"
#include "Eigen/Dense"
#include <cmath>


typedef Eigen::Matrix<long double, 4, 1> Vector_type;
using matrix_type = Eigen::Matrix<long double, 4,4>;
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

  // int    number_of_cells =400000;
  int    number_of_cells;
  double frame_time = 1e1;

  double lambda = 0.0;
  double x_min = 0.0;
  double x_max;
  double T_ignition = 1.0;
  double lambda_max;
  double lambda_min;
  double lambda_run;
  double target_residual = 1e-17;

  double Theta = 1.0;
  double zeta = 0.0;
  double CFL =  1e6;
  double per_FL = 512.0;
  double dx = 1.0/per_FL;
  double domaine_length = 500;
  std::ofstream gnu_input_file;
  gnu_input_file.open("Convergence_Plot.dat", std::ios_base::app);
  gnu_input_file << "#number_of_cells residual time" << std::endl;

  // initial_solution.resize(number_of_cells);
  // manufactured_solution(number_of_cells, initial_solution, x_max, x_min, dx);


  lambda_max = 95600;
  lambda_min = 94500;
  lambda_run = 94550;


  while(CFL < 1e9) {
    // double dx = 1.0/per_FL;
    global_solution_vector_type initial_solution;
    RK4_low_mach_initial_conditions(lambda, number_of_cells, initial_solution, Le, Q_low_mach,
      theta_low_mach, T_ignition, gamma, x_max, mf, dx, domaine_length);
    // double CFL =  number_of_cells/1000;
  // std::string filename = "Movie/Plot_Euler_" + tostring(frame_time) + "_" + tostring(number_of_cells) + "_";
  // std::string filename = "Movie/Test_Implicit_Residual_" + tostring(number_of_cells) + "_";
  // manufactured_solution(number_of_cells, initial_solution, x_max, x_min, dx);

  std::string filename = "Movie/Plot_HLLE_8_" + tostring(per_FL) + "_"
                                       + tostring(domaine_length) + "_"
                                       + tostring(log10(CFL)) + "_"
                                       + tostring(lambda_run) + "_";
  // std::string filename = "Movie/Exact_" + tostring(number_of_cells) + "_";


  std::cout << "//////////////////////" << std::endl;
  std::cout << "Initial Conditions "  << number_of_cells << std::endl;
  std::cout << "//////////////////////" << std::endl;


  // straight_line(number_of_cells, initial_solution, x_max, x_min, mf, gamma);
  // manufactured_solution(number_of_cells, initial_solution, x_max, x_min, dx);
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
  // CFL /= 10;
  bisection_lambda(lambda_min, lambda_max, lambda, check);
//   domaine_length += 2000;
  // per_FL *=2.0;
}

};
