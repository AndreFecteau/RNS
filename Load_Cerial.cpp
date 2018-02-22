#define HYPERBOLIC
#define VISCOUS
#define SOURCE

#include "Solver/Solver.h"
#include "Low_Mach_Solver/RK4_Low_Mach_Solver.h"
#include "Usefull_Headers/Variable_Vector_Isolator.h"
#include "Gnuplot_RNS/Gnuplot_Primitive_Variables.h"
#include "Solver/Implicit_Marching.h"
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

int main(){
  // double Pr = 0.75;
  // double Le = 0.3;
  // double mf = 0.005;
  // double gamma = 1.4;
  // double Q_low_mach = 9.0;
  // double theta_low_mach =500.0/9.0;
  // double Q = Q_low_mach/(mf*mf*(gamma-1));
  // double theta =theta_low_mach/(gamma*mf*mf);
  //
  // // int    number_of_cells =400000;
  // int    number_of_cells;
  // double frame_time = 1e1;
  //
  // double lambda = 0.0;
  // double x_min = 0.0;
  // double x_max;
  // double T_ignition = 1.0;
  // double lambda_max;
  // double lambda_min;
  // double lambda_run;
  // double target_residual = 1e-17;
  //
  // double Theta = 1.0;
  // double zeta = 0.0;
  // double CFL =  1e6;
  // double per_FL = 512.0;
  // double dx = 1.0/per_FL;
  // double domaine_length = 500;
  //
 //  std::string filename = "Movie/Plot_CD_8_" + tostring(per_FL) + "_"
 //                                       + tostring(domaine_length) + "_"
 //                                       + tostring(log10(CFL)) + "_";
 //    global_solution_vector_type initial_solution;
 //    RK4_low_mach_initial_conditions(lambda, number_of_cells, initial_solution, Le, Q_low_mach,
 //      theta_low_mach, T_ignition, gamma, x_max, mf, dx, domaine_length);
 // // initial_solution.resize(number_of_cells);
 // std::cout << number_of_cells << std::endl;
 //  auto solver = Solver<global_solution_vector_type, matrix_type>();
 //  unserialize_to_file(solver, filename);
 //  // // auto explicit_march = explicit_marching_type(Pr, Le, Q, theta, mf, gamma,
 //  // //                       number_of_cells, CFL, (x_max - x_min)/number_of_cells);
 //  auto implicit_march = implicit_marching_type(Pr, Le, Q, theta, mf, gamma,
 //                        number_of_cells, CFL, (x_max - x_min)/number_of_cells, Theta, zeta);
 //  // solver.expand_solution_vector<implicit_marching_type>(1,implicit_march);
 //  solver.solve<implicit_marching_type>(implicit_march, target_residual, frame_time, gamma);
  // // solver.solve<explicit_marching_type>(explicit_march, target_residual, frame_time, gamma);

};
