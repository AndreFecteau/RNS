#define HYPERBOLIC
#define VISCOUS
#define SOURCE
#define IMPLICIT
// #define EXPLICIT
// #define MANUFACTURED
// #define

#include "Solver/Solver.h"
#include "Low_Mach_Solver/RK4_Low_Mach_Solver.h"
#include "Usefull_Headers/Variable_Vector_Isolator.h"
#include "Gnuplot_RNS/Gnuplot_Primitive_Variables.h"
#include "Solver/Implicit_Marching.h"
#include "Solver/Explicit_Marching.h"
#include <iomanip>
// #include "Eigen/Core"
// #include "Eigen/Dense"
// #include <cmath>


typedef Eigen::Matrix<long double, 4, 1> Vector_type;
using matrix_type = Eigen::Matrix<long double, 4,4>;
using global_solution_vector_type = std::vector<Vector_type>;
using solution_vector_type = typename global_solution_vector_type::value_type;
#include "Usefull_Headers/Initial_Conditions.h"

// using implicit_marching_type = Implicit_Marching<global_solution_vector_type, matrix_type>;
// using explicit_marching_type = Explicit_Marching<global_solution_vector_type, matrix_type>;

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
  std::cout << std::setprecision(10);
  double Pr = 0.75;
  double Le = 0.3;
  double mf = 0.005;
  double gamma = 1.4;
  double Q_low_mach = 9.0;
  double theta_low_mach =500.0/9.0;
  double Q = Q_low_mach/(mf*mf*(gamma-1));
  double theta =theta_low_mach/(gamma*mf*mf);


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
  double CFL =  1e3 ;
  double per_FL = 64.0;
  double frame_time = 1e0;
  double dx = 1.0/per_FL;
  double domaine_length = 250;
  int    number_of_cells;
  global_solution_vector_type initial_solution;

  #if defined(MANUFACTURED)
  std::ofstream gnu_input_file;
  gnu_input_file.open("Convergence_Plot.dat", std::ios_base::app);
  gnu_input_file << "#number_of_cells residual time" << std::endl;
  #endif
  // initial_solution.resize(number_of_cells);
  // manufactured_solution(number_of_cells, initial_solution, x_max, x_min, dx);


  lambda_max = 95600;
  lambda_min = 94500;
  lambda_run = 945;


  std::cout << "//////////////////////////" << std::endl;
  std::cout << "Setting Initial Conditions " << std::endl;
  std::cout << "//////////////////////////" << std::endl;
  // straight_line(number_of_cells, initial_solution, x_max, x_min, mf, gamma);
  // case_4(frame_time, number_of_cells, initial_solution, gamma, x_max, x_min);
  // RK4_low_mach_initial_conditions(lambda, number_of_cells, initial_solution, Le, Q_low_mach,
  //                                 theta_low_mach, T_ignition, gamma, x_max, mf, dx, domaine_length);
#if defined(MANUFACTURED)
  manufactured_solution(number_of_cells, initial_solution, x_max, x_min, dx);
#endif

  std::string filename = "Movie/Plot_HLLE_" + tostring(per_FL) + "_"
                                              + tostring(domaine_length) + "_"
                                              + tostring(log10(CFL)) + "_";

  auto solver = Solver<global_solution_vector_type, matrix_type>(initial_solution, filename);

#if defined(EXPLICIT)
  using marching_type = Explicit_Marching<global_solution_vector_type, matrix_type>;
  auto march = marching_type(Pr, Le, Q, theta, mf, gamma,
                            number_of_cells, CFL,
                            (x_max - x_min)/number_of_cells);
#endif
#if defined(IMPLICIT)
  using marching_type = Implicit_Marching<global_solution_vector_type, matrix_type>;
  auto march = marching_type(Pr, Le, Q, theta, mf, gamma,
                            number_of_cells, CFL,
                            (x_max - x_min)/number_of_cells, Theta, zeta);
#endif

  // while(1 < 2) {


  plot<global_solution_vector_type>(filename+"0",
                                    initial_solution, (x_max - x_min)/number_of_cells);


  bool check = solver.solve<marching_type>(march, target_residual, frame_time, gamma, lambda_run);
  // bisection_lambda(lambda_min, lambda_max, lambda_run, check);
  // }

};
