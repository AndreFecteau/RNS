#define HYPERBOLIC
#define VISCOUS
#define SOURCE
#define IMPLICIT
// #define EXPLICIT
// #define MANUFACTURED

#include <iomanip>
#include <fenv.h>
#include "Solver/Solver.h"
#include "Low_Mach_Solver/RK4_Low_Mach_Solver.h"
#include "Usefull_Headers/Variable_Vector_Isolator.h"
#include "Gnuplot_RNS/Gnuplot_Primitive_Variables.h"
#include "Solver/Implicit_Marching.h"
// #include "Solver/Implicit_Marching_4th_Order.h"
#include "Solver/Explicit_Marching.h"

typedef Eigen::Matrix<double, 4, 1> Vector_type;
using matrix_type = Eigen::Matrix<double, 4,4>;
using global_solution_vector_type = std::vector<Vector_type>;
using solution_vector_type = typename global_solution_vector_type::value_type;
#include "Usefull_Headers/Initial_Conditions.h"

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
  // feenableexcept(FE_INVALID | FE_OVERFLOW);
  std::cout << std::setprecision(10);
  double Pr = 0.75;
  double Le = 0.3;
  double Q_low_mach = 9.0;
  double theta_low_mach =500.0/9.0;
  // double Pr = 0.75;
  // double Le = 1.0;
  // double Q_low_mach = 4;
  // double theta_low_mach = 30;
  // double Pr = 0.75;
  // double Le = 1.0;
  // double Q_low_mach = 6;
  // double theta_low_mach = 30;
  double mf = 0.01;
  double gamma = 1.4;
  double Q = Q_low_mach/(mf*mf*(gamma-1));
  double theta =theta_low_mach/(gamma*mf*mf);


  double lambda = 0.0;
  double x_min = 0.0;
  double x_max;
  double T_ignition = 1.0;
  double lambda_max;
  double lambda_min;
  double lambda_run;
  double target_residual = 1e-15;

  double Theta = 1.0;
  double zeta = 0.0;
  double CFL = 5e7;
  double per_FL = 256.0;
  double frame_time = 1e2;
  double dx = 1.0/per_FL;
  double domaine_length = 500;
  int    number_of_cells;
  global_solution_vector_type initial_solution;

  #if defined(MANUFACTURED)
  std::ofstream gnu_input_file;
  gnu_input_file.open("Convergence_Plot.dat", std::ios_base::app);
  gnu_input_file << "#number_of_cells residual time" << std::endl;
  #endif
  // initial_solution.resize(number_of_cells);
  // manufactured_solution(number_of_cells, initial_solution, x_max, x_min, dx);





  std::cout << "//////////////////////////" << std::endl;
  std::cout << "Setting Initial Conditions " << std::endl;
  std::cout << "//////////////////////////" << std::endl;
#if defined(MANUFACTURED)
  manufactured_solution(number_of_cells, initial_solution, x_max, x_min, dx);
#else
RK4_CJ_point(lambda, number_of_cells, initial_solution, Le, Q_low_mach,
             theta_low_mach, T_ignition, gamma, x_max, mf, dx, domaine_length);
// RK4_low_mach_initial_conditions(lambda, number_of_cells, initial_solution, Le, Q_low_mach,
//                                 theta_low_mach, T_ignition, gamma, x_max, mf, dx, domaine_length);
#endif
lambda_max = lambda*1.0001;
lambda_min = lambda*0.9999;
lambda_run = lambda;
  std::string filename = "Movie/Plot11_" + tostring(per_FL) + "_"
                                        + tostring(domaine_length) + "_";
auto solver = Solver<global_solution_vector_type, matrix_type>(initial_solution, filename);

while(mf < 1.0) {
  double Q = Q_low_mach/(mf*mf*(gamma-1));
  double theta =theta_low_mach/(gamma*mf*mf);
  // std::cout << "Q: " << Q << "theta: " << theta << std::endl;

for(size_t i = 0; i < 5; ++i){
#if defined(EXPLICIT)
  using marching_type = Explicit_Marching<global_solution_vector_type, matrix_type>;
  auto march = marching_type(Pr, Le, Q, theta, mf, gamma,
                            number_of_cells, CFL,
                            dx);
#endif
#if defined(IMPLICIT)
  using marching_type = Implicit_Marching<global_solution_vector_type, matrix_type>;
  auto march = marching_type(Pr, Le, Q, theta, mf, gamma,
                            number_of_cells, CFL,
                            dx, Theta, zeta);
#endif

      plot<global_solution_vector_type>(filename+"0",
      initial_solution, (x_max - x_min)/number_of_cells);


      solver.solve<marching_type>(march, target_residual, frame_time, gamma, lambda_run, 1);
}
      solver.set_bound_solution_vector(lambda_run, theta, Q, dx, mf);
      // solver.solve<marching_type>(march, target_residual, frame_time, gamma, lambda_run, 1);
// lambda_max = lambda_run*1.03;
// lambda_min = lambda_run*0.97;
//   while(fabs(lambda_min - lambda_max) > 1e2) {
// #if defined(EXPLICIT)
//   using marching_type = Explicit_Marching<global_solution_vector_type, matrix_type>;
//   march = marching_type(Pr, Le, Q, theta, mf, gamma,
//                             number_of_cells, CFL,
//                             dx);
// #endif
// #if defined(IMPLICIT)
//   using marching_type = Implicit_Marching<global_solution_vector_type, matrix_type>;
//   march = marching_type(Pr, Le, Q, theta, mf, gamma,
//                             number_of_cells, CFL,
//                             dx, Theta, zeta);
// #endif
//     bool check = solver.solve<marching_type>(march, target_residual, frame_time, gamma, lambda_run, 1);
//     bisection_lambda(lambda_min, lambda_max, lambda_run, check);
//     // CFL = 5e6;
//   }
// solver.set_bound_solution_vector(lambda_run, theta, Q, dx, mf);
// // CFL = 2e6;
       double mf_old = mf;
       mf +=0.01;
      //  CFL /=2.0;
       Q = Q_low_mach/(mf*mf*(gamma-1));
       theta =theta_low_mach/(gamma*mf*mf);
       solver.set_new_mf_to_solution_vector(lambda_run, mf_old, mf, Q);
//     std::cout << "//////////////////////////" << std::endl;
//     std::cout << mf << " : " << lambda_run << std::endl;
//     std::cout << "//////////////////////////" << std::endl;
  }

};
