#define HYPERBOLIC
#define VISCOUS
#define SOURCE
// #define MANUFACTURED
#define RECENTER_FLAME
// #define IMPLICIT
// #define EXPLICIT

#include <iomanip>
#include <fenv.h>
#include "Solver/Solver.h"
#include "Low_Mach_Solver/RK4_Low_Mach_Solver.h"
#include "Usefull_Headers/Variable_Vector_Isolator.h"
#include "Gnuplot_RNS/Gnuplot_Primitive_Variables.h"
#include "Solver/Implicit_Marching.h"
// #include "Solver/Implicit_Marching_4th_Order.h"
#include "Solver/Explicit_Marching.h"
#include "Usefull_Headers/Handle_Itterative_Lambda.h"
#include "Physical_Property/Non_Dimensional_Navier_Stokes.h"
#include "Grid/Grid1D.h"
#include "Implicit_Flux_and_Sources/Variable_Implicit_Scheme.h"
#include "Usefull_Headers/Initial_Conditions.h"

int main(){
  std::cout << std::setprecision(10);

  using scalar_type = double;
  using size_type = size_t;
  typedef Eigen::Matrix<scalar_type, 4, 1> Vector_type;
  using matrix_type = Eigen::Matrix<scalar_type, 4,4>;
  using global_solution_vector_type = std::vector<Vector_type>;
//
  using flow_properties_type = Non_Dimensional_Navier_Stokes<scalar_type>;
  using grid_type = Grid1D<scalar_type, size_type, global_solution_vector_type, matrix_type>;
  using flux_type = Variable_Implicit_Scheme<grid_type>;
  using time_stepping_type = Implicit_Marching<grid_type, flow_properties_type>;
  using solver_type = Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>;

  std::cout << "//////////////////////////" << std::endl;
  std::cout << "Setting Initial Conditions " << std::endl;
  std::cout << "//////////////////////////" << std::endl;
  std::string filename;
  solver_type solver;
{
  flow_properties_type flow;
{
  scalar_type Pr                = 0.75;
  scalar_type Le                = 0.3;
  scalar_type Q_low_mach        = 9.0;
  scalar_type beta              = 5;
  scalar_type lambda            = 0.0;
  scalar_type mf                = 0.005;
  scalar_type T_ignition_scalar = 1.02;
  scalar_type gamma             = 1.4;
  scalar_type theta_low_mach    = beta*(1+Q_low_mach)*(1+Q_low_mach)/Q_low_mach;
  flow = flow_properties_type(Pr, Le, Q_low_mach, theta_low_mach, mf, gamma, lambda, T_ignition_scalar);
}

  grid_type grid;
{
  scalar_type x_min = 0.0;
  scalar_type domaine_length = 500;
  scalar_type x_max = x_min + domaine_length;
  scalar_type per_FL = 256.0;
  size_type number_of_cells = domaine_length * per_FL;
  global_solution_vector_type initial_solution = global_solution_vector_type(number_of_cells,
                                                                             Vector_type::Zero());
  grid = grid_type(x_min, x_max, initial_solution);
}

{
  scalar_type Theta = 1.0;
  scalar_type zeta = 0.0;
  scalar_type target_residual = 1e-15;
  scalar_type CFL = 5e8;
  scalar_type frame_time = 5e3;
  RK4_CJ_point(flow, grid);
  load_from_file(flow, grid, "Movie/Plot12_256_500_1260.dat");
  filename = "Movie/Plot_case_1_";
  std::cout << filename << std::endl;
  plot<grid_type>(filename+"0", grid.global_solution_vector, (grid.x_max - grid.x_min)/grid.number_of_cells());
  scalar_type flame_location = 250;
  solver= solver_type(flow, grid, frame_time, target_residual, CFL, Theta, zeta, filename, flame_location);
}
}
  bool old_check1 = 0;
  bool old_check2 = 1;
  bool old_check3 = 0;
  solver.print_stats();
  scalar_type lambda_run = solver.get_lambda();
  // solver.set_lambda(124889);
  scalar_type lambda_max = solver.get_lambda()*1.1;
  scalar_type lambda_min = solver.get_lambda()*0.9;
  int number_of_frames = 400;
  solver.solve(number_of_frames);

  while(fabs(lambda_min - lambda_max) > 1e2) {
    bool check;
    number_of_frames = 100;
    check = solver.solve(number_of_frames);
    scalar_type lambda_run = solver.get_lambda();
    bisection_lambda(lambda_min, lambda_max, lambda_run, check);
    add_lambda_gap(check, old_check1, old_check2, old_check3, lambda_min, lambda_max);
    solver.set_lambda(lambda_run);
  }

  while(fabs(lambda_min - lambda_max) > 1e-2) {
    bool check;
    number_of_frames = 40;
    check = solver.solve(number_of_frames);
    scalar_type lambda_run = solver.get_lambda();
    bisection_lambda(lambda_min, lambda_max, lambda_run, check);
    add_lambda_gap(check, old_check1, old_check2, old_check3, lambda_min, lambda_max);
    solver.set_lambda(lambda_run);
  }

  while(fabs(lambda_min - lambda_max) > 1e-8) {
    bool check;
    number_of_frames = 10;
    check = solver.solve(number_of_frames);
    scalar_type lambda_run = solver.get_lambda();
    bisection_lambda(lambda_min, lambda_max, lambda_run, check);
    add_lambda_gap(check, old_check1, old_check2, old_check3, lambda_min, lambda_max);
    solver.set_lambda(lambda_run);
  }
};
