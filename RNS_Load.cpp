#define HYPERBOLIC
#define VISCOUS
#define SOURCE
#define RECENTER_FLAME
// #define RIGHT_CST_EXTR
// #define LEFT_CST_EXTR
// #define MANUFACTURED

#include <iomanip>
#include <fenv.h>
#include "Solver/Solver.h"
#include "Low_Mach_Solver/RK4_Low_Mach_Solver.h"
#include "Usefull_Headers/Variable_Vector_Isolator.h"
#include "Gnuplot_RNS/Gnuplot_Primitive_Variables.h"
#include "Gnuplot_RNS/Gnuplot_Primitive_Variables_Reduced.h"
#include "Solver/Implicit_Marching.h"
#include "Usefull_Headers/Handle_Itterative_Lambda.h"
#include "Physical_Property/Non_Dimensional_Navier_Stokes.h"
#include "Grid/Grid1D.h"
#include "Implicit_Flux_and_Sources/Implicit_Centered_Difference_2nd_Order.h"
#include "Implicit_Flux_and_Sources/Implicit_HLLE.h"
#include "Usefull_Headers/Initial_Conditions.h"
#include "Usefull_Headers/Read_from_file.h"

int main(){
  std::cout << std::setprecision(20);

  using scalar_type = double;
  using size_type = size_t;
  typedef Eigen::Matrix<scalar_type, 4, 1> Vector_type;
  using matrix_type = Eigen::Matrix<scalar_type, 4,4>;
  using global_solution_vector_type = std::vector<Vector_type>;
//
  using flow_properties_type = Non_Dimensional_Navier_Stokes<scalar_type>;
  using grid_type = Grid1D<scalar_type, size_type, global_solution_vector_type, matrix_type>;
  using flux_type = Implicit_Centered_Difference_2nd_Order<grid_type, flow_properties_type>;
  // using flux_type = Implicit_HLLE<grid_type, flow_properties_type>;
  using time_stepping_type = Implicit_Marching<grid_type, flow_properties_type>;
  using solver_type = Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>;

  std::cout << "//////////////////////////" << std::endl;
  std::cout << "Setting Initial Conditions " << std::endl;
  std::cout << "//////////////////////////" << std::endl;
  solver_type solver;
  // unserialize_to_file(solver, "dat_saves/CJ_point/Pr075_Le03_Q9_b5_g14/D500_F250_R16/Solution");
  // unserialize_to_file(solver, "Movie/Low_Mach_R16_23");
  unserialize_to_file(solver, "dat_saves/CJ_point/Pr075_Le03_Q9_b5_g14/D500_F250_R16/Solution");
  std::cout << "Restarting Simulation With:" << std::endl;
  solver.print_stats();
  // solver.plot_limiter();

  // solver.add_space_in_back(10000);
  // solver.frame_time = 5e5;
  // solver.CFL = 5e8;

  solver.change_filename("Movie/Delete_");
  solver.change_frame_time(1e0);
  solver.change_CFL(1e4);
  // solver.reset_frame_number();
  // solver.refine(64);
  // plot_reduced<grid_type>("Movie/Delete_0", solver.grid, solver.flow, 2500*16);
  // scalar_type lambda_run = solver.get_lambda();
  // solver.set_lambda(lambda_run*1.1);
  // solver.frame_time = 1e-1;
  // solver.CFL = 1;
  int number_of_frames = 10;
  solver.solve(number_of_frames);


  // bool old_check1 = 0;
  // bool old_check2 = 1;
  // bool old_check3 = 0;
// {
//
//   scalar_type lambda_max = solver.get_lambda()*1.005;
//   scalar_type lambda_min = solver.get_lambda()*0.995;
//
//   while(fabs(lambda_min - lambda_max) > 1e2) {
//     bool check;
//     int number_of_frames = 3;
//     check = solver.solve(number_of_frames);
//     scalar_type lambda_run = solver.get_lambda();
//     bisection_lambda(lambda_min, lambda_max, lambda_run, check);
//     add_lambda_gap(check, old_check1, old_check2, old_check3, lambda_min, lambda_max);
//     solver.set_lambda(lambda_run);
// }
//   while(fabs(lambda_min - lambda_max) > 1e-2) {
//     bool check;
//     int number_of_frames = 3;
//     check = solver.solve(number_of_frames);
//     scalar_type lambda_run = solver.get_lambda();
//     bisection_lambda(lambda_min, lambda_max, lambda_run, check);
//     add_lambda_gap(check, old_check1, old_check2, old_check3, lambda_min, lambda_max);
//     solver.set_lambda(lambda_run);
// }
//   while(fabs(lambda_min - lambda_max) > 1e-8) {
//     bool check;
//     int number_of_frames = 3;
//     check = solver.solve(number_of_frames);
//     scalar_type lambda_run = solver.get_lambda();
//     bisection_lambda(lambda_min, lambda_max, lambda_run, check);
//     add_lambda_gap(check, old_check1, old_check2, old_check3, lambda_min, lambda_max);
//     solver.set_lambda(lambda_run);
// }
// }
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// {
//   // solver.add_space_in_back(500);
//   // solver.rename_file("Movie/Case_1_D1000_L250_R256_");
//   solver.reset_frame_number();
//   scalar_type lambda_max = solver.get_lambda()*1.005;
//   scalar_type lambda_min = solver.get_lambda()*0.995;
// //
//   while(fabs(lambda_min - lambda_max) > 1e2) {
//     bool check;
//     int number_of_frames = 20;
//     check = solver.solve(number_of_frames);
//     scalar_type lambda_run = solver.get_lambda();
//     bisection_lambda(lambda_min, lambda_max, lambda_run, check);
//     solver.set_lambda(lambda_run);
// }
//   while(fabs(lambda_min - lambda_max) > 1e-2) {
//     bool check;
//     int number_of_frames = 20;
//     check = solver.solve(number_of_frames);
//     scalar_type lambda_run = solver.get_lambda();
//     bisection_lambda(lambda_min, lambda_max, lambda_run, check);
//     solver.set_lambda(lambda_run);
// }
//   while(fabs(lambda_min - lambda_max) > 1e-8) {
//     bool check;
//     int number_of_frames = 10;
//     check = solver.solve(number_of_frames);
//     scalar_type lambda_run = solver.get_lambda();
//     bisection_lambda(lambda_min, lambda_max, lambda_run, check);
//     solver.set_lambda(lambda_run);
// }
// }
// ///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////////////
// {
//   solver.add_space_in_back(1000);
//   solver.rename_file("Movie/Case_1_D2000_L250_R256_");
//   solver.reset_frame_number();
//   scalar_type lambda_max = solver.get_lambda()*1.005;
//   scalar_type lambda_min = solver.get_lambda()*0.995;
// //
//   while(fabs(lambda_min - lambda_max) > 1e2) {
//     bool check;
//     int number_of_frames = 20;
//     check = solver.solve(number_of_frames);
//     scalar_type lambda_run = solver.get_lambda();
//     bisection_lambda(lambda_min, lambda_max, lambda_run, check);
//     solver.set_lambda(lambda_run);
// }
//   while(fabs(lambda_min - lambda_max) > 1e-2) {
//     bool check;
//     int number_of_frames = 20;
//     check = solver.solve(number_of_frames);
//     scalar_type lambda_run = solver.get_lambda();
//     bisection_lambda(lambda_min, lambda_max, lambda_run, check);
//     solver.set_lambda(lambda_run);
// }
//   while(fabs(lambda_min - lambda_max) > 1e-8) {
//     bool check;
//     int number_of_frames = 10;
//     check = solver.solve(number_of_frames);
//     scalar_type lambda_run = solver.get_lambda();
//     bisection_lambda(lambda_min, lambda_max, lambda_run, check);
//     solver.set_lambda(lambda_run);
// }
// }
// ///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////////////
// {
//   solver.add_space_in_back(1000);
//   solver.rename_file("Movie/Case_1_D3000_L250_R256_");
//   solver.reset_frame_number();
//   scalar_type lambda_max = solver.get_lambda()*1.005;
//   scalar_type lambda_min = solver.get_lambda()*0.995;
// //
//   while(fabs(lambda_min - lambda_max) > 1e2) {
//     bool check;
//     int number_of_frames = 20;
//     check = solver.solve(number_of_frames);
//     scalar_type lambda_run = solver.get_lambda();
//     bisection_lambda(lambda_min, lambda_max, lambda_run, check);
//     solver.set_lambda(lambda_run);
// }
//   while(fabs(lambda_min - lambda_max) > 1e-2) {
//     bool check;
//     int number_of_frames = 20;
//     check = solver.solve(number_of_frames);
//     scalar_type lambda_run = solver.get_lambda();
//     bisection_lambda(lambda_min, lambda_max, lambda_run, check);
//     solver.set_lambda(lambda_run);
// }
//   while(fabs(lambda_min - lambda_max) > 1e-8) {
//     bool check;
//     int number_of_frames = 10;
//     check = solver.solve(number_of_frames);
//     scalar_type lambda_run = solver.get_lambda();
//     bisection_lambda(lambda_min, lambda_max, lambda_run, check);
//     solver.set_lambda(lambda_run);
// }
// }
};
