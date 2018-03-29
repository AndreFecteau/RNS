#ifndef SOLVER_H
#define SOLVER_H

#include <math.h>
#include <chrono>
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>

#include "../Serialization/Serialization_Eigen.h"
#include "../Serialization/Serialize.h"
#include "../Gnuplot_RNS/Gnuplot_Primitive_Variables.h"


template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
class Solver {

using global_solution_vector_type = typename grid_type::global_solution_vector_type;
using solution_vector_type = typename grid_type::global_solution_vector_type::value_type;
using scalar_type = typename grid_type::scalar_type;

 public:
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  Solver() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  Solver(const Solver&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  Solver(Solver&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  Solver& operator=(const Solver&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  Solver& operator=(Solver&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Constructor setting up required inputs.
  Solver(flow_properties_type flow_properties_in, grid_type grid_in, scalar_type frame_time_in,
         scalar_type target_residual_in, scalar_type CFL_in, scalar_type Theta_in,
         scalar_type zeta_in) : flow_properties(flow_properties_in), grid(grid_in),
         frame_time(frame_time_in), target_residual(target_residual_in), CFL(CFL_in),
         time_stepping(time_stepping_type(Theta_in, zeta_in)) {}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Runs the simulation. Outputs frames in a folder "../Movie".
  /// \param marching The type of marching method passed to the solver.
  ///        ex:Implicit or Explicit.
  /// \param target_residual Residual where the solver will stop.
  bool solve(int number_of_frames);

  // void set_new_mf_to_solution_vector(double &lambda, double mf_old, double mf, double Q){
  //   Variable_Vector_Isolator<solution_vector_type> temp_0 = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[0], 1.4);
  //   // lambda *= mf_old*mf_old/(mf*mf);
  //   std::cout << "lambda: " << lambda << " mf: " << mf << std::endl;
  //   for (int i = 0; i < global_solution_vector.size(); ++i){
  //     Variable_Vector_Isolator<solution_vector_type> temp = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i], 1.4);
  //     double rho  = temp.rho();
  //     double u    = temp.u();
  //     double p    = temp.p()/temp_0.p()*(1.0/(mf*mf*1.4));
  //     global_solution_vector[i] << rho, rho * u, p / 0.4 + rho * u * u / 2.0, rho*temp.Y();
  //   }
  // double gamma =1.4;
  // double p_0 = 1/(gamma*mf*mf);
  // double rho_0 = 1.0;
  // double u_0 = 1.0;
  // double zeta =(2*pow(mf,2)*Q + pow(-1 + pow(mf,2),2)*p_0/rho_0*gamma - 2*pow(mf,2)*Q*pow(gamma,2))/(pow(mf,4)*p_0/rho_0*gamma*pow(1 + gamma,2));
  // double p_inf = p_0*((1 + pow(mf,2)*gamma)/(1 + gamma) + pow(mf,2)*gamma*sqrt(zeta));
  // double rho_inf = rho_0/((1 + pow(mf,2)*gamma)/(pow(mf,2)*(1 + gamma)) - sqrt(zeta));
  // double u_inf = u_0*((1 + pow(mf,2)*gamma)/(pow(mf,2)*(1 + gamma)) - sqrt(zeta));
  //   std::cout << "zeta:" << zeta << "p:" << p_inf << " rho: " << rho_inf << " u: " << u_inf << std::endl;
  //
  //     global_solution_vector[0] << 1.0,
  //     1.0 * 1.0,
  //     (1.0 / (gamma * mf * mf)) /
  //     (gamma - 1.0) + 1.0 * 1.0 * 1.0 * 0.5,
  //     1.0 * 1.0;
  //     global_solution_vector[global_solution_vector.size()-1] << rho_inf,
  //     rho_inf * u_inf,
  //     p_inf / (gamma - 1.0) + rho_inf * u_inf * u_inf * 0.5,
  //     rho_inf * 0;
  // }

  // void set_bound_solution_vector(double &lambda, double &theta, double &Q, double &dx, double &mf){
  //   for(int i = 0; i < 500; ++i) {
  //     global_solution_vector[i] = (global_solution_vector[500]);
  //   }
  //   Variable_Vector_Isolator<solution_vector_type> temp_0 = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[0], 1.4);
  //   // lambda /= temp_0.rho()*temp_0.u()*temp_0.u();
  //   // theta /=  (temp_0.u()*temp_0.u());
  //   // Q /=  (temp_0.u()*temp_0.u());
  //   // dx *= temp_0.rho()*temp_0.u();
  //   // std::cout << "theta: " << theta << " Q: " << Q << " mf: " << mf << std::endl;
  //   // std::cout << "rho*u: " << temp_0.rho()*temp_0.u() << std::endl;
  //   // for (int i = 0; i < global_solution_vector.size(); ++i){
  //   //   Variable_Vector_Isolator<solution_vector_type> temp = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i], 1.4);
  //   //   double rho  = temp.rho() / temp_0.rho();
  //   //   double u    = temp.u() / temp_0.u();
  //   //   double p    = temp.p();// / (temp_0.rho()*temp_0.u()*temp_0.u());
  //   //   global_solution_vector[i] << rho, rho * u, p / 0.4 + rho * u * u / 2.0, rho*temp.Y();
  //   // }
  //   // temp_0 = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[0], 1.4);
  //   // mf = temp_0.u()/sqrt(1.4*temp_0.p()/temp_0.rho());
  //   // std::cout << "theta: " << theta << " Q: " << Q << " mf: " << mf << std::endl;
  // }

  void recenter_solution_plus(int position){
    auto global_solution_vector_temp = grid.global_solution_vector;
    position = 0.5*grid.number_of_cells-position;
    for(size_t i = 1; i < position; ++i){
      global_solution_vector_temp[i] = grid.global_solution_vector[0];
    }
    for(size_t i = position; i < grid.number_of_cells-1; ++i){
      global_solution_vector_temp[i] = grid.global_solution_vector[i-position];
    }
    grid.global_solution_vector = global_solution_vector_temp;
  }

  void recenter_solution_minus(int position){
    auto global_solution_vector_temp = grid.global_solution_vector;
    position -= 0.5*grid.number_of_cells;
    for(size_t i = 1; i < grid.number_of_cells - position; ++i){
      global_solution_vector_temp[i] = grid.global_solution_vector[i+position];
    }
    for(size_t i = grid.number_of_cells - position; i < grid.number_of_cells-1; ++i){
      global_solution_vector_temp[i] = grid.global_solution_vector[grid.number_of_cells-1];
    }
    grid.global_solution_vector = global_solution_vector_temp;
  }

  // void reset_close_to_bound(){
  //   auto global_solution_vector_temp = global_solution_vector;
  //   for(size_t i = 1; i < 0.4 * global_solution_vector.size(); ++i){
  //     global_solution_vector_temp[i] = global_solution_vector[0];
  //   }
  //   for(size_t i = global_solution_vector.size()*0.8; i < global_solution_vector.size()-1; ++i){
  //     global_solution_vector_temp[i] = global_solution_vector[global_solution_vector.size()-1];
  //   }
  //   global_solution_vector = global_solution_vector_temp;
  // }
  template<typename Archive>
  void serialize(Archive& archive) {
    archive(flow_properties, grid, frame_time, target_residual, CFL, time_stepping, filename);
  }
 private:
  flow_properties_type flow_properties;
  grid_type grid;
  scalar_type frame_time;
  scalar_type target_residual;
  scalar_type CFL;
  time_stepping_type time_stepping;
  global_solution_vector_type global_solution_vector_backup;
  std::string filename;
  double current_time = 0.0;
  int current_frame = 0;

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  int flame_position_algorithm(double gamma);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  // void manufactured_solution_residual(time_stepping_type march);

};

///////////////////////////////////////////////////////////////////////////////
// Solver
///////////////////////////////////////////////////////////////////////////////
template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
bool Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
solve(int number_of_frames) {

  std::cout << "//////////////////////" << std::endl;
  std::cout << "Lambda = " << flow_properties.lambda << std::endl;
  std::cout << "//////////////////////" << std::endl;

  int old_position;
  double residual = std::numeric_limits<double>::max();
  (void)residual;
  (void)target_residual;
  int i = 0;
    plot<grid_type>(filename + std::to_string(static_cast<int>(2000)), grid.global_solution_vector, grid.dx());
    double frame_time_temp = frame_time;
    double frame_CFL = CFL;

  // while (residual > target_residual){
    while (i < number_of_frames){
    old_position = flame_position_algorithm(flow_properties.gamma);
    auto start = std::chrono::high_resolution_clock::now();
    global_solution_vector_backup = grid.global_solution_vector;
    residual = time_stepping.template timemarch<flux_type>(flow_properties, grid, frame_CFL, frame_time_temp);
    int position = 0;
    auto var_vec = Variable_Vector_Isolator<solution_vector_type>(grid.global_solution_vector[position], flow_properties.gamma);
    while (var_vec.rho() > 0.5) {
    ++position;
    var_vec = Variable_Vector_Isolator<solution_vector_type>(grid.global_solution_vector[position], flow_properties.gamma);
    }
    if(isnan(residual) || residual > 1e10){
      grid.global_solution_vector = global_solution_vector_backup;
      frame_CFL *= 0.5;
      frame_time_temp *= 0.5;
    }else if(position < grid.number_of_cells*0.1 || position > grid.number_of_cells*0.9) {
      grid.global_solution_vector = global_solution_vector_backup;
      frame_time_temp *= 0.5;
    } else {
      current_time += frame_time_temp;
      std::cout << "Time = " <<  current_time << " : ";
      plot<global_solution_vector_type>(filename + std::to_string(static_cast<int>(current_frame)+1), grid.global_solution_vector,grid.dx());
      serialize_to_file(*this, filename + std::to_string(static_cast<int>(current_frame)+1));
      current_frame++;
      auto finish = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = finish - start;
      std::cout << "Elapsed time: " << elapsed.count() << "\n";
      ++i;
      frame_time_temp = frame_time;
      frame_CFL = CFL;
    }
  }

#if defined(MANUFACTURED)
  manufactured_solution_residual<marching_type>(march);
#endif
  int position = 0;
  auto var_vec = Variable_Vector_Isolator<solution_vector_type>(grid.global_solution_vector[position], 1.4);
  while (var_vec.rho() > 0.5) {
  ++position;
  var_vec = Variable_Vector_Isolator<solution_vector_type>(grid.global_solution_vector[position], 1.4);
  }
  std::cout << "position: " << position << std::endl;
  if(position < 0.50*grid.number_of_cells) {
  std::cout << "moved_plus: " <<  std::endl;
    recenter_solution_plus(position);
  }
  if(position > 0.50*grid.number_of_cells) {
  std::cout << "moved_minus: " << std::endl;
    recenter_solution_minus(position);
  }
  if(position < old_position) {
    return 0;
  } else {
    return 1;
  }
}

///////////////////////////////////////////////////////////////////////////////
// Manufactured Solution Residual
///////////////////////////////////////////////////////////////////////////////
// template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
// void Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
// manufactured_solution_residual(time_stepping_type march) {
//   std::cout << march.get_dx() << std::endl;
//   for (size_t i = 0; i < initial_solution.size(); ++i) {
//    double x = (i+0.5)*march.get_dx();
//     initial_solution[i] << (-0.45*tanh(4.0 * x - 10.0) + 0.55),
//                           (-0.45*tanh(4.0 * x - 10.0) + 0.55)*(4.5*tanh(4.0 * x - 10.0) + 5.5),
//                           2.0*tanh(4.0*x - 10.0) + 70000,
//                           (-0.45*tanh(4.0 * x - 10.0) + 0.55)*(-0.5*tanh(x - 8.0/4.0) + 0.5);
//   }
//   double convergence = 0.0;
//   for(size_t i = 0; i < initial_solution.size(); ++i) {
//    for(int j = 0; j < 4; ++j){
//      convergence += std::pow(std::fabs(initial_solution[i][j]-global_solution_vector[i][j]),2)*march.get_dx();
//    }
//   }
//
//   convergence = std::sqrt(convergence);
//   std::cout << "convergence:" << convergence << std::endl;
//   std::ofstream gnu_input_file;
//   gnu_input_file.open("Convergence_Plot.dat", std::ios_base::app);
//   gnu_input_file << global_solution_vector.size() << " " <<  convergence << " " << current_time << std::endl;
// }


///////////////////////////////////////////////////////////////////////////////
// Flame Position Algorithm
///////////////////////////////////////////////////////////////////////////////
template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
int Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
flame_position_algorithm(double gamma) {
  int i = 0;
  auto var_vec = Variable_Vector_Isolator<solution_vector_type>(grid.global_solution_vector[0], gamma);
  while (var_vec.rho() > 0.5) {
  ++i;
  var_vec = Variable_Vector_Isolator<solution_vector_type>(grid.global_solution_vector[i], gamma);
  }
  return i;
}

#endif //#ifndef SOLVER_H
