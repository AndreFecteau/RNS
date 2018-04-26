#ifndef SOLVER_H
#define SOLVER_H

#include <math.h>
#include <chrono>
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>

#include "../Serialization/Serialization_Eigen.h"
#include "../Serialization/Serialize.h"
#include "../Gnuplot_RNS/Gnuplot_Primitive_Variables.h"


template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
class Solver {

using global_solution_vector_type = typename grid_type::global_solution_vector_type;
using solution_vector_type = typename grid_type::global_solution_vector_type::value_type;
using scalar_type = typename grid_type::scalar_type;
using size_type = typename grid_type::size_type;

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
         scalar_type zeta_in, std::string filename_in, scalar_type flame_location_in) :
         frame_time(frame_time_in), CFL(CFL_in), flow(flow_properties_in),
         grid(grid_in), target_residual(target_residual_in),
         time_stepping(time_stepping_type(Theta_in, zeta_in)),
         filename(filename_in), flame_location(flame_location_in) {}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Runs the simulation. Outputs frames in a folder "../Movie".
  /// \param marching The type of marching method passed to the solver.
  ///        ex:Implicit or Explicit.
  /// \param target_residual Residual where the solver will stop.
  bool solve(size_type number_of_frames);

  double get_lambda(){return flow.lambda;}

  void recenter_solution(scalar_type flame_location);

  void add_space_in_back(double space);

  void add_space_in_front(double space);

  void print_stats();

  void rename_file(std::string new_filename) {
    filename = new_filename;
    std::cout << "filename changed to: " << filename << std::endl;
  }

  void reset_frame_number() {
    current_frame = 0;
    std::cout << "Reset Frame Number" << std::endl;
  }

  void set_lambda(scalar_type lambda) {
    flow.lambda = lambda;
    std::cout << "lambda changed to: " << flow.lambda << std::endl;
  }

  void refine(size_type per_flame_length);

  template<typename Archive>
  void serialize(Archive& archive) {
    archive(flow, grid, frame_time, target_residual, CFL, time_stepping, filename, current_time, current_frame, flame_location);
  }

  scalar_type frame_time;
  scalar_type CFL;
 private:
  flow_properties_type flow;
  grid_type grid;
  scalar_type target_residual;
  time_stepping_type time_stepping;
  global_solution_vector_type global_solution_vector_backup;
  scalar_type current_time = 0.0;
  size_type current_frame = 0;
  std::string filename;
  scalar_type flame_location;

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  size_type flame_position_algorithm(double gamma);

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
solve(size_type number_of_frames) {

  std::cout << "//////////////////////" << std::endl;
  std::cout << "Lambda = " << flow.lambda << std::endl;
  std::cout << "//////////////////////" << std::endl;

  size_type old_position = 0;
  double residual = std::numeric_limits<double>::max();
  (void)residual;
  (void)target_residual;
  size_type i = 0;
    plot<grid_type>(filename + std::to_string(static_cast<size_type>(2000)), grid.global_solution_vector, grid.dx());
    double frame_time_temp = frame_time;
    double frame_CFL = CFL;
    size_type position = 0;

  // while (residual > target_residual){
    while (i < number_of_frames){
    old_position = flame_position_algorithm(flow.gamma);
    auto start = std::chrono::high_resolution_clock::now();
    global_solution_vector_backup = grid.global_solution_vector;
    residual = time_stepping.template timemarch<flux_type>(flow, grid, frame_CFL, frame_time_temp);
    position = 0;
    auto var_vec = Variable_Vector_Isolator<grid_type>(grid.global_solution_vector[position], flow.gamma);
    while (var_vec.rho() > 0.5) {
    ++position;
    var_vec = Variable_Vector_Isolator<grid_type>(grid.global_solution_vector[position], flow.gamma);
    }
    if(isnan(residual) || residual > 1e10){
      grid.global_solution_vector = global_solution_vector_backup;
      frame_CFL *= 0.5;
      frame_time_temp *= 0.5;
      std::cout << "." << std::flush;
#if defined(RECENTER_FLAME)
    }else if(position < 100*grid.per_FL() || position > grid.number_of_cells()*0.9) {
      grid.global_solution_vector = global_solution_vector_backup;
      frame_time_temp *= 0.5;
      std::cout << ":" << std::flush;
#endif
    } else {
      current_time += frame_time_temp;
      std::cout << "Frame: " << current_frame << " Frame_time = " <<  frame_time_temp << " Residual: " << residual << std::endl;
      std::cout << "position: " << position - flame_location*grid.per_FL();
      current_frame++;
      plot<grid_type>(filename + std::to_string(static_cast<size_type>(current_frame)), grid.global_solution_vector,grid.dx());
      serialize_to_file(*this, filename + std::to_string(static_cast<size_type>(current_frame)));
      auto finish = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = finish - start;
      std::cout << " Elapsed time: " << elapsed.count() << "\n";
      ++i;
      frame_time_temp = frame_time;
      frame_CFL = CFL;
#if defined(RECENTER_FLAME)
      recenter_solution(flame_location);
#endif

    }
  }

int h = std::max(static_cast<int>(5*grid.per_FL()), static_cast<int>(5*grid.per_FL()-(position - flame_location*grid.per_FL())));
std::cout << h << std::endl;
if(grid.global_solution_vector[h][1] / grid.global_solution_vector[h][0] < 1.0) {
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
//    for(size_type j = 0; j < 4; ++j){
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
typename grid_type::size_type Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
flame_position_algorithm(double gamma) {
  size_type i = 0;
  auto var_vec = Variable_Vector_Isolator<grid_type>(grid.global_solution_vector[0], gamma);
  while (var_vec.rho() > 0.5) {
  ++i;
  var_vec = Variable_Vector_Isolator<grid_type>(grid.global_solution_vector[i], gamma);
  }
  return i;
}

template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
void Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
recenter_solution(scalar_type flame_location) {
  auto global_solution_vector_temp = grid.global_solution_vector;
  size_type position = 0;
  auto var_vec = Variable_Vector_Isolator<grid_type>(grid.global_solution_vector[position], 1.4);
  while (var_vec.rho() > 0.5) {
  ++position;
  var_vec = Variable_Vector_Isolator<grid_type>(grid.global_solution_vector[position], 1.4);
  }
  int delta_position = flame_location*grid.per_FL()-position;
  if(delta_position > 0){
    for(size_type i = 1; i < abs(delta_position); ++i){
      global_solution_vector_temp[i] = grid.global_solution_vector[0];
    }
    for(size_type i = abs(delta_position); i < grid.number_of_cells()-1; ++i){
      global_solution_vector_temp[i] = grid.global_solution_vector[i-delta_position];
    }
    grid.global_solution_vector = global_solution_vector_temp;
  } else {
    for(size_type i = 1; i < grid.number_of_cells() - abs(delta_position); ++i){
      global_solution_vector_temp[i] = grid.global_solution_vector[i+abs(delta_position)];
    }
    for(size_type i = grid.number_of_cells() - abs(delta_position); i < grid.number_of_cells()-1; ++i){
      global_solution_vector_temp[i] = grid.global_solution_vector[grid.number_of_cells()-1];
    }
    grid.global_solution_vector = global_solution_vector_temp;
  }
}

template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
void Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
add_space_in_back(double space) {
  grid.global_solution_vector.resize(grid.number_of_cells() + space/grid.dx());
  grid.x_max += space;
  for(size_t i = grid.number_of_cells() - space/grid.dx(); i < grid.number_of_cells(); ++i) {
    grid.global_solution_vector[i] = grid.global_solution_vector[i-1];
  }
  std::cout << "Added space in domaine, New domaine size: " << grid.domaine_length() << std::endl;
}

template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
void Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
add_space_in_front(double space) {
  global_solution_vector_type temp_global_solution_vector;
  temp_global_solution_vector.resize(grid.number_of_cells() + space/grid.dx());
  grid.x_max += space;
  for(size_t i = 0; i < space/grid.dx(); ++i){
    temp_global_solution_vector[i] = grid.global_solution_vector[0];
  }
  for(size_t i = space/grid.dx(); i < temp_global_solution_vector.size(); ++i) {
    temp_global_solution_vector[i] = grid.global_solution_vector[i - space/grid.dx()];
  }
  grid.global_solution_vector = temp_global_solution_vector;
  std::cout << "Added space in domaine, New domaine size: " << grid.domaine_length() << std::endl;
}

template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
void Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
print_stats() {
    std::cout << "Pr: " << flow.Pr << "\nLe: " << flow.Le << "\ngamma: " << flow.gamma << std::endl;
    std::cout << "mf: " << flow.mf << "\nQ_low_mach: " << flow.Q_low_mach << "\ntheta_low_mach: " << flow.theta_low_mach << std::endl;
    std::cout << "T_ignition: " << flow.T_ignition_scalar << "\nlambda: " << flow.lambda  << std::endl;
    std::cout << "x_min: " << grid.x_min << "\nx_max: " << grid.x_max << "\nnumber_of_cells: " << grid.number_of_cells() << std::endl;
}

template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
void Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
refine(size_type per_flame_length) {
  global_solution_vector_type temp_global_solution_vector;
  if(per_flame_length > grid.per_FL()){
    if(per_flame_length % grid.per_FL() != 0) {
      std::cout << "cant refine to asked amount" << std::endl;
    } else {
      temp_global_solution_vector.resize(grid.domaine_length() * per_flame_length);
      for(size_type i = 0; i < grid.number_of_cells(); ++i) {
        for(size_type j = 0; j < per_flame_length / grid.per_FL(); ++j) {
          temp_global_solution_vector[j] = grid.global_solution_vector[i];
        }
      }
    }
  } else if(per_flame_length < grid.per_FL()) {
    if(grid.per_FL() % per_flame_length != 0) {
      std::cout << "cant refine to asked amount" << std::endl;
    } else {
      global_solution_vector_type temp_global_solution_vector;
      temp_global_solution_vector.resize(grid.domaine_length() * per_flame_length);
      for(size_type i = 0; i < temp_global_solution_vector.size(); ++i) {
        temp_global_solution_vector[i] = grid.global_solution_vector[i * grid.per_FL() / per_flame_length];
      }
      temp_global_solution_vector[0] = grid.global_solution_vector[0];
      temp_global_solution_vector[temp_global_solution_vector.size()-1] = grid.global_solution_vector[grid.number_of_cells()-1];
    }
  } else {
      std::cout << "Something went wrong in refinement" << std::endl;
  }
  grid.global_solution_vector = temp_global_solution_vector;
}

#endif //#ifndef SOLVER_H
