#ifndef SOLVER_H
#define SOLVER_H

#include <math.h>
#include <chrono>
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>

#include "../Serialization/Serialization_Eigen.h"
#include "../Serialization/Serialize.h"
#include "../Usefull_Headers/Math.h"
#include "../Gnuplot_RNS/Gnuplot_Primitive_Variables.h"
#include "../Gnuplot_RNS/Gnuplot_Primitive_Variables_Reduced.h"
#include "../Gnuplot_RNS/Gnuplot_Variables.h"


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
         scalar_type CFL_in, scalar_type Theta_in, scalar_type zeta_in,
         std::string filename_in, scalar_type flame_location_in,
         scalar_type dissipation_magnitude_in = 0.0) :
         flow(flow_properties_in), grid(grid_in), frame_time(frame_time_in),
         CFL(CFL_in),
         time_stepping(time_stepping_type(Theta_in, zeta_in, dissipation_magnitude_in)),
         filename(filename_in), flame_location(flame_location_in) {}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Runs the simulation. Outputs frames in a folder "../Movie".
  /// \param marching The type of marching method passed to the solver.
  ///        ex:Implicit or Explicit.
  /// \param number_of_frames number of frames where the solver will stop.
  bool solve(size_type number_of_frames);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  flow_properties_type get_flow() const {return flow;}

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  grid_type get_grid() const {return grid;}

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  scalar_type get_lambda() const {return flow.lambda;}

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  void recenter_solution(const scalar_type &flame_location);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  void add_space_in_back(scalar_type space);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  void add_space_in_front(scalar_type space);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  void print_stats();

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  void change_filename(std::string new_filename) {
    filename = new_filename;
    std::cout << "filename changed to: " << filename << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  void reset_frame_number() {
    global_current_frame = 0;
    std::cout << "Reset Frame Number" << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  void change_lambda(scalar_type lambda) {
    flow.lambda = lambda;
    std::cout << "lambda changed to: " << flow.lambda << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  void refine(size_type per_flame_length);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  void change_frame_time(scalar_type new_frame_time) {
    frame_time = new_frame_time;
    std::cout << "frame_time changed to: " << frame_time << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  void change_CFL(scalar_type new_CFL) {
    CFL = new_CFL;
    std::cout << "CFL changed to: " << CFL << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  void change_flame_location(scalar_type new_flame_location) {
    flame_location = new_flame_location;
    std::cout << "flame_location changed to: " << flame_location << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  void plot_limiter();

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  void plot_global_solution_vector(std::string plot_name, scalar_type number_of_cells_in_output) {
     plot_reduced<grid_type>(plot_name, grid, flow, number_of_cells_in_output);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  template<typename Archive>
  void serialize(Archive& archive) {
    archive(flow, grid, frame_time, CFL, time_stepping, filename, current_time, global_current_frame, flame_location);
  }

 private:
  flow_properties_type flow;
  grid_type grid;
  scalar_type frame_time;
  scalar_type CFL;
  time_stepping_type time_stepping;
  global_solution_vector_type global_solution_vector_backup;
  scalar_type current_time = 0.0;
  size_type global_current_frame = 0;
  std::string filename;
  scalar_type flame_location;

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  size_type flame_position_algorithm();

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  void manufactured_solution_residual();

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  bool flame_left_domaine(const size_type &position, scalar_type &frame_time_temp);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  bool solution_is_unstable(const scalar_type &residual, scalar_type &frame_CFL, scalar_type &frame_time_temp);

};

///////////////////////////////////////////////////////////////////////////////
// Solver
///////////////////////////////////////////////////////////////////////////////
template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
bool Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
solve(const size_type number_of_frames) {

  std::cout << "//////////////////////" << std::endl;
  std::cout << "Solver"<< std::endl;
  std::cout << "//////////////////////" << std::endl;

  scalar_type residual = std::numeric_limits<scalar_type>::max();
  size_type frame_counter = 0;
  scalar_type frame_time_temp = frame_time;
  scalar_type frame_CFL = CFL;
  size_type position = 0;
  plot<grid_type>(filename + std::to_string(static_cast<size_type>(global_current_frame)), grid.global_solution_vector,grid.dx());
  while (frame_counter < number_of_frames) {
    auto start = std::chrono::high_resolution_clock::now();

    position = 0;
    global_solution_vector_backup = grid.global_solution_vector;

    residual = time_stepping.template timemarch<flux_type>(flow, grid, frame_CFL, frame_time_temp);

    position = flame_position_algorithm();
    if (solution_is_unstable(residual, frame_CFL, frame_time_temp)) {
#if defined(RECENTER_FLAME)
    } else if (flame_left_domaine(position, frame_time_temp)) {
#endif
    } else {
      current_time += frame_time_temp;
      ++global_current_frame;
      ++frame_counter;

      plot<grid_type>(filename + std::to_string(static_cast<size_type>(global_current_frame)), grid.global_solution_vector,grid.dx());
      serialize_to_file(*this, filename + std::to_string(static_cast<size_type>(global_current_frame)));

      auto finish = std::chrono::high_resolution_clock::now();
      std::chrono::duration<scalar_type> elapsed = finish - start;
      std::cout << "Frame: " << global_current_frame << " Frame_time = " <<  frame_time_temp << " Residual: " << residual << std::endl;
      std::cout << "position: " << position - flame_location*grid.per_FL();
      std::cout << " Elapsed time: " << elapsed.count() << "\n";

#if defined(RECENTER_FLAME)
      recenter_solution(flame_location);
#endif
    frame_time_temp = frame_time;
    frame_CFL = CFL;
    }
#if defined(MANUFACTURED)
    manufactured_solution_residual();
#endif
  }
  const size_type h = std::max(static_cast<int>(5*grid.per_FL()), static_cast<int>(5*grid.per_FL()-(position - flame_location*grid.per_FL())));
  if(grid.global_solution_vector[h][1] / grid.global_solution_vector[h][0] < 1.0) {
    return 0;
  } else {
    return 1;
  }
}

///////////////////////////////////////////////////////////////////////////////
// Manufactured Solution Residual
///////////////////////////////////////////////////////////////////////////////
template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
void Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
manufactured_solution_residual() {
  global_solution_vector_type initial_solution = global_solution_vector_type(grid.number_of_cells(),
  solution_vector_type::Zero());
  for (size_t i = 0; i < grid.number_of_cells(); ++i) {
   scalar_type x = (i+0.5)*grid.dx();
    initial_solution[i] << (-0.45*tanh(4.0 * x - 10.0) + 0.55),
                          (-0.45*tanh(4.0 * x - 10.0) + 0.55)*(4.5*tanh(4.0 * x - 10.0) + 5.5),
                          2.0*tanh(4.0*x - 10.0) + 70000,
                          (-0.45*tanh(4.0 * x - 10.0) + 0.55)*(-0.5*tanh(x - 8.0/4.0) + 0.5);
  }
  scalar_type convergence = 0.0;
  for(size_t i = 0; i < initial_solution.size(); ++i) {
   for(size_type j = 0; j < 4; ++j){
     convergence += std::pow(std::fabs(initial_solution[i][j]-grid.global_solution_vector[i][j]),2)*grid.dx();
   }
  }

  convergence = std::sqrt(convergence);
  std::cout << "convergence: " << convergence << std::endl;
  std::ofstream gnu_input_file;
  gnu_input_file.open("Convergence_Plot.dat", std::ios_base::app);
  gnu_input_file << grid.number_of_cells() << " " <<  convergence << " " << current_time << std::endl;
  gnu_input_file.close();
}


///////////////////////////////////////////////////////////////////////////////
// Flame Position Algorithm
///////////////////////////////////////////////////////////////////////////////
template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
typename grid_type::size_type Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
flame_position_algorithm() {
  size_type position = 0;
  auto var_vec = Variable_Vector_Isolator<grid_type>(grid.global_solution_vector[0], flow.gamma);
  while (var_vec.Y() > 0.5) {
  ++position;
  var_vec = Variable_Vector_Isolator<grid_type>(grid.global_solution_vector[position], flow.gamma);
  }
  return position;
}

///////////////////////////////////////////////////////////////////////////////
// Recenter Flame
///////////////////////////////////////////////////////////////////////////////
template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
void Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
recenter_solution(const scalar_type &flame_location) {
  auto global_solution_vector_temp = grid.global_solution_vector;
  size_type position = 0;
  auto var_vec = Variable_Vector_Isolator<grid_type>(grid.global_solution_vector[position], 1.4);
  while (var_vec.Y() > 0.5) {
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

///////////////////////////////////////////////////////////////////////////////
// add_space_in_back
///////////////////////////////////////////////////////////////////////////////
template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
void Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
add_space_in_back(scalar_type space) {
  grid.global_solution_vector.resize(grid.number_of_cells() + space/grid.dx());
  grid.x_max += space;
  for(size_t i = grid.number_of_cells() - space/grid.dx(); i < grid.number_of_cells(); ++i) {
    grid.global_solution_vector[i] = grid.global_solution_vector[i-1];
  }
  std::cout << "Added space in domaine, New domaine size: " << grid.domaine_length() << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// add_space_in_front
///////////////////////////////////////////////////////////////////////////////
template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
void Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
add_space_in_front(scalar_type space) {
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

///////////////////////////////////////////////////////////////////////////////
// print_stats
///////////////////////////////////////////////////////////////////////////////
template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
void Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
print_stats() {
    std::cout << "Pr: " << flow.Pr << "\nLe: " << flow.Le << "\ngamma: " << flow.gamma << std::endl;
    std::cout << "mf: " << flow.mf << "\nQ_low_mach: " << flow.Q_low_mach << "\ntheta_low_mach: " << flow.theta_low_mach << std::endl;
    std::cout << "T_ignition: " << flow.T_ignition_scalar << "\nlambda: " << flow.lambda  << std::endl;
    std::cout << "x_min: " << grid.x_min << "\nx_max: " << grid.x_max << "\nnumber_of_cells: " << grid.number_of_cells() << std::endl;
    std::cout <<"Q: " << flow.Q() << "\ntheta: " << flow.theta() << std::endl;

}

///////////////////////////////////////////////////////////////////////////////
// refine
///////////////////////////////////////////////////////////////////////////////
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
          temp_global_solution_vector[i*per_flame_length / grid.per_FL()+j] = grid.global_solution_vector[i];
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
  std::cout << "Number_cells: " << grid.number_of_cells() << std::endl;
  plot<grid_type>(filename + std::to_string(static_cast<size_type>(2000)), grid.global_solution_vector, grid.dx());
}

///////////////////////////////////////////////////////////////////////////////
// plot_limiter
///////////////////////////////////////////////////////////////////////////////
template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
void Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
plot_limiter() {
  global_solution_vector_type min_mod(grid.number_of_cells());
  global_solution_vector_type van_albada(grid.number_of_cells());
  global_solution_vector_type solution_vector_rho(grid.number_of_cells());
  for (size_type j = 0; j < grid.number_of_cells(); ++j){
    auto var_vec_m = Variable_Vector_Isolator<grid_type>(grid.global_solution_vector[std::max(0,static_cast<int>(j)-1)], flow.gamma);
    auto var_vec = Variable_Vector_Isolator<grid_type>(grid.global_solution_vector[j], flow.gamma);
    auto var_vec_p = Variable_Vector_Isolator<grid_type>(grid.global_solution_vector[std::min(j+1, grid.number_of_cells() -1)], flow.gamma);
    const solution_vector_type U = var_vec.w();
    const solution_vector_type Ul =var_vec_m.w();
    const solution_vector_type Ur = var_vec_p.w();

    const solution_vector_type a = (U - Ul) / grid.dx();
    const solution_vector_type b = (Ur - U) / grid.dx();
    for(size_type i = 0; i < 4; ++i) {
      min_mod[j][i] = math::sign(a[i])*std::max(0.0,static_cast<scalar_type>(std::min(fabs(a[i]), math::sign(a[i])*b[i])));
    }

    double epsilon = 1.0e-6;
    van_albada[j] =  a.array() * b.array() * (a.array() + b.array()) / (a.array() * a.array() + b.array() * b.array() + epsilon);
    for (int i = 0; i < 4; ++i) {
      if (a[i] / b[i] <= 0.0 || b[i] == 0) {
        van_albada[j][i] = 0.0;
      }
    }
   solution_vector_rho[j] = grid.global_solution_vector[j] / grid.global_solution_vector[j][0];
  }
  plot_v<global_solution_vector_type>("min_mod", min_mod, grid.dx());
  plot_v<global_solution_vector_type>("van_albada", van_albada, grid.dx());
  plot_v<global_solution_vector_type>("solution_vector", grid.global_solution_vector, grid.dx());
  plot_v<global_solution_vector_type>("solution_vector_rho", solution_vector_rho, grid.dx());
}

///////////////////////////////////////////////////////////////////////////////
// flame_left_domaine
///////////////////////////////////////////////////////////////////////////////
template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
bool Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
flame_left_domaine(const size_type &position, scalar_type &frame_time_temp) {
  bool check = position < 100*grid.per_FL() || position > grid.number_of_cells()*0.9;
  if(check) {
    grid.global_solution_vector = global_solution_vector_backup;
    frame_time_temp *= 0.5;
    std::cout << ":" << std::flush;
  }
  return check;
}

///////////////////////////////////////////////////////////////////////////////
// solution_is_unstable
///////////////////////////////////////////////////////////////////////////////
template <typename flow_properties_type, typename grid_type, typename flux_type, typename time_stepping_type>
bool Solver<flow_properties_type, grid_type, flux_type, time_stepping_type>::
solution_is_unstable(const scalar_type &residual, scalar_type &frame_CFL, scalar_type &frame_time_temp) {
  bool check = isnan(residual) || residual > 1e10;
  if(check){
    grid.global_solution_vector = global_solution_vector_backup;
    frame_CFL *= 0.5;
    frame_time_temp *= 0.5;
    std::cout << "." << std::flush;
  }
  return check;
}

#endif //#ifndef SOLVER_H
