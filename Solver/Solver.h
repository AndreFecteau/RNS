#ifndef SOLVER_H
#define SOLVER_H

#include <chrono>
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>

#include "../Serialization/Serialization_Eigen.h"
#include "../Serialization/Serialize.h"
#include "../Gnuplot_RNS/Gnuplot_Primitive_Variables.h"


template <typename global_solution_vector_type, typename matrix_type>
class Solver {

using solution_vector_type = typename global_solution_vector_type::value_type;

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
  Solver(global_solution_vector_type initial_solution_in,
         std::string filename_in) : global_solution_vector(initial_solution_in),
         filename(filename_in), initial_solution((initial_solution_in)) {
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Runs the simulation. Outputs frames in a folder "../Movie".
  /// \param marching The type of marching method passed to the solver.
  ///        ex:Implicit or Explicit.
  /// \param target_residual Residual where the solver will stop.
  template <typename marching_type>
  bool solve(marching_type marching, double target_residual, double frame_time, double gamma, double lambda);

  template<typename Archive>
  void serialize(Archive& archive) {
    archive(global_solution_vector, current_time, filename, current_frame);
  }

 private:
  global_solution_vector_type global_solution_vector;
  std::string filename;
  global_solution_vector_type initial_solution;
  double current_time = 0.0;
  int current_frame = 0;

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  int flame_position_algorithm(double gamma);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  template <typename marching_type>
  void manufactured_solution_residual(marching_type march);
};

///////////////////////////////////////////////////////////////////////////////
// Solver
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
template <typename marching_type>
bool Solver<global_solution_vector_type, matrix_type>::solve(marching_type march,
                                                             double target_residual,
                                                             double frame_time,
                                                             double gamma,
                                                             double lambda) {

  std::cout << "//////////////////////" << std::endl;
  std::cout << "Lambda = " << lambda << std::endl;
  std::cout << "//////////////////////" << std::endl;

  int old_position;
  double residual = std::numeric_limits<double>::max();
  (void)residual;
  (void)target_residual;
  int i = 0;

  // while (i < 5){
  while (residual > target_residual){
    old_position = flame_position_algorithm(gamma);
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Time = " <<  current_time << std::endl;
    residual = march.timemarch(frame_time, global_solution_vector, lambda);
    current_time += frame_time;
    plot<global_solution_vector_type>(filename + std::to_string(static_cast<int>(current_frame)+1), global_solution_vector, march.get_dx());
    current_frame++;
    serialize_to_file(*this, filename);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
    ++i;
  }

#if defined(MANUFACTURED)
  manufactured_solution_residual<marching_type>(march);
#endif
  int position = 0;
  auto var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[position], 1.4);
  while (var_vec.rho() > 0.5) {
  ++position;
  var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[position], 1.4);
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
template <typename global_solution_vector_type, typename matrix_type>
template <typename marching_type>
void Solver<global_solution_vector_type, matrix_type>::manufactured_solution_residual(marching_type march){
  for (size_t i = 0; i < initial_solution.size(); ++i) {
   double x = (i+0.5)*march.get_dx();
   initial_solution[i] << (-0.45*tanh(4.0 * x - 10.0) + 0.55),
                         (-0.45*tanh(4.0 * x - 10.0) + 0.55)*(4.5*tanh(4.0 * x - 10.0) + 5.5),
                         2.0*tanh(4.0*x - 10.0) + 70000,
                         (-0.45*tanh(4.0 * x - 10.0) + 0.55)*(-0.5*tanh(x - 8.0/4.0) + 0.5);
  }
  double convergence = 0.0;
  for(size_t i = 0; i < initial_solution.size(); ++i) {
   for(int j = 0; j < 4; ++j){
     convergence += std::pow(std::fabs(initial_solution[i][j]-global_solution_vector[i][j]),2)*march.get_dx();
   }
  }

  convergence = std::sqrt(convergence);

  std::ofstream gnu_input_file;
  gnu_input_file.open("Convergence_Plot.dat", std::ios_base::app);
  gnu_input_file << global_solution_vector.size() << " " <<  convergence << " " << current_time << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
// Flame Position Algorithm
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
int Solver<global_solution_vector_type, matrix_type>::flame_position_algorithm(double gamma) {
  int i = 0;
  auto var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[0], gamma);
  while (var_vec.rho() > 0.5) {
  ++i;
  var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i], gamma);
  }
  return i;
}

#endif //#ifndef SOLVER_H
