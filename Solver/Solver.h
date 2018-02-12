#ifndef SOLVER_H
#define SOLVER_H

#include <math.h>
#include <vector>
#include <chrono>
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include "../Gnuplot_RNS/Gnuplot_Primitive_Variables.h"
#include "../Serialization/Serialization_Eigen.h"
#include "../Serialization/Serialize.h"
#include "../Usefull_Headers/Variable_Vector_Isolator.h"


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
  Solver(global_solution_vector_type initial_solution_in, double Lambda_in,
         std::string filename_in) : global_solution_vector(initial_solution_in),
         lambda(Lambda_in), filename(filename_in), initial_solution((initial_solution_in)) {
    std::cout << "//////////////////////" << std::endl;
    std::cout << "Solver, Lambda = " << lambda << std::endl;
    std::cout << "//////////////////////" << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Runs the simulation. Outputs frames in a folder "../Movie".
  /// \param marching The type of marching method passed to the solver.
  ///        ex:Implicit or Explicit.
  /// \param target_residual Residual where the solver will stop.
  template <typename marching_type>
  bool solve(marching_type marching, double target_residual, double frame_time, double gamma);

  template<typename Archive>
  void serialize(Archive& archive) {
    archive(global_solution_vector, current_time, lambda,
            filename, old_position, current_frame);
  }
  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
template <typename marching_type>
  void expand_solution_vector(int size, marching_type march) {
    global_solution_vector_type new_solution_vector(global_solution_vector.size()*size);
    for(size_t i = 0; i < global_solution_vector.size(); ++i) {
      solution_vector_type du;
      if (i < global_solution_vector.size() - 1){
        du = (global_solution_vector[i+1] - global_solution_vector[i])/size;
      }
      for (size_t j = 0; j < size; ++j) {
        new_solution_vector[i*size+j] = global_solution_vector[i] + du*j;
      }
    }
    // new_solution_vector[new_solution_vector.size()-1] = global_solution_vector[global_solution_vector.size()-1];
    global_solution_vector.resize(new_solution_vector.size());
    global_solution_vector = new_solution_vector;
    // plot<global_solution_vector_type>(filename + std::to_string(static_cast<int>(lambda)) + "_" + std::to_string(static_cast<int>(current_frame)+1), global_solution_vector, march.get_dx());

  }
 private:
   global_solution_vector_type global_solution_vector;
   global_solution_vector_type initial_solution;
   std::string filename;
   double lambda;
   double current_time = 0.0;
   int current_frame = 0;
   int old_position;

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
   int flame_position_algorithm(global_solution_vector_type global_solution_vector, double gamma);

};

///////////////////////////////////////////////////////////////////////////////
// Solver
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
template <typename marching_type>
bool Solver<global_solution_vector_type, matrix_type>::solve(marching_type march,
                                                             double target_residual,
                                                             double frame_time,
                                                             double gamma) {
  // plot<global_solution_vector_type>(filename + std::to_string(static_cast<int>(lambda)) + "_0" ), global_solution_vector, march.get_dx());
  double residual = 100000.0;
  int i = 0;
  // while (i < 10){
  while (residual > target_residual){
    old_position = flame_position_algorithm(global_solution_vector, gamma);
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

  for (size_t i = 0; i < initial_solution.size(); ++i) {
    double dx = march.get_dx();
    double x = (i+0.5)*dx;
    initial_solution[i] << (-0.45*tanh(4.0 * x - 10.0) + 0.55),
                          (-0.45*tanh(4.0 * x - 10.0) + 0.55)*(4.5*tanh(4.0 * x - 10.0) + 5.5),
                          2.0*tanh(4.0*x - 10.0) + 70000,
                          (-0.45*tanh(4.0 * x - 10.0) + 0.55)*(-0.5*tanh(x - 8.0/4.0) + 0.5);
  }
  double convergence = 0.0;
  for(size_t i = 0; i < initial_solution.size(); ++i) {
    for(int j = 0; j < 4; ++j){
      convergence += std::fabs(initial_solution[i][j]-global_solution_vector[i][j]);
    }
  }

  std::ofstream gnu_input_file;
  gnu_input_file.open("Convergence_Plot.dat", std::ios_base::app);
  gnu_input_file << global_solution_vector.size() << " " <<  convergence << " " << current_time << std::endl;

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

template <typename global_solution_vector_type, typename matrix_type>
int Solver<global_solution_vector_type, matrix_type>::flame_position_algorithm(global_solution_vector_type global_solution_vector, double gamma) {
  int i = 0;
  auto var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[0], gamma);
  while (var_vec.rho() > 0.5) {
  ++i;
  var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i], gamma);
  }
  return i;
}

#endif //#ifndef SOLVER_H
