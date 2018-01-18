#ifndef SOLVER_H
#define SOLVER_H

#include <math.h>
#include <vector>
#include <chrono>
#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include "../Usefull_Headers/Gnuplot_Primitive_Variables.h"
#include "../Serialization/Serialization_Eigen.h"
#include "../Serialization/Serialize.h"


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
         std::string filename_in) :
          global_solution_vector(initial_solution_in), lambda(Lambda_in), filename(filename_in)
  {
    std::cout << "//////////////////////" << std::endl;
    std::cout << "Solver, Lambda = " << lambda << std::endl;
    std::cout << "//////////////////////" << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Runs the simulation. Outputs frames in a folder "../Movie".
  template <typename marching_type>
  bool solve(marching_type marching, double target_residual, double frame_time, double gamma);

  template<typename Archive>
  void serialize(Archive& archive) {
    archive(global_solution_vector, current_time, lambda,
            filename, old_position, current_frame);
  }

 private:
   global_solution_vector_type global_solution_vector;
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
  double residual = 100000.0;
  old_position = flame_position_algorithm(global_solution_vector, gamma);
  // while (residual > target_residual){
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Time = " <<  current_time << std::endl;
    residual = march.timemarch(frame_time, global_solution_vector, lambda);
    current_time += frame_time;
    plot<global_solution_vector_type>(filename + std::to_string(static_cast<int>(lambda)) + "_" + std::to_string(static_cast<int>(current_frame)+1), global_solution_vector, march.get_dx());
    current_frame++;
    serialize_to_file(*this, filename);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
  // }

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
