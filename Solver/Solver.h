#ifndef SOLVER_H
#define SOLVER_H

#include <math.h>
#include <vector>
#include <chrono>
#include<cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include "Explicit_Marching.h"
#include "../Usefull_Headers/Gnuplot_Primitive_Variables.h"
#include "../Serialization/Serialization_Eigen.h"
#include "../Serialization/Serialize.h"


template <typename global_solution_vector_type>
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
  Solver(global_solution_vector_type initial_solution_in, double Pr, double Le,
        double Q, double theta, double mf, double Lambda_in, double Lambda_min_in, double Lambda_max_in, double gamma,
        double number_of_cells, double CFL, double x_max, double x_min, double final_time_in,
        double frames_in, std::string filename_in) :
          global_solution_vector(initial_solution_in), lambda(Lambda_in), lambda_min(Lambda_min_in),
          lambda_max(Lambda_max_in), final_time(final_time_in), frames(frames_in), filename(filename_in),
          explicit_march(Explicit_Marching<global_solution_vector_type>(Pr, Le, Q, theta, mf, gamma,
          number_of_cells, CFL, (x_max - x_min)/number_of_cells)) {
    plot<global_solution_vector_type>(filename + std::to_string(static_cast<int>(Lambda_in)) + "_0",
                                      global_solution_vector, (x_max - x_min)/number_of_cells);
    old_position = flame_position_algorithm(initial_solution_in, gamma);
    std::cout << "//////////////////////" << std::endl;
    std::cout << "Solver, Lambda = " << lambda << std::endl;
    std::cout << "//////////////////////" << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Runs the simulation. Outputs frames in a folder "../Movie".
  void solve();

  template<typename Archive>
  void serialize(Archive& archive) {
    archive(global_solution_vector, current_time, lambda, lambda_min,
            lambda_max, final_time, frames, filename, explicit_march, old_position, curr_frame);
  }

 private:
   global_solution_vector_type global_solution_vector;
   double current_time = 0.0;
   double lambda, lambda_min, lambda_max;
   double final_time;
   int frames;
   std::string filename;
   Explicit_Marching<global_solution_vector_type> explicit_march;
   int old_position;
   int curr_frame = 0;

   void bisection_lambda(double& lambda_min, double& lambda_max, double& lambda_run, int& sum_initial, int& sum);

   int flame_position_algorithm(global_solution_vector_type global_solution_vector, double gamma);
};


///////////////////////////////////////////////////////////////////////////////
// Solver
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type>
void Solver<global_solution_vector_type>::solve() {
  double time_per_frame = final_time / frames;

  while (curr_frame < frames){
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Time = " <<  current_time << std::endl;
    explicit_march.timemarch(time_per_frame, global_solution_vector, lambda);
    current_time += time_per_frame;
    plot<global_solution_vector_type>(filename + std::to_string(static_cast<int>(lambda)) + "_" + std::to_string(static_cast<int>(curr_frame)+1), global_solution_vector, explicit_march.get_dx());
    curr_frame++;
    serialize_to_file(*this, filename);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
  }

  int position = 0;
  auto var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[position], 1.4);
  while (var_vec.rho() > 0.5) {
  // std::cout << global_solution_vector[i][0] << std::endl;
  ++position;
  var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[position], 1.4);
  }
  bisection_lambda(lambda_min, lambda_max, lambda, old_position, position);
}

template <typename global_solution_vector_type>
void Solver<global_solution_vector_type>::bisection_lambda(double& lambda_min, double& lambda_max, double& lambda_run, int& sum_initial, int& sum) {
  if (sum > sum_initial){
    lambda_min = lambda_run;
    lambda_run = (lambda_min + lambda_max)*0.5;
  } else {
    lambda_max = lambda_run;
    lambda_run = (lambda_min + lambda_max)*0.5;
  }
}

template <typename global_solution_vector_type>
int Solver<global_solution_vector_type>::flame_position_algorithm(global_solution_vector_type global_solution_vector, double gamma) {
  int i = 0;
  auto var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[0], gamma);
  while (var_vec.rho() > 0.5) {
  // std::cout << global_solution_vector[i][0] << std::endl;
  ++i;
  var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i], gamma);
  }
  std::cout << "initial:" << i << std::endl;
  return i;
}

#endif //#ifndef SOLVER_H
