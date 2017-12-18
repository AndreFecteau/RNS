#ifndef SOLVER_H
#define SOLVER_H

#include <math.h>
#include <vector>
#include "Marching.h"
#include "../Usefull_Headers/Gnuplot_Primitive_Variables.h"

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
        double Q, double theta, double mf, double Lambda_in, double gamma,
        double number_of_cells, double CFL, double x_max, double x_min, double final_time_in,
        double frames_in, std::string filename_in) :
          global_solution_vector(initial_solution_in), Lambda(Lambda_in),
          final_time(final_time_in), dx((x_max - x_min)/number_of_cells), frames(frames_in), filename(filename_in),
          march(Marching<global_solution_vector_type>(Pr, Le, Q, theta, mf, Lambda, gamma, number_of_cells, CFL, dx)) {
    plot<global_solution_vector_type>(filename + std::to_string(static_cast<int>(Lambda)) + "_0", global_solution_vector, dx);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Runs the simulation. Outputs frames in a folder "../Movie".
  double solve();

 private:
   global_solution_vector_type global_solution_vector;
   double current_time = 0.0;
   const double Lambda;
   const double final_time;
   const double dx;
   const int frames;
   std::string filename;
   Marching<global_solution_vector_type> march;
};

///////////////////////////////////////////////////////////////////////////////
// Solver
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type>
double Solver<global_solution_vector_type>::solve() {
  double time_per_frame = final_time / frames;
  for (int i = 0; i < frames; ++i){
    std::cout << "Time = " <<  current_time << std::endl;
    march.timemarch(time_per_frame, global_solution_vector);
    current_time += time_per_frame;
    plot<global_solution_vector_type>(filename + std::to_string(static_cast<int>(Lambda)) + "_" + std::to_string(static_cast<int>(i)+1), global_solution_vector, dx);
  }

  double sum = 0.0;
  for (size_t i = 0; i < global_solution_vector.size(); ++i){
    Variable_Vector_Isolator<solution_vector_type> var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i], 1.4);
    sum += var_vec.T();
  }
  return sum;

  //
  // double max_p = 0;
  // for (size_t i = 1000; i < global_solution_vector.size(); ++i) {
  //   Variable_Vector_Isolator<solution_vector_type> var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i], 1.4);
  //   if (max_p < var_vec.p()) {
  //     max_p = var_vec.p();
  //   }
  // }
  // return max_p;
}

#endif //#ifndef SOLVER_H
