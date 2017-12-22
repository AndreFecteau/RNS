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
          global_solution_vector(initial_solution_in), final_time(final_time_in), frames(frames_in),
          filename(filename_in),
          march(Marching<global_solution_vector_type>(Pr, Le, Q, theta, mf, Lambda_in, gamma,
          number_of_cells, CFL, (x_max - x_min)/number_of_cells)) {
    plot<global_solution_vector_type>(filename + std::to_string(static_cast<int>(Lambda_in)) + "_0",
                                      global_solution_vector, (x_max - x_min)/number_of_cells);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Runs the simulation. Outputs frames in a folder "../Movie".
  double solve();

 private:
   global_solution_vector_type global_solution_vector;
   double current_time = 0.0;
   const double final_time;
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
    plot<global_solution_vector_type>(filename + std::to_string(static_cast<int>(march.get_Lambda())) + "_" + std::to_string(static_cast<int>(i)+1), global_solution_vector, march.get_dx());
  }

  int i = 0;
  auto var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[0], 1.4);
  while (var_vec.rho() > 0.5) {
  // std::cout << global_solution_vector[i][0] << std::endl;
  ++i;
  var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i], 1.4);
  }
  std::cout << "i:" << i << std::endl;
  return i;
}

#endif //#ifndef SOLVER_H
