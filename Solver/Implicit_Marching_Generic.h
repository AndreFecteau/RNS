#ifndef IMPLICIT_MARCHING_H
#define IMPLICIT_MARCHING_H

#include <omp.h>
#include <math.h>
#include <limits>
#include <vector>
#include "Eigen/Core"
#include "Eigen/Dense"
// #include "../Implicit_Flux_and_Sources/Implicit_Euler.h"
#include "../Implicit_Flux_and_Sources/Variable_Implicit_Scheme.h"
// #include "../Usefull_Headers/Block_Triagonal_Matrix_Inverse.h"
#include "../Matrix_Inverse/Gaussian_Block_Triagonal_Matrix_Inverse.h"
#include "../Usefull_Headers/Variable_Vector_Isolator.h"


template <typename global_solution_vector_type, typename matrix_type>
class Implicit_Marching {

  /////////////////////////////////////////////////////////////////////////
  /// \brief type for individual cell solution vector.
  using solution_vector_type = typename global_solution_vector_type::value_type;

 public:
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  Implicit_Marching() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  Implicit_Marching(const Implicit_Marching&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  Implicit_Marching(Implicit_Marching&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  Implicit_Marching& operator=(const Implicit_Marching&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  Implicit_Marching& operator=(Implicit_Marching&&) = default;

  Implicit_Marching(double Pr_in, double Le_in, double Q_in, double theta_in, double mf_in,
           double gamma_in, double number_of_cells_in,  double CFL_in, double dx_in,
           double Theta_in, double zeta_in) :
           Pr(Pr_in), Le(Le_in), Q(Q_in), theta(theta_in), mf(mf_in), gamma(gamma_in),
           CFL(CFL_in), number_of_cells(number_of_cells_in), dx(dx_in), Theta(Theta_in),
          zeta(zeta_in) {}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Execute step in time.
  /// \param time_frame Time to stop the time marching.
  /// \param global_solution_vector vector containing cell states from all the cells.
  double timemarch(double time_frame, global_solution_vector_type &global_solution_vector, double lambda);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  double get_dx() {
    return dx;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  template<typename Archive>
  void serialize(Archive& archive) {
    archive(Pr, Le, Q, theta, mf, gamma, CFL, number_of_cells, dx, dt, zeta, Theta);
  }

 private:
  HLLE<global_solution_vector_type> hyperbolic_flux;
  double Pr;
  double Le;
  double Q;
  double theta;
  double mf;
  double gamma;
  double CFL;
  double Theta;
  double zeta;
  int number_of_cells;
  double dx;
  double dt = 0.0;
  double residual;
  global_solution_vector_type delta_global_solution_vector_past = global_solution_vector_type(number_of_cells - 2,
  solution_vector_type::Zero());;


  /////////////////////////////////////////////////////////////////////////
  /// \brief Calculates maximum stable timestep.
  /// \param global_solution_vector vector containing cell states from all the cells.
  double calculate_dt(const global_solution_vector_type &global_solution_vector);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Calculates wavespeed for timestep control.
  /// \param global_solution_vector vector containing cell states from all the cells.
  double lambda_eigenvalue(const global_solution_vector_type &global_solution_vector);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Calculates the variable for timestep control from second order
  ///        derivatives.
  /// \param global_solution_vector vector containing cell states from all the cells.
  double K_value(const global_solution_vector_type &global_solution_vector);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  // template <typename T> int sign(T val) {return (T(0) < val) - (val < T(0));}

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  double squaredNorm(solution_vector_type solution_vector){
    return sqrt(solution_vector[0]*solution_vector[0] +
                solution_vector[1]*solution_vector[1] +
                solution_vector[2]*solution_vector[2] +
                solution_vector[3]*solution_vector[3]);
  }

  solution_vector_type numerical_dissipation(const global_solution_vector_type &global_solution_vector, const int i);

};


///////////////////////////////////////////////////////////////////////////////
// TimeMarch
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
double Implicit_Marching<global_solution_vector_type, matrix_type>::timemarch(double time_frame,
                             global_solution_vector_type &global_solution_vector, double Lambda) {
  double current_time = 0.0;
  std::vector<matrix_type>    mid(global_solution_vector.size()-2);
  std::vector<matrix_type>    top(global_solution_vector.size()-2);
  std::vector<matrix_type>    bot(global_solution_vector.size()-2);
  global_solution_vector_type rhs(global_solution_vector.size()-2);
  global_solution_vector_type delta_global_solution_vector;
#pragma omp parallel
  {
  while (current_time < time_frame){

#pragma omp single
{
  dt = calculate_dt(global_solution_vector);
  if(current_time + dt > time_frame) {
    dt = time_frame - current_time;
  }
}
#pragma omp single
{
  int i = 1;
  mid[i-1] = create_mid_band_matrix<solution_vector_type, matrix_type>
                                 (global_solution_vector[i-1], global_solution_vector[i],
                                  global_solution_vector[i+1], gamma, Pr, Le, Q, Lambda,
                                  theta, dx, dt, zeta, Theta);
  bot[i-1] = create_bot_band_matrix<solution_vector_type, matrix_type>
                                 (global_solution_vector[i-1], global_solution_vector[i],
                                  global_solution_vector[i+1], gamma, Pr, Le, Q, Lambda,
                                  theta, dx, dt, zeta, Theta);
  top[i-1] = create_top_band_matrix<solution_vector_type, matrix_type>
                                 (global_solution_vector[i-1], global_solution_vector[i],
                                  global_solution_vector[i+1], gamma, Pr, Le, Q, Lambda,
                                  theta, dx, dt, zeta, Theta);
  rhs[i-1] = create_rhs_vector<solution_vector_type, matrix_type>
                                 (global_solution_vector[i-1], global_solution_vector[i],
                                  global_solution_vector[i+1], gamma, Pr, Le, Q, Lambda,
                                  theta, dx, dt, zeta, Theta, delta_global_solution_vector_past[i-1]);
}
#pragma omp for
  for(size_t i = 2; i < global_solution_vector.size()-2; ++i) {
    mid[i-1] = create_mid_band_matrix<solution_vector_type, matrix_type>
                                   (global_solution_vector[i-1], global_solution_vector[i],
                                    global_solution_vector[i+1], gamma, Pr, Le, Q, Lambda,
                                    theta, dx, dt, zeta, Theta);
    bot[i-1] = create_bot_band_matrix<solution_vector_type, matrix_type>
                                   (global_solution_vector[i-1], global_solution_vector[i],
                                    global_solution_vector[i+1], gamma, Pr, Le, Q, Lambda,
                                    theta, dx, dt, zeta, Theta);
    top[i-1] = create_top_band_matrix<solution_vector_type, matrix_type>
                                   (global_solution_vector[i-1], global_solution_vector[i],
                                    global_solution_vector[i+1], gamma, Pr, Le, Q, Lambda,
                                    theta, dx, dt, zeta, Theta);
    rhs[i-1] = create_rhs_vector<solution_vector_type, matrix_type>
                                   (global_solution_vector[i-1], global_solution_vector[i],
                                    global_solution_vector[i+1], gamma, Pr, Le, Q, Lambda,
                                    theta, dx, dt, zeta, Theta, delta_global_solution_vector_past[i-1]);
    // rhs[i-1] +=0.1* numerical_dissipation(global_solution_vector, i);
    // if (mid[i-1].determinant() < (bot[i-1].determinant() + top[i-1].determinant())) {
    //   std::cout << "Matrix inverse will not work." << std::endl;
    // }
  }
#pragma omp single
{
  int i = global_solution_vector.size()-2;
  mid[i-1] = create_mid_band_matrix<solution_vector_type, matrix_type>
                                 (global_solution_vector[i-1], global_solution_vector[i],
                                  global_solution_vector[i+1], gamma, Pr, Le, Q, Lambda,
                                  theta, dx, dt, zeta, Theta);
  bot[i-1] = create_bot_band_matrix<solution_vector_type, matrix_type>
                                 (global_solution_vector[i-1], global_solution_vector[i],
                                  global_solution_vector[i+1], gamma, Pr, Le, Q, Lambda,
                                  theta, dx, dt, zeta, Theta);
  top[i-1] = create_top_band_matrix<solution_vector_type, matrix_type>
                                 (global_solution_vector[i-1], global_solution_vector[i],
                                  global_solution_vector[i+1], gamma, Pr, Le, Q, Lambda,
                                  theta, dx, dt, zeta, Theta);
  rhs[i-1] = create_rhs_vector<solution_vector_type, matrix_type>
                                 (global_solution_vector[i-1], global_solution_vector[i],
                                  global_solution_vector[i+1], gamma, Pr, Le, Q, Lambda,
                                  theta, dx, dt, zeta, Theta, delta_global_solution_vector_past[i-1]);
}

#pragma omp single
  mid[global_solution_vector.size()-3] +=
  create_top_band_matrix<solution_vector_type, matrix_type>
                        (global_solution_vector[global_solution_vector.size()-3],
                         global_solution_vector[global_solution_vector.size()-2],
                         global_solution_vector[global_solution_vector.size()-1],
                         gamma, Pr, Le, Q, Lambda, theta, dx, dt, zeta, Theta);

#pragma omp single
  delta_global_solution_vector = block_triagonal_matrix_inverse<matrix_type, solution_vector_type>(mid, top, bot, rhs);

#pragma omp for
  for (int i = 1; i < number_of_cells-1; ++i) {
    global_solution_vector[i] += delta_global_solution_vector[i-1];
  }
  #pragma omp single
  global_solution_vector[global_solution_vector.size()-1] = global_solution_vector[global_solution_vector.size()-2];

#pragma omp single
  current_time += dt;
  }
  residual = 0.0;
#pragma omp for reduction(+: residual)
  for (int i = 1; i < number_of_cells-1; ++i) {
    residual += delta_global_solution_vector[i-1].squaredNorm() * dx / dt;
  }
  }
  std::cout << "residual: " << residual << std::endl;
  return residual;
}

///////////////////////////////////////////////////////////////////////////////
// Calculate dt;
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
double Implicit_Marching<global_solution_vector_type, matrix_type>::calculate_dt(const global_solution_vector_type &global_solution_vector) {
  double dt1 = CFL * dx / lambda_eigenvalue(global_solution_vector);
  double dt2 = CFL * dx*dx / (K_value(global_solution_vector));
  return std::min(dt1, dt2);
}

///////////////////////////////////////////////////////////////////////////////
// WaveSpeed
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
double Implicit_Marching<global_solution_vector_type, matrix_type>::lambda_eigenvalue(const global_solution_vector_type &global_solution_vector){
  double lambda = 0.0;
  for (int i = 0; i < number_of_cells; ++i) {
    Variable_Vector_Isolator<solution_vector_type> var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i], gamma);
    if (lambda < std::fabs(var_vec.u()) + sqrt(gamma * var_vec.p()/var_vec.rho())) {
      lambda = std::fabs(var_vec.u()) + sqrt(gamma*var_vec.p()/var_vec.rho());
    }
  }
  return lambda;
}

///////////////////////////////////////////////////////////////////////////////
// K value
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
double Implicit_Marching<global_solution_vector_type, matrix_type>::K_value(const global_solution_vector_type &global_solution_vector) {
  double min_rho = std::numeric_limits<double>::max();
  for (int i = 0; i < number_of_cells; ++i) {
    Variable_Vector_Isolator<solution_vector_type> var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i], gamma);
    if (min_rho > var_vec.rho()) {
      min_rho = var_vec.rho();
    }
  }
  return 4.0 * std::max(std::max(Pr,Le), gamma / (gamma - 1.0)) / min_rho;
}

///////////////////////////////////////////////////////////////////////////////
// Add Numerical Dissipation
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
typename global_solution_vector_type::value_type Implicit_Marching<global_solution_vector_type, matrix_type>::numerical_dissipation(const global_solution_vector_type &global_solution_vector, const int i) {
            return -0.2/(1.0+zeta)/8.0*(global_solution_vector[i+2]-4.0*global_solution_vector[i+1] + 6.0*global_solution_vector[i] - 4.0*global_solution_vector[i-1]+global_solution_vector[i-2]);
}

#endif //#ifndef IMPLICIT_MARCHING_H
