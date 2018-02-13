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
// #include "../Implicit_Flux_and_Sources/Variable_Implicit_Scheme_HLLE.h"
// #include "../Usefull_Headers/Block_Triagonal_Matrix_Inverse.h"
#include "../Matrix_Inverse/Gaussian_Block_Triagonal_Matrix_Inverse.h"
#include "../Usefull_Headers/Variable_Vector_Isolator.h"


#include "../Explicit_Flux_and_Sources/HLLE.h"
#include "../Explicit_Flux_and_Sources/Centered_Difference.h"
#include "../Explicit_Flux_and_Sources/Sources.h"
#include "Explicit_Marching.h"


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
  double Pr, Le, Q, theta;
  double mf, gamma, CFL;
  int number_of_cells;
  double dx;
  double Theta, zeta;
  double dt = 0.0;
  double residual;
  int count = 0;


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

  solution_vector_type numerical_dissipation(const global_solution_vector_type &global_solution_vector, const int i, const double omega);

  solution_vector_type manufactured_residual(const double lambda, const int i);

  template <typename T>
  double Power(const T num, const int expo) {return pow(num, expo);}

  double Power(const char num, const double expo) {return exp(expo);}

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
  global_solution_vector_type delta_global_solution_vector = global_solution_vector_type(number_of_cells-2,
  solution_vector_type::Zero());;

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
#pragma omp for
  for(int i = 1; i < static_cast<int>(global_solution_vector.size()-1); ++i) {
    auto matrix_entries = Implicit_Matrix_Entries<global_solution_vector_type, matrix_type>
                                       (global_solution_vector[std::max(i-2,0)], global_solution_vector[i-1], global_solution_vector[i],
                                        global_solution_vector[i+1],global_solution_vector[std::min(i+2,number_of_cells-1)],
                                        delta_global_solution_vector[i-1],
                                        gamma, Pr, Le, Q, Lambda,
                                        theta, dx, dt, zeta, Theta);
    mid[i-1] = matrix_entries.mid_matrix();
    bot[i-1] = matrix_entries.bot_matrix();
    top[i-1] = matrix_entries.top_matrix();
    rhs[i-1] = matrix_entries.rhs_matrix();
    rhs[i-1] += manufactured_residual(Lambda, i)*dt;

    // rhs[i-1] += numerical_dissipation(global_solution_vector, i, 0.9);
    // rhs[i-1] += numerical_dissipation(global_solution_vector, i, 0.1);
    // rhs[i-1] += numerical_dissipation(global_solution_vector, i, 0.01);
  }

// Implicit Boundary Conditions
#pragma omp single
  mid[global_solution_vector.size()-3] += top[global_solution_vector.size()-3];

#pragma omp single
  delta_global_solution_vector = block_triagonal_matrix_inverse<matrix_type, solution_vector_type>(mid, top, bot, rhs);

#pragma omp for
  for (int i = 1; i < number_of_cells-1; ++i) {
    if(current_time == 0.0) {
      global_solution_vector[i] += delta_global_solution_vector[i-1]*(1);
    } else {
      global_solution_vector[i] += delta_global_solution_vector[i-1];
    }
  }

//Explicit Boundary Conditions
// #pragma omp single
//   {
//     global_solution_vector[global_solution_vector.size()-1] = global_solution_vector[global_solution_vector.size()-2];
//   }

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
  ++count;
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
typename global_solution_vector_type::value_type Implicit_Marching<global_solution_vector_type, matrix_type>::
numerical_dissipation(const global_solution_vector_type &global_solution_vector, const int i, const double omega) {
            return -omega/(1.0+zeta)/8.0*(global_solution_vector[std::min(i+2,number_of_cells-1)] -
                                        4.0*global_solution_vector[std::min(i+1,number_of_cells-1)] +
                                        6.0*global_solution_vector[i] -
                                        4.0*global_solution_vector[std::max(i-1,0)] +
                                        global_solution_vector[std::max(i-2,0)]);
}

///////////////////////////////////////////////////////////////////////////////
// Manufactured Solution
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
typename global_solution_vector_type::value_type Implicit_Marching<global_solution_vector_type, matrix_type>::manufactured_residual(const double lambda, const int i) {
    solution_vector_type temp;
    solution_vector_type man_sol;
    man_sol << 0,0,0,0;
    double x = dx*(i+0.5);
    ///////////////////////////////////////////////////////////////////////////////
    // Full
    ///////////////////////////////////////////////////////////////////////////////
  //   temp << (81*Power(1.0/cosh(10 - 4*x),2)*tanh(10 - 4*x))/5.,(Power(1.0/cosh(10 - 4*x),2)*
  //     (2947 - 769*gamma - 6*(-891 + 297*gamma + 1280*Pr)*tanh(10 - 4*x) + 2187*(-3 + gamma)*Power(tanh(10 - 4*x),2)))/40.,
  //  (36*Power(1.0/cosh(10 - 4*x),4)*(-100791541*gamma - 15972*Pr + 9801*(3*gamma - 4*Pr)*tanh(10 - 4*x) + 8019*(3*gamma - 4*Pr)*Power(tanh(10 - 4*x),2) +
  //        2187*(3*gamma - 4*Pr)*Power(tanh(10 - 4*x),3)))/Power(11 + 9*tanh(10 - 4*x),3) -
  //   (Power(1.0/cosh(10 - 4*x),2)*(-121*(11979 + 50389781*gamma) + (6147901762*gamma + 7986*(-297 + 640*Pr))*tanh(10 - 4*x) -
  //        54*(75686097*gamma - 121*(297 + 640*Pr))*Power(tanh(10 - 4*x),2) - 972*(-3267 + 387*gamma + 3520*Pr)*Power(tanh(10 - 4*x),3) +
  //        2187*(-297 + 1257*gamma - 1280*Pr)*Power(tanh(10 - 4*x),4) + 1062882*(-1 + gamma)*Power(tanh(10 - 4*x),5)))/(40.*Power(11 + 9*tanh(10 - 4*x),2))
  //     - (exp((8*theta*(11 + 9*tanh(10 - 4*x)))/
  //        ((-1 + gamma)*(-11198669 - 769*tanh(10 - 4*x) - 891*Power(tanh(10 - 4*x),2) + 729*Power(tanh(10 - 4*x),3))))*lambda*Q*(11 + 9*tanh(10 - 4*x))*
  //      (1 + tanh(2 - x)))/40.,(Power(1.0/cosh(2 - x),2)*(-11 + 9*tanh(10 - 4*x))*(11 + 9*tanh(10 - 4*x)) + (80*Power(1.0/cosh(2 - x),2)*tanh(2 - x))/Le +
  //     36*Power(1.0/cosh(10 - 4*x),2)*(-11 + 9*tanh(10 - 4*x))*(1 + tanh(2 - x)) +
  //     2*exp((8*theta*(11 + 9*tanh(10 - 4*x)))/
  //        ((-1 + gamma)*(-11198669 - 769*tanh(10 - 4*x) - 891*Power(tanh(10 - 4*x),2) + 729*Power(tanh(10 - 4*x),3))))*lambda*(11 + 9*tanh(10 - 4*x))*
  //      (1 + tanh(2 - x)) + 36*Power(1.0/cosh(10 - 4*x),2)*(11 + 9*tanh(10 - 4*x))*(1 + tanh(2 - x)))/80.;

    // man_sol += temp;
    ///////////////////////////////////////////////////////////////////////////////
    // Hyperbolic
    ///////////////////////////////////////////////////////////////////////////////
    temp << (81*Power(1/cosh(10 - 4*x),2)*tanh(10 - 4*x))/5.,(Power(1/cosh(10 - 4*x),2)*
      (2947 - 769*gamma - 1782*(-3 + gamma)*tanh(10 - 4*x) + 2187*(-3 + gamma)*Power(tanh(10 - 4*x),2)))/40.,
   (Power(1/cosh(10 - 4*x),2)*(11979 + 50389781*gamma - 2880*gamma*tanh(10 - 4*x) + 24057*(-1 + gamma)*Power(tanh(10 - 4*x),2) -
        13122*(-1 + gamma)*Power(tanh(10 - 4*x),3)))/40.,(Power(1/cosh(2 - x),2)*(-121 + 81*Power(tanh(10 - 4*x),2)) +
      648*Power(1/cosh(10 - 4*x),2)*tanh(10 - 4*x)*(1 + tanh(2 - x)))/80.;

    man_sol += temp;
    ///////////////////////////////////////////////////////////////////////////////
    // Viscous
    ///////////////////////////////////////////////////////////////////////////////
    temp << 0,-192*Pr*Power(1/cosh(10 - 4*x),2)*tanh(10 - 4*x),(36*Power(1/cosh(10 - 4*x),4)*
       (-100791541*gamma - 15972*Pr + 9801*(3*gamma - 4*Pr)*tanh(10 - 4*x) + 8019*(3*gamma - 4*Pr)*Power(tanh(10 - 4*x),2) +
         2187*(3*gamma - 4*Pr)*Power(tanh(10 - 4*x),3)) - 8*Power(1/cosh(10 - 4*x),2)*tanh(10 - 4*x)*
       (554287591*gamma + 175692*Pr + 18*(25188901*gamma + 15972*Pr)*tanh(10 - 4*x) + 48114*(3*gamma - 4*Pr)*Power(tanh(10 - 4*x),3) +
         19683*(3*gamma - 4*Pr)*Power(tanh(10 - 4*x),4)))/Power(11 + 9*tanh(10 - 4*x),3),(Power(1/cosh(2 - x),2)*tanh(2 - x))/Le;

    man_sol += temp;
    ///////////////////////////////////////////////////////////////////////////////
    // Source
    ///////////////////////////////////////////////////////////////////////////////
    temp << 0.,0.,-(exp((8*theta*(11 + 9*tanh(10 - 4*x)))/
         ((-1 + gamma)*(-11198669 - 769*tanh(10 - 4*x) - 891*Power(tanh(10 - 4*x),2) + 729*Power(tanh(10 - 4*x),3))))*lambda*Q*(11 + 9*tanh(10 - 4*x))*
       (1 + tanh(2 - x)))/40.,(exp((8*theta*(11 + 9*tanh(10 - 4*x)))/
        ((-1 + gamma)*(-11198669 - 769*tanh(10 - 4*x) - 891*Power(tanh(10 - 4*x),2) + 729*Power(tanh(10 - 4*x),3))))*lambda*(11 + 9*tanh(10 - 4*x))*
      (1 + tanh(2 - x)))/40.;

    man_sol += temp;

   return man_sol;
}

#endif //#ifndef IMPLICIT_MARCHING_H
