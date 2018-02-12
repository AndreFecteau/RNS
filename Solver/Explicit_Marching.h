#ifndef EXPLICIT_MARCHING_H
#define EXPLICIT_MARCHING_H

#include <omp.h>
#include <math.h>
#include <limits>
#include <vector>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "../Explicit_Flux_and_Sources/HLLE.h"
#include "../Explicit_Flux_and_Sources/Centered_Difference.h"
#include "../Explicit_Flux_and_Sources/Sources.h"
#include "../Serialization/Serialization_Eigen.h"

template <typename global_solution_vector_type, typename matrix_type>
class Explicit_Marching {

  /////////////////////////////////////////////////////////////////////////
  /// \brief type for individual cell solution vector.
  using solution_vector_type = typename global_solution_vector_type::value_type;

 public:
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  Explicit_Marching() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  Explicit_Marching(const Explicit_Marching&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  Explicit_Marching(Explicit_Marching&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  Explicit_Marching& operator=(const Explicit_Marching&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  Explicit_Marching& operator=(Explicit_Marching&&) = default;

  Explicit_Marching(double Pr_in, double Le_in, double Q_in, double theta_in, double mf_in,
           double gamma_in, double number_of_cells_in,  double CFL_in, double dx_in) :
           Pr(Pr_in), Le(Le_in), Q(Q_in), theta(theta_in), mf(mf_in), gamma(gamma_in),
           CFL(CFL_in), number_of_cells(number_of_cells_in), dx(dx_in) {}

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
  Centered_Difference<global_solution_vector_type> parabolic_flux;
  Sources<global_solution_vector_type> sources;
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
  double dt =0.0;
  double residual;


  /////////////////////////////////////////////////////////////////////////
  /// \brief Calculates maximum stable timestep.
  /// \param global_solution_vector vector containing cell states from all the cells.
  void calculate_dt(const global_solution_vector_type &global_solution_vector);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  solution_vector_type minmod_limiter(const solution_vector_type &Ul, const solution_vector_type &U, const solution_vector_type &Ur);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  solution_vector_type vanalbada_limiter(const solution_vector_type &Ul, const solution_vector_type &U, const solution_vector_type &Ur);

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
  /// \brief Time marches one step of time dt.
  /// \param global_solution_vector vector containing cell states from all the cells.
  // void calculate_next_step(global_solution_vector_type &global_solution_vector);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Modifies cells corresponding to boundary conditions.
  /// \param global_solution_vector vector containing cell states from all the cells.
  void boundary_conditions(global_solution_vector_type &global_solution_vector);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  template <typename T> int sign(T val) {return (T(0) < val) - (val < T(0));}

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  // void get_residual(solution_vector_type solution_vector) {
  //   residual = solution_vector.squaredNorm()
  // }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  double squaredNorm(solution_vector_type solution_vector){
    return sqrt(solution_vector[0]*solution_vector[0] +
                solution_vector[1]*solution_vector[1] +
                solution_vector[2]*solution_vector[2] +
                solution_vector[3]*solution_vector[3]);
  }

  solution_vector_type manufactured_residual(const double lambda, const int i);

};


///////////////////////////////////////////////////////////////////////////////
// TimeMarch
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
double Explicit_Marching<global_solution_vector_type, matrix_type>::timemarch(double time_frame, global_solution_vector_type &global_solution_vector, double lambda) {
  double current_time = 0.0;
  auto phi = global_solution_vector_type(global_solution_vector.size(),
  solution_vector_type::Zero());;
  auto hyperbolic_flux_vector = global_solution_vector_type(global_solution_vector.size(),
  solution_vector_type::Zero());
  auto hyperbolic_flux_future_vector = global_solution_vector_type(global_solution_vector.size(),
  solution_vector_type::Zero());
  auto global_solution_vector_future = global_solution_vector_type(global_solution_vector.size(),
  solution_vector_type::Zero());
  double residual = 0.0;

#pragma omp parallel
  {

  while (current_time < time_frame){
  residual = 0.0;
  auto global_flux_future_vector = global_solution_vector_type(global_solution_vector.size(),
  solution_vector_type::Zero());
  auto global_flux_vector = global_solution_vector_type(global_solution_vector.size(),
  solution_vector_type::Zero());

#pragma omp single
    {
    calculate_dt(global_solution_vector);
    if(current_time + dt > time_frame) {
      dt = time_frame - current_time;
    }
    }

#pragma omp for private (hyperbolic_flux)
    for (int i = 1; i < number_of_cells-1; ++i) {
      auto var_vec_l = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i-1], gamma);
      auto var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i], gamma);
      auto var_vec_r = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i+1], gamma);
      phi[i] = vanalbada_limiter(var_vec_l.w(),var_vec.w(), var_vec_r.w());
      solution_vector_type Wl = var_vec_l.w().array() + phi[i-1].array() *dx / 2.0;
      solution_vector_type Wr = var_vec.w().array() - phi[i].array() * dx / 2.0;
      hyperbolic_flux_vector[i] = hyperbolic_flux.flux(Wl, Wr, gamma);
    }
#pragma omp for private (parabolic_flux, sources)
    for (int i = 1; i < number_of_cells-2; ++i) {
      global_flux_vector[i] += (hyperbolic_flux_vector[i] - hyperbolic_flux_vector[i+1]) / dx * dt ;
      global_flux_vector[i] += parabolic_flux.flux(global_solution_vector[i-1], global_solution_vector[i], global_solution_vector[i+1], gamma, Le, Pr, dx) * dt;
      global_flux_vector[i] += sources.flux(global_solution_vector[i], gamma, Q, lambda, theta) * dt;
      global_flux_vector[i] += manufactured_residual(lambda, i)*dt;
    }

#pragma omp single
    global_solution_vector_future = global_solution_vector;
#pragma omp for reduction (+:residual)
    for (int i = 1; i < number_of_cells-2; ++i) {
      global_solution_vector_future[i] += global_flux_vector[i];
    }

#pragma omp single
    boundary_conditions(global_solution_vector_future);

////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma omp for private (hyperbolic_flux)
    for (int i = 1; i < number_of_cells-1; ++i) {
      auto var_vec_l = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector_future[i-1], gamma);
      auto var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector_future[i], gamma);
      auto var_vec_r = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector_future[i+1], gamma);
      phi[i] = vanalbada_limiter(var_vec_l.w(),var_vec.w(), var_vec_r.w());
      solution_vector_type Ul = var_vec_l.w().array() + phi[i-1].array() *dx / 2.0;
      solution_vector_type Ur = var_vec.w().array() - phi[i].array() * dx / 2.0;
      hyperbolic_flux_future_vector[i] = hyperbolic_flux.flux(Ul, Ur, gamma);
    }

#pragma omp for private (parabolic_flux, sources)
    for (int i = 1; i < number_of_cells-2; ++i) {
      global_flux_future_vector[i] += (hyperbolic_flux_future_vector[i] - hyperbolic_flux_future_vector[i+1]) / dx * dt ;
      global_flux_future_vector[i] += parabolic_flux.flux(global_solution_vector_future[i-1], global_solution_vector_future[i], global_solution_vector_future[i+1], gamma, Le, Pr, dx) * dt;
      global_flux_future_vector[i] += sources.flux(global_solution_vector_future[i], gamma, Q, lambda, theta) * dt;
      global_flux_future_vector[i] += manufactured_residual(lambda, i)*dt;
    }

#pragma omp for reduction (+:residual)
    for (int i = 1; i < number_of_cells-2; ++i) {
      global_solution_vector[i] += (global_flux_vector[i]+global_flux_future_vector[i])*0.5;
      residual += squaredNorm(global_flux_vector[i]) * dx / dt;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma omp single
    boundary_conditions(global_solution_vector);

#pragma omp single
    current_time += dt;
  }
#pragma omp single
  std::cout << "residual: " << residual << std::endl;
  }
  return residual;
}

///////////////////////////////////////////////////////////////////////////////
// Calculate dt;
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
void Explicit_Marching<global_solution_vector_type, matrix_type>::calculate_dt(const global_solution_vector_type &global_solution_vector) {
  double dt1 = CFL * dx / lambda_eigenvalue(global_solution_vector);
  double dt2 = CFL * dx*dx / (K_value(global_solution_vector));
  dt = std::min(dt1, dt2);
}

///////////////////////////////////////////////////////////////////////////////
// MinMod;
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
typename global_solution_vector_type::value_type Explicit_Marching<global_solution_vector_type, matrix_type>::minmod_limiter(const typename global_solution_vector_type::value_type &Ul,
                                                                                                       const typename global_solution_vector_type::value_type &U,
                                                                                                       const typename global_solution_vector_type::value_type &Ur) {
  solution_vector_type phi;
  solution_vector_type a = (U - Ul) / dx;
  solution_vector_type b = (Ur - U) / dx;
  for(int i = 0; i < 4; ++i) {
  phi[i] = sign(a[i])*std::max(0.0,std::min(fabs(a[i]), sign(a[i])*b[i]));
  }
  return phi;
}

///////////////////////////////////////////////////////////////////////////////
// VanAlbata;
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
typename global_solution_vector_type::value_type Explicit_Marching<global_solution_vector_type, matrix_type>::vanalbada_limiter(const typename global_solution_vector_type::value_type &Ul,
                                                                                                       const typename global_solution_vector_type::value_type &U,
                                                                                                       const typename global_solution_vector_type::value_type &Ur) {
  solution_vector_type a = (U - Ul) / dx;
  solution_vector_type b = (Ur - U) / dx;
  double epsilon = 1.0e-6;
  solution_vector_type phi =  a.array() * b.array() * (a.array() + b.array()) / (a.array() * a.array() + b.array() * b.array() + epsilon);
  for (int i = 0; i < 3; ++i) {
    if (a[i] / b[i] <= 0.0 || b[i] == 0) {
      phi[i] = 0.0;
    }
  }
  return phi;

}

///////////////////////////////////////////////////////////////////////////////
// Boundary Conditions
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
void Explicit_Marching<global_solution_vector_type, matrix_type>::boundary_conditions(global_solution_vector_type &global_solution_vector) {
    // global_solution_vector[number_of_cells - 1] << global_solution_vector[number_of_cells - 3];
    // global_solution_vector[number_of_cells - 2] << global_solution_vector[number_of_cells - 3];
    // global_solution_vector[0] << global_solution_vector[1];
}

///////////////////////////////////////////////////////////////////////////////
// WaveSpeed
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
double Explicit_Marching<global_solution_vector_type, matrix_type>::lambda_eigenvalue(const global_solution_vector_type &global_solution_vector){
  double cst = 0.0;
  for (int i = 0; i < number_of_cells; ++i) {
    Variable_Vector_Isolator<solution_vector_type> var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i], gamma);
    if (cst < std::fabs(var_vec.u()) + sqrt(gamma * var_vec.p()/var_vec.rho())) {
      cst = std::fabs(var_vec.u()) + sqrt(gamma*var_vec.p()/var_vec.rho());
    }
  }
  return cst;
}

///////////////////////////////////////////////////////////////////////////////
// K value
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
double Explicit_Marching<global_solution_vector_type, matrix_type>::K_value(const global_solution_vector_type &global_solution_vector) {
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
// Manufactured Solution
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type, typename matrix_type>
typename global_solution_vector_type::value_type Explicit_Marching<global_solution_vector_type, matrix_type>::manufactured_residual(const double lambda, const int i) {
      solution_vector_type temp;
      solution_vector_type man_sol;
      man_sol << 0,0,0,0;
      double x = dx*(i+0.5);
      // std::cout << "x: " << x << std::endl;
      ///////////////////////////////////////////////////////////////////////////////
      // Full
      ///////////////////////////////////////////////////////////////////////////////
    //   temp << -2*(10 + cos(x))*sin(x),(-449 + 149*gamma)*sin(x) + (cos(x)*(8*Pr + 9*(-3 + gamma)*(20 + cos(x))*sin(x)))/6.,
    //  (4*Pr*cos(x)*(10 + cos(x)))/3. - (lambda*Q*Power(10 + cos(x),2))/
    //    exp((theta*(10 + cos(x)))/((-1 + gamma)*(10000 + cos(x) - Power(10 + cos(x),3)/2.))) -
    //   (gamma*(131970 + 882204*cos(x) - 19800*cos(2*x) + 1803*cos(3*x) + 70*cos(4*x) + cos(5*x)))/(8.*Power(10 + cos(x),3)) +
    //   ((10 + cos(x))*(-2 + (-1 + gamma)*(-2 + 3*Power(10 + cos(x),2)))*sin(x))/2. -
    //   (10000 + cos(x) + (-1 + gamma)*(10000 + cos(x) - Power(10 + cos(x),3)/2.))*sin(x) - (4*Pr*Power(sin(x),2))/3.,
    //  cos(x)/Le + (lambda*Power(10 + cos(x),2))/exp((theta*(10 + cos(x)))/((-1 + gamma)*(10000 + cos(x) - Power(10 + cos(x),3)/2.))) -
    //   3*Power(10 + cos(x),2)*sin(x);

      // man_sol += temp;
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

#endif //#ifndef EXPLICIT_MARCHING_H
