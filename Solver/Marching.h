#ifndef MARCHING_H
#define MARCHING_H

#include <omp.h>
#include <math.h>
#include <limits>
#include <vector>
#include "../Flux_and_Sources/HLLE.h"
#include "../Flux_and_Sources/Centered_Difference.h"
#include "../Flux_and_Sources/Sources.h"
#include "../Serialization/Serialization_Eigen.h"

template <typename global_solution_vector_type>
class Marching {

  /////////////////////////////////////////////////////////////////////////
  /// \brief type for individual cell solution vector.
  using solution_vector_type = typename global_solution_vector_type::value_type;

 public:
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  Marching() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  Marching(const Marching&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  Marching(Marching&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  Marching& operator=(const Marching&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  Marching& operator=(Marching&&) = default;

  Marching(double Pr_in, double Le_in, double Q_in, double theta_in, double mf_in,
           double gamma_in, double number_of_cells_in,  double CFL_in, double dx_in) :
           Pr(Pr_in), Le(Le_in), Q(Q_in), theta(theta_in), mf(mf_in), gamma(gamma_in),
           CFL(CFL_in), number_of_cells(number_of_cells_in), dx(dx_in) {}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Execute step in time.
  /// \param time_frame Time to stop the time marching.
  /// \param global_solution_vector vector containing cell states from all the cells.
  void timemarch(double time_frame, global_solution_vector_type &global_solution_vector, double lambda);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  double get_dx() {
    return dx;
  }
  template<typename Archive>
  void serialize(Archive& archive) {
    archive(Pr, Le, Q, theta, mf, gamma, CFL, number_of_cells, dx, dt);
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
  int number_of_cells;
  double dx;
  double dt;

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
};


///////////////////////////////////////////////////////////////////////////////
// TimeMarch
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type>
void Marching<global_solution_vector_type>::timemarch(double time_frame, global_solution_vector_type &global_solution_vector, double lambda) {
  double current_time = 0.0;
  auto phi = global_solution_vector_type(global_solution_vector.size(),
  solution_vector_type::Zero());;
  auto hyperbolic_flux_vector = global_solution_vector_type(global_solution_vector.size(),
  solution_vector_type::Zero());
  double residual = 0.0;
#pragma omp parallel num_threads (12)
  {
  #pragma omp single
  std::cout << "Number of threads being used: " << omp_get_num_threads() << std::endl;

  while (current_time < time_frame){
  residual = 0.0;
  auto global_flux_vector = global_solution_vector_type(global_solution_vector.size(),
  solution_vector_type::Zero());

#pragma omp single
    {
    calculate_dt(global_solution_vector);
    if(current_time + dt > time_frame) {
      dt = time_frame - current_time;
    }
    }
#pragma omp for
    for (int i = 1; i < number_of_cells-1; ++i) {
      phi[i] = vanalbada_limiter(global_solution_vector[i-1], global_solution_vector[i], global_solution_vector[i+1]);
    }

#pragma omp for private (hyperbolic_flux)
    for (int i = 1; i < number_of_cells; ++i) {
      solution_vector_type Ul;
      solution_vector_type Ur;
      Ul = global_solution_vector[i-1] + phi[i-1] * dx * 0.5;
      Ur = global_solution_vector[i] - phi[i] * dx * 0.5;
      hyperbolic_flux_vector[i] = hyperbolic_flux.flux(Ul, Ur, gamma);
    }

#pragma omp for private (parabolic_flux, sources)
    for (int i = 1; i < number_of_cells-1; ++i) {
      global_flux_vector[i] += (hyperbolic_flux_vector[i] - hyperbolic_flux_vector[i+1]) / dx * dt ;
      global_flux_vector[i] += parabolic_flux.flux(global_solution_vector[i-1], global_solution_vector[i], global_solution_vector[i+1], gamma, Le, Pr, dx) * dt;
      global_flux_vector[i] += sources.flux(global_solution_vector[i], gamma, Q, lambda, theta) * dt;
      global_flux_vector[i][4] /= dt;
    }

#pragma omp for reduction (+:residual)
    for (int i = 1; i < number_of_cells-1; ++i) {

      global_solution_vector[i][4] = 0.0;
      global_solution_vector[i] += global_flux_vector[i];
      residual += squaredNorm(global_flux_vector[i]) * dx / dt;
    }

#pragma omp single
    boundary_conditions(global_solution_vector);

if defined(FIX){
#pragma omp for
    for (int i = 1; i < number_of_cells - 2; ++i) {
      Variable_Vector_Isolator<solution_vector_type> var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i], gamma);
      double rho_local = var_vec.rho();
      double u_local = var_vec.u();
      double T_local = var_vec.T();
      double Y_local = var_vec.Y();
      rho_local = std::min(rho_local, 1.0);
      u_local = std::max(u_local,1.0);
      u_local = std::min(u_local, 10.0);
      T_local = std::max(T_local, 1.0 / (gamma*mf*mf));
      global_solution_vector[i] <<  rho_local,
                                    rho_local * u_local,
                                    rho_local*T_local/(gamma - 1.0) + rho_local * u_local * u_local * 0.5,
                                    rho_local*Y_local,
                                    global_solution_vector[i][4];
    }
}
#pragma omp single
    current_time += dt;
  }
  #pragma omp single
  std::cout << "residual: " << residual << std::endl;
  }
}

///////////////////////////////////////////////////////////////////////////////
// Calculate dt;
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type>
void Marching<global_solution_vector_type>::calculate_dt(const global_solution_vector_type &global_solution_vector) {
  double dt1 = CFL * dx / lambda_eigenvalue(global_solution_vector);
  double dt2 = CFL * dx*dx / (K_value(global_solution_vector));
  // std::cout << "Hyperbolic: " << dt1 << " Viscous: " << dt2 << std::endl;
  dt = std::min(dt1, dt2);
}

///////////////////////////////////////////////////////////////////////////////
// MinMod;
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type>
typename global_solution_vector_type::value_type Marching<global_solution_vector_type>::minmod_limiter(const typename global_solution_vector_type::value_type &Ul,
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
template <typename global_solution_vector_type>
typename global_solution_vector_type::value_type Marching<global_solution_vector_type>::vanalbada_limiter(const typename global_solution_vector_type::value_type &Ul,
                                                                                                       const typename global_solution_vector_type::value_type &U,
                                                                                                       const typename global_solution_vector_type::value_type &Ur) {
  solution_vector_type phi;
  solution_vector_type a = (U - Ul) / dx;
  solution_vector_type b = (Ur - U) / dx;
  for(int i = 0; i < 4; ++i) {
  phi[i] = (a[i] * b[i]) / (a[i]*a[i] + b[i]*b[i] + 0.0001) * (a[i]+b[i]);
  }
  return phi;
}
///////////////////////////////////////////////////////////////////////////////
// Calculate Next Step
///////////////////////////////////////////////////////////////////////////////
// template <typename global_solution_vector_type>
// void Marching<global_solution_vector_type>::calculate_next_step(global_solution_vector_type &global_solution_vector) {
//   auto flux_vector = global_solution_vector_type(global_solution_vector.size(),
//   solution_vector_type::Zero());
//   solution_vector_type flux_left;// = hyperbolic_flux.flux(global_solution_vector[0], global_solution_vector[1], gamma);
//   solution_vector_type flux_right;
//
//   // double t1 = omp_get_wtime();
//   #pragma omp parallel for private (flux_left, flux_right, hyperbolic_flux, parabolic_flux, sources)
//   for (int i = 1; i < number_of_cells -2; ++i) {
//     solution_vector_type phi = minmod_limiter(global_solution_vector[i], global_solution_vector[i+1], global_solution_vector[i+2]);
//     solution_vector_type Ul;
//     solution_vector_type Ur;
//     for(int j = 0; j < 4; ++j) {
//     Ul[j] = phi[j] * (global_solution_vector[i][j] - global_solution_vector[i-1][j]) * 0.5 + global_solution_vector[i][j];
//     Ur[j] = phi[j] * (global_solution_vector[i+1][j] - global_solution_vector[i+2][j]) * 0.5 + global_solution_vector[i+1][j];
//     }
//
//     flux_right = hyperbolic_flux.flux(Ul, Ur, gamma);
//
//     flux_vector[i+1] += flux_right / dx;
//     flux_vector[i] -= (flux_right)/dx;
//     flux_vector[i] += parabolic_flux.flux(global_solution_vector[i-1], global_solution_vector[i], global_solution_vector[i+1], gamma, Le, Pr, dx);
//     flux_vector[i] += sources.flux(global_solution_vector[i], gamma, Q, Lambda, theta) * dx;
//
//
//   }
//   // #pragma omp parallel for //reduction (+:global_solution_vector)
//   for (int i = 2; i < number_of_cells -2; ++i) {
//     global_solution_vector[i] += flux_vector[i] * dt;
//     global_solution_vector[i][4] = flux_vector[i][4];
//   }
//   // std::cout << omp_get_wtime() - t1 << std::endl;
//   // t_total += omp_get_wtime() - t1;
//
//   boundary_conditions(global_solution_vector);
// }

///////////////////////////////////////////////////////////////////////////////
// Boundary Conditions
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type>
void Marching<global_solution_vector_type>::boundary_conditions(global_solution_vector_type &global_solution_vector) {
  for(int i = 0; i < 2; ++i) {
    global_solution_vector[number_of_cells - 2 + i] << global_solution_vector[number_of_cells - 3];
  }
}

///////////////////////////////////////////////////////////////////////////////
// WaveSpeed
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type>
double Marching<global_solution_vector_type>::lambda_eigenvalue(const global_solution_vector_type &global_solution_vector){
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
template <typename global_solution_vector_type>
double Marching<global_solution_vector_type>::K_value(const global_solution_vector_type &global_solution_vector) {
  double min_rho = std::numeric_limits<double>::max();
  for (int i = 0; i < number_of_cells; ++i) {
    Variable_Vector_Isolator<solution_vector_type> var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i], gamma);
    if (min_rho > var_vec.rho()) {
      min_rho = var_vec.rho();
    }
  }
  return 4.0 * std::max(std::max(Pr,Le), gamma / (gamma - 1.0)) / min_rho;
}


#endif //#ifndef MARCHING_H
