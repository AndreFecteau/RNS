#ifndef EXPLICIT_MARCHING_EULER_H
#define EXPLICIT_MARCHING_EULER_H

#include <omp.h>
#include <math.h>
#include <limits>
#include <vector>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "../Explicit_Flux_and_Sources/HLLE_Euler.h"
#include "../Serialization/Serialization_Eigen.h"
#include "../Usefull_Headers/Variable_Vector_Isolator_Euler.h"
#include "../Usefull_Headers/Gnuplot_Primitive_Variables0.h"

template <typename global_solution_vector_type>
class Explicit_Marching_Euler {

  /////////////////////////////////////////////////////////////////////////
  /// \brief type for individual cell solution vector.
  using solution_vector_type = typename global_solution_vector_type::value_type;

 public:
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  Explicit_Marching_Euler() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  Explicit_Marching_Euler(const Explicit_Marching_Euler&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  Explicit_Marching_Euler(Explicit_Marching_Euler&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  Explicit_Marching_Euler& operator=(const Explicit_Marching_Euler&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  Explicit_Marching_Euler& operator=(Explicit_Marching_Euler&&) = default;

  Explicit_Marching_Euler(double gamma_in, double number_of_cells_in,  double CFL_in, double dx_in) :
           gamma(gamma_in), CFL(CFL_in), number_of_cells(number_of_cells_in), dx(dx_in) {}

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
    archive(gamma, CFL, number_of_cells, dx, dt);
  }

 private:
  HLLE_Euler<global_solution_vector_type> hyperbolic_flux;
  double gamma;
  double CFL;
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
  solution_vector_type vanalbada_limiter(const solution_vector_type Ul, const solution_vector_type U, const solution_vector_type Ur);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Calculates wavespeed for timestep control.
  /// \param global_solution_vector vector containing cell states from all the cells.
  double lambda_eigenvalue(const global_solution_vector_type &global_solution_vector);

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
  double squaredNorm(solution_vector_type solution_vector){
    return sqrt(solution_vector[0]*solution_vector[0] +
                solution_vector[1]*solution_vector[1] +
                solution_vector[2]*solution_vector[2]);
  }
};


///////////////////////////////////////////////////////////////////////////////
// TimeMarch
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type>
double Explicit_Marching_Euler<global_solution_vector_type>::timemarch(double time_frame, global_solution_vector_type &global_solution_vector, double lambda) {
  double current_time = 0.0;
  auto phi = global_solution_vector_type(global_solution_vector.size(),
  solution_vector_type::Zero());
  auto phi_futire = global_solution_vector_type(global_solution_vector.size(),
  solution_vector_type::Zero());
  auto hyperbolic_flux_vector = global_solution_vector_type(global_solution_vector.size(),
  solution_vector_type::Zero());
  auto hyperbolic_flux_future_vector = global_solution_vector_type(global_solution_vector.size(),
  solution_vector_type::Zero());
  auto global_solution_vector_future = global_solution_vector_type(global_solution_vector.size(),
  solution_vector_type::Zero());
  double residual = 0.0;

// #pragma omp parallel
  {
  while (current_time < time_frame){
  residual = 0.0;
  auto global_flux_vector = global_solution_vector_type(global_solution_vector.size(),
  solution_vector_type::Zero());

// #pragma omp single
    {
    calculate_dt(global_solution_vector);
    if(current_time + dt > time_frame) {
      dt = time_frame - current_time;
    }
    }
    // for (int i = 1; i < number_of_cells-1; ++i) {
    //   auto var_vec_l = Variable_Vector_Isolator_Euler<solution_vector_type>(global_solution_vector[i-1], gamma);
    //   auto var_vec = Variable_Vector_Isolator_Euler<solution_vector_type>(global_solution_vector[i], gamma);
    //   auto var_vec_r = Variable_Vector_Isolator_Euler<solution_vector_type>(global_solution_vector[i+1], gamma);
    // }
    phi[0] << 0,0,0;
// #pragma omp for private (hyperbolic_flux)
    for (int i = 1; i < number_of_cells-1; ++i) {
      auto var_vec_l = Variable_Vector_Isolator_Euler<solution_vector_type>(global_solution_vector[i-1], gamma);
      auto var_vec = Variable_Vector_Isolator_Euler<solution_vector_type>(global_solution_vector[i], gamma);
      auto var_vec_r = Variable_Vector_Isolator_Euler<solution_vector_type>(global_solution_vector[i+1], gamma);
      phi[i] = vanalbada_limiter(var_vec_l.w(),var_vec.w(), var_vec_r.w());
      // phi[i] = vanalbada_limiter(var_vec_l.w(),var_vec.w(), var_vec_r.w());
      solution_vector_type Wl = var_vec_l.w().array() + phi[i-1].array() *dx / 2.0;
      solution_vector_type Wr = var_vec.w().array() - phi[i].array() * dx / 2.0;
      hyperbolic_flux_vector[i] = hyperbolic_flux.flux(Wl, Wr, gamma);
    }
    global_solution_vector_future = global_solution_vector;
    for (int i = 1; i < number_of_cells-2; ++i) {
      global_solution_vector_future[i] =  global_solution_vector[i] + (hyperbolic_flux_vector[i] - hyperbolic_flux_vector[i+1]) / dx * dt ;
    }
    // for (int i = 1; i < number_of_cells-1; ++i) {
    //   auto var_vec_l = Variable_Vector_Isolator_Euler<solution_vector_type>(global_solution_vector_future[i-1], gamma);
    //   auto var_vec = Variable_Vector_Isolator_Euler<solution_vector_type>(global_solution_vector_future[i], gamma);
    //   auto var_vec_r = Variable_Vector_Isolator_Euler<solution_vector_type>(global_solution_vector_future[i+1], gamma);
    // }
    for (int i = 1; i < number_of_cells-1; ++i) {
      auto var_vec_l = Variable_Vector_Isolator_Euler<solution_vector_type>(global_solution_vector_future[i-1], gamma);
      auto var_vec = Variable_Vector_Isolator_Euler<solution_vector_type>(global_solution_vector_future[i], gamma);
      auto var_vec_r = Variable_Vector_Isolator_Euler<solution_vector_type>(global_solution_vector_future[i+1], gamma);
      phi_futire[i] = vanalbada_limiter(var_vec_l.w(),var_vec.w(), var_vec_r.w());
      // phi[i] = vanalbada_limiter(var_vec_l.w(),var_vec.w(), var_vec_r.w());
      solution_vector_type Ul = var_vec_l.w().array() + phi_futire[i-1].array() *dx / 2.0;
      solution_vector_type Ur = var_vec.w().array() - phi_futire[i].array() * dx / 2.0;
      hyperbolic_flux_future_vector[i] = hyperbolic_flux.flux(Ul, Ur, gamma);
    }

// #pragma omp for
    for (int i = 1; i < number_of_cells-2; ++i) {
      global_solution_vector[i] += (hyperbolic_flux_vector[i] - hyperbolic_flux_vector[i+1] +
                                    hyperbolic_flux_future_vector[i] - hyperbolic_flux_future_vector[i+1]) / (2.0*dx) * dt ;
      // residual += squaredNorm(global_flux_vector[i]) * dx / dt;
    }

// #pragma omp single
    // boundary_conditions(global_solution_vector);

// #pragma omp single
    current_time += dt;
  }
  // #pragma omp single
  std::cout << "residual: " << residual << std::endl;
  }
  return residual;
}

///////////////////////////////////////////////////////////////////////////////
// Calculate dt;
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type>
void Explicit_Marching_Euler<global_solution_vector_type>::calculate_dt(const global_solution_vector_type &global_solution_vector) {
  dt = CFL * dx / lambda_eigenvalue(global_solution_vector);
}

///////////////////////////////////////////////////////////////////////////////
// MinMod;
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type>
typename global_solution_vector_type::value_type Explicit_Marching_Euler<global_solution_vector_type>::minmod_limiter(
                                                            const typename global_solution_vector_type::value_type &Wl,
                                                            const typename global_solution_vector_type::value_type &W,
                                                            const typename global_solution_vector_type::value_type &Wr) {

  solution_vector_type phi;
  solution_vector_type a = (W - Wl) / dx;
  solution_vector_type b = (Wr - W) / dx;
  solution_vector_type dudt = (Wr - Wl) / (2*dx);

  for(int i = 0; i < 3; ++i) {
  // phi[i] = sign(a[i])*std::max(0.0,std::min(fabs(a[i]), sign(a[i])*b[i]));
  phi[i] = sign(a[i])*std::fmax(0.0,std::fmin(fabs(a[i]), sign(a[i])*b[i]));
  }
  // for (int i = 0; i < 3; ++i) {
  //   // phi[i] = (fabs(a[i] * b[i]) + a[i] * b[i])/(a[i] + b[i]);
  //   if (a[i] + b[i] = 0.0) {
  //     phi[i] = 0.0;
  //   }
  // }
  return phi;
}

///////////////////////////////////////////////////////////////////////////////
// VanAlbata;
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type>
typename global_solution_vector_type::value_type Explicit_Marching_Euler<global_solution_vector_type>::vanalbada_limiter(
                                                            const typename global_solution_vector_type::value_type Wl,
                                                            const typename global_solution_vector_type::value_type W,
                                                            const typename global_solution_vector_type::value_type Wr) {
  solution_vector_type a = (W - Wl) / (dx);
  solution_vector_type b = (Wr - W) / (dx);
  double epsilon = 1.0e-6;
  solution_vector_type phi =  a.array() * b.array() * (a.array() + b.array()) / (a.array() * a.array() + b.array() * b.array() + epsilon);
  // solution_vector_type phi =(a.array()+b.array()) * (a.array() * b.array() / (a.array() * a.array() + b.array() * b.array() + epsilon));
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
template <typename global_solution_vector_type>
void Explicit_Marching_Euler<global_solution_vector_type>::boundary_conditions(global_solution_vector_type &global_solution_vector) {
  for(int i = 0; i < 2; ++i) {
    global_solution_vector[number_of_cells - 2 + i] << global_solution_vector[number_of_cells - 3];
  }
}

///////////////////////////////////////////////////////////////////////////////
// WaveSpeed
///////////////////////////////////////////////////////////////////////////////
template <typename global_solution_vector_type>
double Explicit_Marching_Euler<global_solution_vector_type>::lambda_eigenvalue(const global_solution_vector_type &global_solution_vector){
  double cst = 0.0;
  for (int i = 0; i < number_of_cells; ++i) {
    Variable_Vector_Isolator_Euler<solution_vector_type> var_vec = Variable_Vector_Isolator_Euler<solution_vector_type>(global_solution_vector[i], gamma);
    if (cst < std::fabs(var_vec.u()) + sqrt(gamma * var_vec.p()/var_vec.rho())) {
      cst = std::fabs(var_vec.u()) + sqrt(gamma*var_vec.p()/var_vec.rho());
    }
  }
  return cst;
}

#endif //#ifndef EXPLICIT_MARCHING_EULER_H