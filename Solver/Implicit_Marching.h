#ifndef IMPLICIT_MARCHING_H
#define IMPLICIT_MARCHING_H

#include <math.h>
// #include "../Implicit_Flux_and_Sources/Variable_Implicit_Scheme.h"
// #include "../Implicit_Flux_and_Sources/Variable_Implicit_Scheme_4th_Order.h"
// #include "../Implicit_Flux_and_Sources/Variable_Implicit_Scheme_HLLE.h"
// #include "../Matrix_Inverse/Gaussian_Block_Triagonal_Matrix_Inverse.h"
#include "../Matrix_Inverse/Thomas_Block_Triagonal_Matrix_Inverse.h"
#include "../Usefull_Headers/Variable_Vector_Isolator.h"

template <typename grid_type, typename flow_properties_type>
class Implicit_Marching {

  using global_solution_vector_type = typename grid_type::global_solution_vector_type;
  using solution_vector_type = typename global_solution_vector_type::value_type;
  using matrix_type = typename grid_type::matrix_type;
  using scalar_type = typename grid_type::scalar_type;
  using size_type = typename grid_type::size_type;
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

  Implicit_Marching(scalar_type Theta_in, scalar_type zeta_in) : Theta(Theta_in), zeta(zeta_in) {}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Execute step in time.
  /// \param frame_time Time to stop the time marching.
  /// \param global_solution_vector vector containing cell states from all the cells.
  template <typename flux_type>
  scalar_type timemarch(flow_properties_type flow, grid_type& grid, scalar_type CFL, scalar_type frame_time);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  // scalar_type get_dx() {return dx;}

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  // void reduce_CFL() {CFL /= 2.0; std::cout << "Reduced CFL: " << CFL << std::endl;}

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  template<typename Archive>
  void serialize(Archive& archive) {
    archive(zeta, Theta);
  }

 private:
  scalar_type Theta;
  scalar_type zeta;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Calculates maximum stable timestep.
  /// \param global_solution_vector vector containing cell states from all the cells.
  scalar_type calculate_dt(const grid_type& grid,  const flow_properties_type& flow, const scalar_type& CFL);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Calculates wavespeed for timestep control.
  /// \param global_solution_vector vector containing cell states from all the cells.
  scalar_type lambda_eigenvalue(const grid_type& grid, const flow_properties_type& flow);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Calculates the variable for timestep control from second order
  ///        derivatives.
  /// \param global_solution_vector vector containing cell states from all the cells.
  scalar_type K_value(const grid_type& grid, const flow_properties_type& flow);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  scalar_type squaredNorm(const solution_vector_type solution_vector){
    return sqrt(solution_vector[0]*solution_vector[0] +
                solution_vector[1]*solution_vector[1] +
                solution_vector[2]*solution_vector[2] +
                solution_vector[3]*solution_vector[3]);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  solution_vector_type numerical_dissipation( const grid_type& grid, const size_type i, const scalar_type omega);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  solution_vector_type manufactured_residual(const scalar_type lambda, const size_type i, grid_type grid, flow_properties_type flow);

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  template <typename T>
  scalar_type Power(const T num, const size_type expo) {return pow(num, expo);}
};


///////////////////////////////////////////////////////////////////////////////
// TimeMarch
///////////////////////////////////////////////////////////////////////////////
template <typename grid_type, typename flow_properties_type>
template <typename flux_type>
typename grid_type::scalar_type Implicit_Marching<grid_type, flow_properties_type>::
timemarch(flow_properties_type flow,
          grid_type& grid, scalar_type CFL,
          scalar_type frame_time) {

  scalar_type current_time = 0.0;
  scalar_type dt = 0.0;
  scalar_type residual;
  std::vector<matrix_type>    mid(grid.global_solution_vector.size()-2);
  std::vector<matrix_type>    top(grid.global_solution_vector.size()-2);
  std::vector<matrix_type>    bot(grid.global_solution_vector.size()-2);
  global_solution_vector_type rhs(grid.global_solution_vector.size()-2);
  global_solution_vector_type delta_global_solution_vector = global_solution_vector_type(grid.global_solution_vector.size()-2,
  solution_vector_type::Zero());

#pragma omp parallel
  {
  while (current_time < frame_time){

#pragma omp single
{
  dt = calculate_dt(grid, flow, CFL);
  if(current_time + dt > frame_time) {
    dt = frame_time - current_time;
  }
}

#pragma omp for
  for(size_type i = 1; i < grid.number_of_cells - 1; ++i) {
    auto matrix_entries = flux_type(grid.global_solution_vector[std::max(i-2,static_cast<size_type>(0))],
                                    grid.global_solution_vector[i-1],
                                    grid.global_solution_vector[i],
                                    grid.global_solution_vector[i+1],
                                    grid.global_solution_vector[std::min(i+2,grid.number_of_cells-1)],
                                    delta_global_solution_vector[i-1],
                                    flow.gamma, flow.Pr, flow.Le, flow.Q(), flow.lambda,
                                    flow.theta(), grid.dx(), dt, zeta, Theta, flow.mf, flow.T_ignition());
    mid[i-1] = matrix_entries.mid_matrix();
    bot[i-1] = matrix_entries.bot_matrix();
    top[i-1] = matrix_entries.top_matrix();
    rhs[i-1] = matrix_entries.rhs_matrix();
#if defined(MANUFACTURED)
    rhs[i-1] += manufactured_residual(lambda, i)*dt/(1+zeta);
#endif
    rhs[i-1] += numerical_dissipation(grid, i, 0.9);
    // rhs[i-1] += numerical_dissipation(grid.global_solution_vector, i, 0.5);
    // rhs[i-1] += numerical_dissipation(grid.global_solution_vector, i, 0.01);
  }

// Implicit Boundary Conditions
#pragma omp single
{
#if !defined(MANUFACTURED)
  // mid[grid.global_solution_vector.size()-3] += top[grid.global_solution_vector.size()-3];
  // mid[0] += bot[0];
#endif
}

#pragma omp single
  delta_global_solution_vector = block_triagonal_matrix_inverse<matrix_type, solution_vector_type>(mid, top, bot, rhs);

#pragma omp for
  for (size_type i = 1; i < grid.number_of_cells-1; ++i) {
    if(current_time == 0.0) {
      grid.global_solution_vector[i] += delta_global_solution_vector[i-1]*(1.0+zeta);
    } else {
      grid.global_solution_vector[i] += delta_global_solution_vector[i-1];
    }
  }

// Explicit Boundary Conditions
#pragma omp single
  {
#if !defined(MANUFACTURED)
    // grid.global_solution_vector[grid.global_solution_vector.size()-1] = grid.global_solution_vector[grid.global_solution_vector.size()-2];
    // grid.global_solution_vector[0] = grid.global_solution_vector[1];
#endif
  }

#pragma omp single
  current_time += dt;
// #pragma omp single
//   std::cout << current_time << std::endl;
  }
  residual = 0.0;
#pragma omp for reduction(+: residual)
  for (size_type i = 1; i < grid.number_of_cells-1; ++i) {
    residual += delta_global_solution_vector[i-1].squaredNorm() * grid.dx() / dt;
  }
  }
  std::cout << "residual: " << residual << " : " << CFL << std::endl;
  return residual;
}

///////////////////////////////////////////////////////////////////////////////
// Calculate dt;
///////////////////////////////////////////////////////////////////////////////
template <typename grid_type, typename flow_properties_type>
typename grid_type::scalar_type Implicit_Marching<grid_type, flow_properties_type>::
calculate_dt(const grid_type& grid, const flow_properties_type& flow, const scalar_type& CFL) {
  scalar_type dt1 = CFL * grid.dx() / lambda_eigenvalue(grid, flow);
  scalar_type dt2 = CFL * grid.dx()*grid.dx() / (K_value(grid, flow));

  if(isnan(grid.global_solution_vector[1][2])){
    dt1 = 1e4;
    dt2 = 1e4;
  }

  // std::cout << dt1 << " : " << dt2 << std::endl;
  return std::min(fabs(dt1), fabs(dt2));
}

///////////////////////////////////////////////////////////////////////////////
// WaveSpeed
///////////////////////////////////////////////////////////////////////////////
template <typename grid_type, typename flow_properties_type>
typename grid_type::scalar_type Implicit_Marching<grid_type, flow_properties_type>::
lambda_eigenvalue(const grid_type& grid, const flow_properties_type& flow){
  scalar_type wavespeed = 0.0;
  for (size_type i = 0; i < grid.number_of_cells; ++i) {
    Variable_Vector_Isolator<grid_type> var_vec = Variable_Vector_Isolator<grid_type>(grid.global_solution_vector[i], flow.gamma);
    if (wavespeed < std::fabs(var_vec.u()) + sqrt(flow.gamma * var_vec.p()/var_vec.rho())) {
      wavespeed = std::fabs(var_vec.u()) + sqrt(flow.gamma*var_vec.p()/var_vec.rho());
    }
  }
  return wavespeed;
}

///////////////////////////////////////////////////////////////////////////////
// K value
///////////////////////////////////////////////////////////////////////////////
template <typename grid_type, typename flow_properties_type>
typename grid_type::scalar_type Implicit_Marching<grid_type, flow_properties_type>::
K_value(const grid_type& grid, const flow_properties_type& flow) {
  scalar_type min_rho = std::numeric_limits<scalar_type>::max();
  for (size_type i = 0; i < grid.number_of_cells; ++i) {
    auto var_vec = Variable_Vector_Isolator<grid_type>(grid.global_solution_vector[i], flow.gamma);
    if (min_rho > var_vec.rho()) {
      min_rho = var_vec.rho();
    }
  }
  return 4.0 * std::max(std::max(flow.Pr, flow.Le), flow.gamma / (flow.gamma - 1.0)) / min_rho;
}

///////////////////////////////////////////////////////////////////////////////
// Add Numerical Dissipation
///////////////////////////////////////////////////////////////////////////////
template <typename grid_type, typename flow_properties_type>
typename grid_type::global_solution_vector_type::value_type Implicit_Marching<grid_type, flow_properties_type>::
numerical_dissipation(const grid_type &grid, const size_type i, const scalar_type omega) {
            return -omega/(1.0+zeta)/8.0*(grid.global_solution_vector[std::min(i+2,grid.number_of_cells-1)] -
                                        4.0*grid.global_solution_vector[std::min(i+1,grid.number_of_cells-1)] +
                                        6.0*grid.global_solution_vector[i] -
                                        4.0*grid.global_solution_vector[std::max(i-1,static_cast<size_type>(0))] +
                                        grid.global_solution_vector[std::max(i-2,static_cast<size_type>(0))]);
}

///////////////////////////////////////////////////////////////////////////////
// Manufactured Solution
///////////////////////////////////////////////////////////////////////////////
template <typename grid_type, typename flow_properties_type>
typename grid_type::global_solution_vector_type::value_type Implicit_Marching<grid_type, flow_properties_type>::manufactured_residual(const scalar_type lambda, const size_type i, grid_type grid, flow_properties_type flow) {
    solution_vector_type temp;
    solution_vector_type man_sol;
    scalar_type Pr = flow.Pr;
    scalar_type Le = flow.Le;
    scalar_type Q = flow.Q();
    scalar_type theta = flow.theta();
    scalar_type gamma = flow.gamma;

    man_sol << 0,0,0,0;
    scalar_type x = grid.dx()*(i+0.5);
    ///////////////////////////////////////////////////////////////////////////////
    // Hyperbolic
    ///////////////////////////////////////////////////////////////////////////////
#if defined(HYPERBOLIC)
    temp << (81*Power(1/cosh(10 - 4*x),2)*tanh(10 - 4*x))/5.,(Power(1/cosh(10 - 4*x),2)*
      (2947 - 769*gamma - 1782*(-3 + gamma)*tanh(10 - 4*x) + 2187*(-3 + gamma)*Power(tanh(10 - 4*x),2)))/40.,
   (Power(1/cosh(10 - 4*x),2)*(11979 + 50389781*gamma - 2880*gamma*tanh(10 - 4*x) + 24057*(-1 + gamma)*Power(tanh(10 - 4*x),2) -
        13122*(-1 + gamma)*Power(tanh(10 - 4*x),3)))/40.,(Power(1/cosh(2 - x),2)*(-121 + 81*Power(tanh(10 - 4*x),2)) +
      648*Power(1/cosh(10 - 4*x),2)*tanh(10 - 4*x)*(1 + tanh(2 - x)))/80.;

    man_sol += temp;
#endif
    ///////////////////////////////////////////////////////////////////////////////
    // Viscous
    ///////////////////////////////////////////////////////////////////////////////
#if defined(VISCOUS)
    temp << 0,-192*Pr*Power(1/cosh(10 - 4*x),2)*tanh(10 - 4*x),(36*Power(1/cosh(10 - 4*x),4)*
       (-100791541*gamma - 15972*Pr + 9801*(3*gamma - 4*Pr)*tanh(10 - 4*x) + 8019*(3*gamma - 4*Pr)*Power(tanh(10 - 4*x),2) +
         2187*(3*gamma - 4*Pr)*Power(tanh(10 - 4*x),3)) - 8*Power(1/cosh(10 - 4*x),2)*tanh(10 - 4*x)*
       (554287591*gamma + 175692*Pr + 18*(25188901*gamma + 15972*Pr)*tanh(10 - 4*x) + 48114*(3*gamma - 4*Pr)*Power(tanh(10 - 4*x),3) +
         19683*(3*gamma - 4*Pr)*Power(tanh(10 - 4*x),4)))/Power(11 + 9*tanh(10 - 4*x),3),(Power(1/cosh(2 - x),2)*tanh(2 - x))/Le;

    man_sol += temp;
#endif
    ///////////////////////////////////////////////////////////////////////////////
    // Source
    ///////////////////////////////////////////////////////////////////////////////
#if defined(SOURCE)
    temp << 0.,0.,-(exp((8*theta()*(11 + 9*tanh(10 - 4*x)))/
         ((-1 + gamma)*(-11198669 - 769*tanh(10 - 4*x) - 891*Power(tanh(10 - 4*x),2) + 729*Power(tanh(10 - 4*x),3))))*lambda*Q()*(11 + 9*tanh(10 - 4*x))*
       (1 + tanh(2 - x)))/40.,(exp((8*theta()*(11 + 9*tanh(10 - 4*x)))/
        ((-1 + gamma)*(-11198669 - 769*tanh(10 - 4*x) - 891*Power(tanh(10 - 4*x),2) + 729*Power(tanh(10 - 4*x),3))))*lambda*(11 + 9*tanh(10 - 4*x))*
      (1 + tanh(2 - x)))/40.;

    man_sol += temp;
#endif

   return man_sol;
}

#endif //#ifndef IMPLICIT_MARCHING_H
