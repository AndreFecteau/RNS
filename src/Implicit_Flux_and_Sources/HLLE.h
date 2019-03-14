#ifndef HLLE_H
#define HLLE_H
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Header containg the definition of the class for a
///         HLLE.h
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
  /// \brief  Header containg the HLLE solver for first order terms.
///////////////////////////////////////////////////////////////////////////////
#include <math.h>
#include <algorithm>
#include "../Usefull_Headers/Variable_Vector_Isolator.h"

template <typename global_solution_vector_type>
class HLLE {
 public:
   using solution_vector_type = typename global_solution_vector_type::value_type;
   using scalar_type = typename global_solution_vector_type::value_type::value_type;
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  HLLE() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  HLLE(const HLLE&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  HLLE(HLLE&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  HLLE& operator=(const HLLE&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  HLLE& operator=(HLLE&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the flux of the given cells.
  /// \return Flux between cells.
  solution_vector_type flux(const solution_vector_type &variable_vector_left,
                            const solution_vector_type &variable_vector_right,
                            const scalar_type &gamma_in) {
    gamma = gamma_in;
    rho_right =variable_vector_right[0];
    u_right  = variable_vector_right[1];
    p_right  = variable_vector_right[2];
    Y_right  = variable_vector_right[3];
    rho_left =  variable_vector_left[0];
    u_left   =  variable_vector_left[1];
    p_left   =  variable_vector_left[2];
    Y_left   =  variable_vector_left[3];

    if (lambda_left() > 0) {
      return F_left();
    } else if (lambda_right() < 0) {
      return F_right();
    } else {
      return  (lambda_right() * F_left() - lambda_left() * F_right()) / (lambda_right() -
              lambda_left()) + lambda_right()*lambda_left() / (lambda_right() -
              lambda_left()) * (solution_vector_right() - solution_vector_left());
    }
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the enthalpy for the right state.
  solution_vector_type solution_vector_right() {
    solution_vector_type temp; temp << rho_right, rho_right*u_right, p_right / (gamma-1.0) + 0.5*rho_right*u_right*u_right, rho_right*Y_right;
    return temp;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the enthalpy for the right state.
  solution_vector_type solution_vector_left() {
    solution_vector_type temp; temp << rho_left, rho_left*u_left, p_left / (gamma-1.0) + 0.5*rho_left*u_left*u_left, rho_left*Y_left;
    return temp;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the enthalpy for the right state.
  scalar_type h_right() {
    // return (0.5 * rho_right * u_right * u_right + p_right / (gamma-1) + p_right) / rho_right;
    return gamma * p_right / (rho_right * (gamma - 1.0)) + u_right * u_right * 0.5;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the enthalpy for the left state.
  scalar_type h_left() {
    // return (0.5 * rho_left * u_left * u_left + p_left / (gamma-1) + p_left) / rho_left;
    return gamma * p_left / (rho_left * (gamma - 1.0)) + u_left * u_left * 0.5;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the average rho state.
  scalar_type rho_rhoavg() {
    return sqrt(rho_left * rho_right);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the average rho state of the velocity.
  scalar_type u_rhoavg() {
    return (sqrt(rho_left) * u_left + sqrt(rho_right) * u_right) /
           (sqrt(rho_left) + sqrt(rho_right));
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the average rho state of the enthalpy.
  scalar_type h_rhoavg() {
    return (sqrt(rho_left) * h_left() + sqrt(rho_right) * h_right()) /
           (sqrt(rho_left) + sqrt(rho_right));
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the average rho state of the temperature.
  scalar_type p_rhoavg() {
    return (h_rhoavg() - u_rhoavg()*u_rhoavg() * 0.5) * (gamma - 1.0) / gamma * rho_rhoavg();
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the left flux.
  solution_vector_type F_left() {
    return create_flux_vector(rho_left, u_left, p_left, Y_left);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the right flux.
  solution_vector_type F_right() {
    return create_flux_vector(rho_right, u_right, p_right, Y_right);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate an arbitrary flux.
  solution_vector_type create_flux_vector(scalar_type rho, scalar_type u, scalar_type p, scalar_type Y) {
    solution_vector_type partial_flux;
    partial_flux << rho * u, rho * u * u + p, u * (gamma * p / (gamma - 1.0) +
                    rho * u * u * 0.5),rho * u * Y;
    return partial_flux;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the corrected lambda_left.
  scalar_type lambda_left() {
    scalar_type lambda_1 = u_left - sqrt(gamma*p_left/rho_left);
    // scalar_type lambda_2 = u_rhoavg() - sqrt((gamma - 1) * (h_rhoavg() - 0.5 * u_rhoavg()*u_rhoavg()));
    scalar_type lambda_2 = u_rhoavg() - sqrt(gamma*p_rhoavg()/rho_rhoavg());
    return std::min(lambda_1,lambda_2);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the corrected lambda_right.
  scalar_type lambda_right() {
    scalar_type lambda_1 = u_right + sqrt(gamma*p_right/rho_right);
    // scalar_type lambda_2 = u_rhoavg() + sqrt((gamma - 1) * (h_rhoavg() - 0.5 * u_rhoavg()*u_rhoavg()));
    scalar_type lambda_2 = u_rhoavg() + sqrt(gamma*p_rhoavg()/rho_rhoavg());
    return std::max(lambda_1,lambda_2);
  }

 private:
  scalar_type rho_right;
  scalar_type u_right;
  scalar_type p_right;
  scalar_type rho_left;
  scalar_type u_left;
  scalar_type p_left;
  scalar_type Y_right;
  scalar_type Y_left;
  scalar_type gamma;
};

#endif //#ifndef HLLE_H
