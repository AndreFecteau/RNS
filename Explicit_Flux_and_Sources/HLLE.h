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
  solution_vector_type flux(solution_vector_type variable_vector_left,
                            solution_vector_type variable_vector_right,
                            double gamma_in) {
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
    solution_vector_type temp; temp << rho_right, rho_right*u_right, p_right / (gamma-1) + 0.5*rho_right*u_right*u_right, rho_right*Y_right;
    return temp;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the enthalpy for the right state.
  solution_vector_type solution_vector_left() {
    solution_vector_type temp; temp << rho_left, rho_left*u_left, p_left / (gamma-1) + 0.5*rho_left*u_left*u_left, rho_left*Y_left;
    return temp;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the enthalpy for the right state.
  double h_right() {
    // return (0.5 * rho_right * u_right * u_right + p_right / (gamma-1) + p_right) / rho_right;
    return gamma * p_right / (rho_right * (gamma - 1.0)) + u_right * u_right * 0.5;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the enthalpy for the left state.
  double h_left() {
    // return (0.5 * rho_left * u_left * u_left + p_left / (gamma-1) + p_left) / rho_left;
    return gamma * p_left / (rho_left * (gamma - 1.0)) + u_left * u_left * 0.5;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the average rho state.
  double rho_rhoavg() {
    return sqrt(rho_left * rho_right);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the average rho state of the velocity.
  double u_rhoavg() {
    return (sqrt(rho_left) * u_left + sqrt(rho_right) * u_right) /
           (sqrt(rho_left) + sqrt(rho_right));
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the average rho state of the enthalpy.
  double h_rhoavg() {
    return (sqrt(rho_left) * h_left() + sqrt(rho_right) * h_right()) /
           (sqrt(rho_left) + sqrt(rho_right));
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the average rho state of the temperature.
  double p_rhoavg() {
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
  solution_vector_type create_flux_vector(double rho, double u, double p, double Y) {
    solution_vector_type partial_flux;
    partial_flux << rho * u, rho * u * u + p, u * (gamma * p / (gamma - 1.0) +
                    rho * u * u * 0.5),rho * u * Y;
    return partial_flux;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the corrected lambda_left.
  double lambda_left() {
    double lambda_1 = u_left - sqrt(gamma*p_left/rho_left);
    // double lambda_2 = u_rhoavg() - sqrt((gamma - 1) * (h_rhoavg() - 0.5 * u_rhoavg()*u_rhoavg()));
    double lambda_2 = u_rhoavg() - sqrt(gamma*p_rhoavg()/rho_rhoavg());
    return std::min(lambda_1,lambda_2);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the corrected lambda_right.
  double lambda_right() {
    double lambda_1 = u_right + sqrt(gamma*p_right/rho_right);
    // double lambda_2 = u_rhoavg() + sqrt((gamma - 1) * (h_rhoavg() - 0.5 * u_rhoavg()*u_rhoavg()));
    double lambda_2 = u_rhoavg() + sqrt(gamma*p_rhoavg()/rho_rhoavg());
    return std::max(lambda_1,lambda_2);
  }

 private:
  double rho_right;
  double u_right;
  double p_right;
  double rho_left;
  double u_left;
  double p_left;
  double Y_right;
  double Y_left;
  double gamma;
};

#endif //#ifndef HLLE_H