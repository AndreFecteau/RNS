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
  solution_vector_type flux(solution_vector_type solution_vector_left,
                            solution_vector_type solution_vector_right,
                            double gamma_in) {
    gamma = gamma_in;
    left_var = Variable_Vector_Isolator<solution_vector_type>(solution_vector_left, gamma);
    right_var = Variable_Vector_Isolator<solution_vector_type>(solution_vector_right, gamma);

    if (lambda_left() > 0) {
      return F_left();
    } else if (lambda_right() < 0) {
      return F_right();
    } else {
      return  (lambda_right() * F_left() - lambda_left() * F_right()) / (lambda_right() -
              lambda_left()) + lambda_right()*lambda_left() / (lambda_right() -
              lambda_left()) * (solution_vector_right - solution_vector_left);
    }
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the enthalpy for the right state.
  double h_right() {
    return gamma * right_var.p() / (right_var.rho() * (gamma - 1.0)) + right_var.u() * right_var.u() * 0.5;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the enthalpy for the left state.
  double h_left() {
    return gamma * left_var.p() / (left_var.rho() * (gamma - 1.0)) + left_var.u() * left_var.u() * 0.5;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the average rho state.
  double rho_rhoavg() {
    return sqrt(left_var.rho() * right_var.rho());
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the average rho state of the velocity.
  double u_rhoavg() {
    return (sqrt(left_var.rho()) * left_var.u() + sqrt(right_var.rho()) * right_var.u()) /
           (sqrt(left_var.rho()) + sqrt(right_var.rho()));
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the average rho state of the enthalpy.
  double h_rhoavg() {
    return (sqrt(left_var.rho()) * h_left() + sqrt(right_var.rho()) * h_right()) /
           (sqrt(left_var.rho()) + sqrt(right_var.rho()));
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the average rho state of the temperature.
  double p_rhoavg() {
    return (h_rhoavg() - u_rhoavg()*u_rhoavg() * 0.5) * (gamma - 1.0) / gamma * rho_rhoavg();
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the left flux.
  solution_vector_type F_left() {
    return create_flux_vector(left_var.rho(), left_var.u(), left_var.p(), left_var.Y());
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the right flux.
  solution_vector_type F_right() {
    return create_flux_vector(right_var.rho(), right_var.u(), right_var.p(), right_var.Y());
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate an arbitrary flux.
  solution_vector_type create_flux_vector(double rho, double u, double p, double Y) {
    solution_vector_type partial_flux;
    partial_flux << rho * u, rho * u * u + p, u * (gamma * p / (gamma - 1.0) +
                    rho * u * u * 0.5), rho * u * Y;
    return partial_flux;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the corrected lambda_left.
  double lambda_left() {
    double lambda_1 = left_var.u() - sqrt(gamma*left_var.p()/left_var.rho());
    double lambda_2 = u_rhoavg() - sqrt(gamma*p_rhoavg()/rho_rhoavg());
    return std::min(lambda_1,lambda_2);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the corrected lambda_right.
  double lambda_right() {
    double lambda_1 = right_var.u() + sqrt(gamma*right_var.p()/right_var.rho());
    double lambda_2 = u_rhoavg() + sqrt(gamma*p_rhoavg()/rho_rhoavg());
    return std::max(lambda_1,lambda_2);
  }

 private:
  Variable_Vector_Isolator<solution_vector_type> left_var;
  Variable_Vector_Isolator<solution_vector_type> right_var;
  double gamma;
};

#endif //#ifndef HLLE_H
