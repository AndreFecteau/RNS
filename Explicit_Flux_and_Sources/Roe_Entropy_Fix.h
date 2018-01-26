#ifndef ROE_ENTROPY_FIX_H
#define ROE_ENTROPY_FIX_H
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Header containg the definition of the class for a
///         Roe_Entropy_Fix.h
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
  /// \brief  Header containg the Roe_Entropy_Fix solver for first order terms.
///////////////////////////////////////////////////////////////////////////////
#include <math.h>
#include <algorithm>
#include "../Usefull_Headers/Variable_Vector_Isolator.h"

template <typename global_solution_vector_type>
class Roe_Entropy_Fix {
 public:
   using solution_vector_type = typename global_solution_vector_type::value_type;
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  Roe_Entropy_Fix() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  Roe_Entropy_Fix(const Roe_Entropy_Fix&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  Roe_Entropy_Fix(Roe_Entropy_Fix&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  Roe_Entropy_Fix& operator=(const Roe_Entropy_Fix&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  Roe_Entropy_Fix& operator=(Roe_Entropy_Fix&&) = default;

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
    solution_vector_type lambda;
    solution_vector_type del;
    del = del_function();
    if (del[0]/2.0 <= fabs(u_tilda - a_tilda)){
      lambda[0] = fabs(u_tilda - a_tilda);
    } else {
      lambda[0] = pow(u_tilda - a_tilda,2)/del[0] + del[0] / 4.0;
    }
    lambda[1] = fabs(u_tilda);

    if (del[2]/2.0 <= fabs(u_tilda + a_tilda)){
      lambda[2] = fabs(u_tilda + a_tilda);
    } else {
      lambda[2] = pow((u_tilda + a_tilda),2)/del[2] + del[2] / 4.0;// * fabs(u_tilda + a_tilda);
    }

    if (del[3]/2.0 <= fabs(u_tilda - a_tilda)){
      lambda[3] = fabs(u_tilda - a_tilda);
    } else {
      lambda[3] = pow(u_tilda - a_tilda,2)/del[3] + del[3] / 4.0;
    }

  Eigen::Matrix3d eigenvalues; eigenvalues << lambda[0], 0.0, 0.0, 0.0,
                                              0.0, lambda[1], 0.0, 0.0,
                                              0.0, 0.0, lambda[2], 0.0,
                                              0.0, 0.0, 0.0, lambda[3];
                                              
  Eigen::Matrix3d right_eigenvectors; right_eigenvectors << 1.0, 1.0, 1.0,
                                                            u_tilda - a_tilda, u_tilda, u_tilda + a_tilda,
                                                            H_tilda - u_tilda*a_tilda, 0.5 * pow(u_tilda,2), H_tilda + u_tilda*a_tilda;
  Eigen::Vector3d delta_U; delta_U << rho_right - rho_left, rho_right*u_right - rho_left*u_left, (p_right / (gamma- 1.0) + rho_right*pow(u_right,2)/2.0) - (p_left/(gamma- 1) + rho_left*pow(u_left,2)/2.0);
  Eigen::Matrix3d left_eigenvectors = right_eigenvectors.inverse();
  delta_F = right_eigenvectors * (eigenvalues * (left_eigenvectors * delta_U));
  F_left << rho_left * u_left, rho_left * pow(u_left, 2) + p_left, u_left * (gamma * p_left/(gamma - 1)+rho_left * pow(u_left, 2) / 2.0);
  F_right << rho_right * u_right, rho_right * pow(u_right, 2) + p_right, u_right * (gamma * p_right/(gamma - 1)+rho_right * pow(u_right, 2) / 2.0);
}

  // delta_F = 0,0,0,((rho_left - rho_right)*Power(u_rhoavg(),2)*(2*Sqrt(gamma*p_rhoavg()*rho_rhoavg()) - (-1 + gamma)*rho_rhoavg()*u_rhoavg())*Y_rhoavg())/(4.*gamma*p_rhoavg()),
  //  (-2*(-1 + gamma)*(rho_left*u_left - rho_right*u_right)*rho_rhoavg()*Y_rhoavg())/(gamma*p_rhoavg()),
  //  -(((-1 + gamma)*(rho_left*u_left - rho_right*u_right)*rho_rhoavg()*Power(u_rhoavg(),3))/(gamma*p_rhoavg())),
  //  ((rho_left*u_left - rho_right*u_right)*u_rhoavg()*(Sqrt(gamma*p_rhoavg()*rho_rhoavg()) + (-1 + gamma)*rho_rhoavg()*u_rhoavg())*Y_rhoavg())/(2.*gamma*p_rhoavg()),0,
  //  (0.5*(2.*p_left - 2.*p_right + (-1. + 1.*gamma)*rho_left*Power(u_left,2) + (1. - 1.*gamma)*rho_right*Power(u_right,2))*Power(rho_rhoavg(),2)*
  //     (-1.*gamma*p_rhoavg() + Sqrt(gamma*p_rhoavg()*rho_rhoavg())*u_rhoavg()))/Power(gamma*p_rhoavg()*rho_rhoavg(),1.5),
  //  (-0.25*(2.*p_left - 2.*p_right + (-1. + 1.*gamma)*rho_left*Power(u_left,2) + (1. - 1.*gamma)*rho_right*Power(u_right,2))*rho_rhoavg()*Power(u_rhoavg(),2)*
  //     (Sqrt(gamma*p_rhoavg()*rho_rhoavg()) - 1.*rho_rhoavg()*u_rhoavg())*(-1.*gamma*p_rhoavg() + Sqrt(gamma*p_rhoavg()*rho_rhoavg())*u_rhoavg()))/(Power(gamma*p_rhoavg()*rho_rhoavg(),1.5)*Y_rhoavg()),
  //  (0.25*(2.*p_left - 2.*p_right + (-1. + 1.*gamma)*rho_left*Power(u_left,2) + (1. - 1.*gamma)*rho_right*Power(u_right,2))*rho_rhoavg()*
  //     (-1.*gamma*p_rhoavg() + Sqrt(gamma*p_rhoavg()*rho_rhoavg())*u_rhoavg())*(-1.*gamma*p_rhoavg() + 1.*(-1. + 1.*gamma)*Sqrt(gamma*p_rhoavg()*rho_rhoavg())*u_rhoavg() +
  //       (0.5 - 0.5*gamma)*rho_rhoavg()*Power(u_rhoavg(),2)))/((-1. + gamma)*Power(gamma*p_rhoavg()*rho_rhoavg(),1.5)),
  //  (-0.25*(2.*p_left - 2.*p_right + (-1. + 1.*gamma)*rho_left*Power(u_left,2) + (1. - 1.*gamma)*rho_right*Power(u_right,2))*Power(rho_rhoavg(),2)*
  //     (-1.*gamma*p_rhoavg() + Sqrt(gamma*p_rhoavg()*rho_rhoavg())*u_rhoavg())*Y_rhoavg())/Power(gamma*p_rhoavg()*rho_rhoavg(),1.5),
  //  ((-(rho_left*Y_left) + rho_right*Y_right)*((gamma*p_rhoavg())/Sqrt(gamma*p_rhoavg()*rho_rhoavg()) + u_rhoavg()))/Y_rhoavg(),0,0,0;

   return 0.5 * (F_left() - F_right()) - 0.5*delta_F;
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
  /// \brief Function to calculate the mach number for the right state.
  double a_right() {
    return sqrt(gamma * p_right / rho_right);
  }
  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the mach number for the left state.
  double a_left() {
    return sqrt(gamma * p_left / rho_left);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the enthalpy for the right state.
  double h_right() {
    return gamma * p_right / (rho_right * (gamma - 1.0)) + u_right * u_right * 0.5;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the enthalpy for the left state.
  double h_left() {
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
  /// \brief Function to calculate the average rho state of the Y.
  double Y_rhoavg() {
    return (sqrt(rho_left) * Y_left + sqrt(rho_right) * Y_right) /
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
  /// \brief Function to calculate the.
  solution_vector_type del_function() {
    solution_vector_type del;
    del << std::max(0.0, 4.0 * ((u_right - a_right()) - (u_left - a_left())),
        << std::max(0.0, 4.0 * (u_right - u_left),
        << std::max(0.0, 4.0 * ((u_right + a_right()) - (u_left + a_left()));
    return del;
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

 p_rhoavg()ivate:
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

#endif //#ifndef ROE_ENTROPY_FIX_H
