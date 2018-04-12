#ifndef RK4_LOW_MACH_SOLVER_H
#define RK4_LOW_MACH_SOLVER_H
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Header containg the definition of the class for a
///         RK4_Low_Mach_Solver
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Gnuplot_RK4.h"

class RK4_Low_Mach_Solver{

///////////////////////////////////////////////////////////////////////////////
  /// \brief  Header containg a solution to the time invarient low_mach
  ///         Navier-Stokes Equations using Runge-Kutta 4. Based on "Sharpe, G.J.,
  ///         2003, Linear stability of planar premixed flames: Reactive
  ///         Navier-Stokes equations with finite activation energy and arbitrary
  ///         Lewis number. Combustion Theory and Modelling."
///////////////////////////////////////////////////////////////////////////////

public:
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  RK4_Low_Mach_Solver() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  RK4_Low_Mach_Solver(const RK4_Low_Mach_Solver&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  RK4_Low_Mach_Solver(RK4_Low_Mach_Solver&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  RK4_Low_Mach_Solver& operator=(const RK4_Low_Mach_Solver&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  RK4_Low_Mach_Solver& operator=(RK4_Low_Mach_Solver&&) = default;

  RK4_Low_Mach_Solver(double le_in, double Q_in, double theta_in, double T_ignition_in) :
    number_of_nodes(100000), Le(le_in), Q(Q_in), theta(theta_in), T_ignition(T_ignition_in) {

    if (Q < 4 || Q > 9) {
      std::cout << "Make sure that Q :" << Q << " is the correct value." << std::endl;
    }
    if (theta > 200) {
      std::cout << "Make sure that theta :" << theta << " is the correct value." << std::endl;
    }
    if (Le < 0.3 || Le > 1.8) {
      std::cout << "Make sure that Le: " << Le << " is the correct value" << std::endl;
    }
    if (T_ignition < 1 || T_ignition > 2) {
      std::cout << "Make sure that T_ignition: " << T_ignition << " is the correct value" << std::endl;
    }
    resize_arrays();
    T_final = 1.0 + Q;
    delta_T = (T_initial-T_final)/(number_of_nodes-1);
    make_low_mach_reactive_solution();
    re_position_positive();
    // export_data_to_dat();
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Returns Y at a location.
  double get_Y(double location) {return interpolate(Y_x, location);}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Returns Z at a location.
  double get_Z(double location) {return interpolate(Z_x, location);}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Returns T at a location.
  double get_T(double location) {return interpolate(T_x, location);}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Returns U at a location.
  double get_U(double location) {return get_T(location);}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Returns rho at a location.
  double get_rho(double location) {return 1 / get_T(location);}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Returns lambda at a location.
  double get_lambda() {return lambda;}

  /////////////////////////////////////////////////////////////////////////
  /// \brief returns the length of the domaine.
  double length();

  /////////////////////////////////////////////////////////////////////////
  /// \brief Positions the solution starting at x = 0;
  void re_position_positive();

  /////////////////////////////////////////////////////////////////////////
  /// \brief Creates a .dat file containing solution to
  ///        the low_mach solution.
  void export_data_to_dat();
private:

  /////////////////////////////////////////////////////////////////////////
  /// \brief Gets the index from a location.
  int get_i(double location);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Gets a value between nodes.
  double interpolate(std::vector<double> T, double location);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Heaviside function for temperature.
  double heavyside_T(double x);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Resize all variable arrays.
  void resize_arrays();

  /////////////////////////////////////////////////////////////////////////
  /// \brief Build internal structure of the flame in temperature and
  ///        distance space.
  /// \param Le     Lewis Number.
  /// \param Q      Heat of reaction.
  /// \param theta  Activation Energy.
  /// \param T_ignition_in    Actication Temperature
  void make_low_mach_reactive_solution();

  /////////////////////////////////////////////////////////////////////////
  /// \brief Calculates the value of B, see documentation for info on B;
  /// \return double B;
  double B() {return -lambda * exp(-theta/(1+Q)) / (1+Q);}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Calculates the value of h, see documentation for info on h;
  /// \return double h;
  double h() {return (Le - std::sqrt(Le*Le - 4*Le*B()))/2;}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Sets up boundary conditions.
  void setup_boundary_conditions();

  /////////////////////////////////////////////////////////////////////////
  /// \brief One itteration of runge_kutta_4.
  //TODO:(Andre) Seperate Distance space for only the final solution.
  void runge_kutta_4(int i);

  /////////////////////////////////////////////////////////////////////////
  /// \brief  Modifies lambda, lambda_min and lambda_max using bisection to
  ///         find the actual lambda.
  void bisection_lambda();

  double lambda_min = 0.0;
  double lambda = 3.0e6;
  double lambda_max = 3.0e12;
  double accuracy = 1.0e-19;
  double wb = 1.0e-4;
  double T_initial = 1.0;
  int number_of_nodes;
  std::vector<double> Z;
  std::vector<double> T;
  std::vector<double> Y;
  std::vector<double> x;
  std::vector<double> Y_x;
  std::vector<double> T_x;
  std::vector<double> Z_x;
  double Le;
  double Q;
  double theta;
  double T_final;
  double T_ignition;
  double delta_T;
  double delta_x;
  double dydT;
  double dzdT;
  double dxdT;
};

///////////////////////////////////////////////////////////////////////////////
//Gets a value between nodes.
///////////////////////////////////////////////////////////////////////////////
double RK4_Low_Mach_Solver::interpolate(std::vector<double> T, double location) {
  int i = get_i(location);
  return (T[i+1] - T[i]) * (location - x[i]) / (x[i+1] - x[i]) + T[i];
}

///////////////////////////////////////////////////////////////////////////////
//Gets the index from a location.
///////////////////////////////////////////////////////////////////////////////
int RK4_Low_Mach_Solver::get_i(double location) {
  int i = 0;
  while(((x[i] - location) * (x[i+1] - location)) > 0){
    ++i;
  }
  return i;
}

///////////////////////////////////////////////////////////////////////////////
// Creates a .dat file containing solution to the low_mach solution.
///////////////////////////////////////////////////////////////////////////////
void RK4_Low_Mach_Solver::export_data_to_dat() {
  std::vector<std::vector<double>> data(8);
  data = {T, x, Y, Z, Y_x, T_x, Z_x};
  std::string variables = "# T x Y Z Y_x T_x Z_x";
  gnuplot_RK4<double>("RK4", data, variables);
}

///////////////////////////////////////////////////////////////////////////////
// Positions the solution starting at x = 0;
///////////////////////////////////////////////////////////////////////////////
void RK4_Low_Mach_Solver::re_position_positive() {
  double change_in_x = *std::min_element(std::begin(x), std::end(x));
  for (size_t i = 0; i < x.size(); ++i) {
    x[i] -= change_in_x;
  }
}


///////////////////////////////////////////////////////////////////////////////
// return length of domaine.
///////////////////////////////////////////////////////////////////////////////
double RK4_Low_Mach_Solver::length() {
  auto it_max = std::max_element(std::begin(x), std::end(x));
  auto it_min = std::min_element(std::begin(x), std::end(x));
  return *it_max - *it_min;
}

///////////////////////////////////////////////////////////////////////////////
// Resize the arrays for the number of nodes in the solution.
///////////////////////////////////////////////////////////////////////////////
void RK4_Low_Mach_Solver::resize_arrays() {
  Y.resize(number_of_nodes-1);
  Z.resize(number_of_nodes-1);
  T.resize(number_of_nodes-1);
  x.resize(number_of_nodes-1);
  Y_x.resize(number_of_nodes-1);
  T_x.resize(number_of_nodes-1);
  Z_x.resize(number_of_nodes-1);
}

///////////////////////////////////////////////////////////////////////////////
// Runs all required functions to solve the pde's
///////////////////////////////////////////////////////////////////////////////
void RK4_Low_Mach_Solver::make_low_mach_reactive_solution() {
  while (Z[number_of_nodes-2] < 1.0-accuracy || Z[number_of_nodes-2] > 1.0 + accuracy){
    setup_boundary_conditions();
    for(int i = 0; i < number_of_nodes-2; ++i) {
      runge_kutta_4(i);
    }
    bisection_lambda();
  }
}

///////////////////////////////////////////////////////////////////////////////
// setup all boundary conditions for the shooting method
///////////////////////////////////////////////////////////////////////////////
void RK4_Low_Mach_Solver::setup_boundary_conditions() {
  T[0] = 1.0 + Q;
  Y[0] = h()*(1-h())/(Q*B()) * wb;
  Z[0] = (1 - h()) / Q * wb;
  x[0] = 0;
  Y_x[0] = Y[0];
  T_x[0] = T[0];
  Z_x[0] = Z[0];
}

///////////////////////////////////////////////////////////////////////////////
// Runge kutta 4
///////////////////////////////////////////////////////////////////////////////
void RK4_Low_Mach_Solver::runge_kutta_4(int i) {
  dydT = (Le*(Y[i]-Z[i]) / (T[i] - 1 + Q*(Z[i] - 1)));
  dzdT = -(lambda * Y[i] * exp(-theta/T[i]) * heavyside_T(T[i] - T_ignition) / (T[i]*(T[i] - 1 + Q*(Z[i] - 1))));
  dxdT = 1/(T[i]-1+Q*(Z[i]-1));
  delta_x = dxdT*delta_T;

  T[i+1] = T[i] + delta_T;
  Y[i+1] = Y[i] + dydT*delta_T;
  Z[i+1] = Z[i] + dzdT*delta_T;

  x[i+1] = x[i] + dxdT*delta_T;
  Y_x[i+1] = Y_x[i] + dydT/dxdT*delta_x;
  T_x[i+1] = T_x[i] + 1/dxdT*delta_x;
  Z_x[i+1] = Z_x[i] + dzdT/dxdT*delta_x;
}

///////////////////////////////////////////////////////////////////////////////
// Bisection for flame speed
///////////////////////////////////////////////////////////////////////////////
inline void RK4_Low_Mach_Solver::bisection_lambda() {
  if (Z[number_of_nodes-2] < 1.0){
    lambda_min = lambda;
    lambda = lambda + (lambda_max-lambda) / 2;
  } else if (Z[number_of_nodes-2] > 1.0){
    lambda_max = lambda;
    lambda = lambda - (lambda - lambda_min)/2;
  }
}

///////////////////////////////////////////////////////////////////////////////
// heavyside_T function
///////////////////////////////////////////////////////////////////////////////
inline double RK4_Low_Mach_Solver::heavyside_T(double x) {
  return static_cast<double>(x >= 0.0);
}

#endif //#endif RK4_Low_Mach_Solver_H
