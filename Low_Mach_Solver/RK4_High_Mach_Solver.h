#ifndef RK4_HIGH_MACH_SOLVER_H
#define RK4_HIGH_MACH_SOLVER_H
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Header containg the definition of the class for a
///         RK4_High_Mach_Solver
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Gnuplot_RK4.h"

class RK4_High_Mach_Solver{

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
  RK4_High_Mach_Solver() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  RK4_High_Mach_Solver(const RK4_High_Mach_Solver&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  RK4_High_Mach_Solver(RK4_High_Mach_Solver&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  RK4_High_Mach_Solver& operator=(const RK4_High_Mach_Solver&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  RK4_High_Mach_Solver& operator=(RK4_High_Mach_Solver&&) = default;

  RK4_High_Mach_Solver(long double le_in, long double Q_in, long double theta_in, long double gamma_in, long double mf_in, long double Pr_in) :
    number_of_nodes(1e7), domaine_length(50), Le(le_in), Q(Q_in), theta(theta_in), gamma(gamma_in), mf(mf_in), Pr(Pr_in) {
    delta_x = domaine_length/static_cast<long double>(number_of_nodes);
    std::cout << "delta_x: " << delta_x << std::endl;
    if (Q < 4 || Q > 9) {
      std::cout << "Make sure that Q :" << Q << " is the correct value." << std::endl;
    }
    if (theta > 200) {
      std::cout << "Make sure that theta :" << theta << " is the correct value." << std::endl;
    }
    if (Le < 0.3 || Le > 1.8) {
      std::cout << "Make sure that Le: " << Le << " is the correct value" << std::endl;
    }
    std::cout << "passed warnings" << std::endl;
    resize_arrays();
    std::cout << "resized arrays" << std::endl;
    make_low_mach_reactive_solution();
    // re_position_positive();
    export_data_to_dat();
  }

  // /////////////////////////////////////////////////////////////////////////
  // /// \brief Returns Y at a location.
  // long double get_Y(long double location) {return interpolate(Y_x, location);}
  //
  // /////////////////////////////////////////////////////////////////////////
  // /// \brief Returns Z at a location.
  // long double get_Z(long double location) {return interpolate(Z_x, location);}
  //
  // /////////////////////////////////////////////////////////////////////////
  // /// \brief Returns T at a location.
  // long double get_T(long double location) {return interpolate(T_x, location);}
  //
  // /////////////////////////////////////////////////////////////////////////
  // /// \brief Returns U at a location.
  // long double get_U(long double location) {return get_T(location);}
  //
  // /////////////////////////////////////////////////////////////////////////
  // /// \brief Returns rho at a location.
  // long double get_rho(long double location) {return 1 / get_T(location);}
  //
  // /////////////////////////////////////////////////////////////////////////
  // /// \brief Returns lambda at a location.
  // long double get_lambda() {return lambda;}
  //
  // /////////////////////////////////////////////////////////////////////////
  // /// \brief returns the length of the domaine.
  // long double length();
  //
  // /////////////////////////////////////////////////////////////////////////
  // /// \brief Positions the solution starting at x = 0;
  // void re_position_positive();
  //
  // /////////////////////////////////////////////////////////////////////////
  // /// \brief Creates a .dat file containing solution to
  // ///        the low_mach solution.
  void export_data_to_dat();
private:

  /////////////////////////////////////////////////////////////////////////
  /// \brief Gets the index from a location.
  // int get_i(long double location);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Gets a value between nodes.
  // long double interpolate(std::vector<long double> T, long double location);

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
  /// \brief Sets up boundary conditions.
  void setup_boundary_conditions();

long double dudx_(long double Theta_, long double U_, long double dYdx_, long double Y_);

long double dThetadx_(long double Theta_, long double U_, long double dYdx_, long double Y_);

long double ddyddx_(long double Theta_, long double U_, long double dYdx_, long double Y_);

long double dydx_(long double Theta_, long double U_, long double dYdx_, long double Y_);

  /////////////////////////////////////////////////////////////////////////
  /// \brief One itteration of runge_kutta_4.
  //TODO:(Andre) Seperate Distance space for only the final solution.
  void runge_kutta_4(int i);

  /////////////////////////////////////////////////////////////////////////
  /// \brief  Modifies lambda, lambda_min and lambda_max using bisection to
  ///         find the actual lambda.
  void bisection_lambda(long double Y);

  long double U2();
  long double Theta2();
  long double Y1();
  // long double eps();

  long double lambda_min = 3e2;
  long double lambda = 6.01477275136390624829800799489021301269531250000000e+07;
  long double lambda_max = 3e010;
  int number_of_nodes;
  std::vector<long double> Y;
  std::vector<long double> U;
  std::vector<long double> Theta;
  std::vector<long double> dYdx;
  long double Le;
  long double Q;
  long double theta;
  long double Pr;
  long double delta_x;
  long double domaine_length;
  long double gamma;
  long double mf;
  long double accuracy = 1e-20;
};

///////////////////////////////////////////////////////////////////////////////
//Gets a value between nodes.
///////////////////////////////////////////////////////////////////////////////
// long double RK4_High_Mach_Solver::interpolate(std::vector<long double> T, long double location) {
//   int i = get_i(location);
//   return (T[i+1] - T[i]) * (location - x[i]) / (x[i+1] - x[i]) + T[i];
// }

///////////////////////////////////////////////////////////////////////////////
//Gets the index from a location.
///////////////////////////////////////////////////////////////////////////////
// int RK4_High_Mach_Solver::get_i(long double location) {
//   int i = 0;
//   while(((x[i] - location) * (x[i+1] - location)) > 0){
//     ++i;
//   }
//   return i;
// }

///////////////////////////////////////////////////////////////////////////////
// Creates a .dat file containing solution to the low_mach solution.
///////////////////////////////////////////////////////////////////////////////
void RK4_High_Mach_Solver::export_data_to_dat() {
  std::vector<std::vector<long double>> data;
  std::vector<std::vector<long double>> data_usual;
  std::vector<long double> x(number_of_nodes);
  for(int i = 0; i < number_of_nodes; ++i){
    x[i] = delta_x*i;
  }
  std::vector<long double> rho(number_of_nodes);
  for(int i = 0; i < number_of_nodes; ++i){
    rho[i] = 1.0/U[i];
  }
  std::vector<long double> T(number_of_nodes);
  for(int i = 0; i < number_of_nodes; ++i){
    T[i] = Theta[i]/(gamma*mf*mf);
  }
  std::vector<long double> p(number_of_nodes);
  for(int i = 0; i < number_of_nodes; ++i){
    p[i] = T[i]*rho[i];
  }

  data = {x, Y, U, Theta, dYdx};
  data_usual = {x,rho, U, p, T};
  std::string variables = "# x Y U Theta dydx";
  std::string variables_usual = "# x rho U p T";
  gnuplot_RK4<long double>("RK4", data, variables);
  gnuplot_RK4<long double>("RK4_Usual", data_usual, variables_usual);
}

///////////////////////////////////////////////////////////////////////////////
// Positions the solution starting at x = 0;
///////////////////////////////////////////////////////////////////////////////
// void RK4_High_Mach_Solver::re_position_positive() {
//   long double change_in_x = *std::min_element(std::begin(x), std::end(x));
//   for (size_t i = 0; i < x.size(); ++i) {
//     x[i] -= change_in_x;
//   }
// }


///////////////////////////////////////////////////////////////////////////////
// return length of domaine.
///////////////////////////////////////////////////////////////////////////////
// long double RK4_High_Mach_Solver::length() {
//   auto it_max = std::max_element(std::begin(x), std::end(x));
//   auto it_min = std::min_element(std::begin(x), std::end(x));
//   return *it_max - *it_min;
// }

///////////////////////////////////////////////////////////////////////////////
// Resize the arrays for the number of nodes in the solution.
///////////////////////////////////////////////////////////////////////////////
void RK4_High_Mach_Solver::resize_arrays() {
  Y.resize(number_of_nodes-1);
  U.resize(number_of_nodes-1);
  Theta.resize(number_of_nodes-1);
  dYdx.resize(number_of_nodes-1);
}

///////////////////////////////////////////////////////////////////////////////
// Runs all required functions to solve the pde's
///////////////////////////////////////////////////////////////////////////////
void RK4_High_Mach_Solver::make_low_mach_reactive_solution() {
  // while((lambda_max - lambda_min) > 1e-20){

    setup_boundary_conditions();
    for(int i = 0; i < number_of_nodes-2; ++i) {
      // runge_kutta_4(0);
      runge_kutta_4(i);
    }
    int temp = 0;
    for(int i = 0; i < number_of_nodes-2; ++i) {
      if((Y[i] < 0 || Y[i] > 9.1) && temp == 0){
        temp = 1;
        bisection_lambda(Y[i]);
      }
    }
    if(temp == 0){
      bisection_lambda(100.0);
    }
  // }
}

///////////////////////////////////////////////////////////////////////////////
// setup all boundary conditions for the shooting method
///////////////////////////////////////////////////////////////////////////////
long double RK4_High_Mach_Solver::U2(){

  long double zeta = (1.0 + mf*mf*mf*mf - 2.0 * mf*mf*(1.0+Q+gamma*Q)) / (pow(pow(mf, 2)*(1.0+gamma),2));
  long double rho_inf = 1.0/((1.0 + pow(mf,2)*gamma)/(pow(mf,2)*(1.0 + gamma)) - sqrt(zeta));

  return 1.0/rho_inf;

}

long double RK4_High_Mach_Solver::Theta2() {
  long double p_0 = 1/(gamma*mf*mf);
  long double zeta = (1.0 + mf*mf*mf*mf - 2.0 * mf*mf*(1.0+Q+gamma*Q)) / (pow(pow(mf, 2)*(1.0+gamma),2));
  long double rho_inf = 1.0/((1.0 + pow(mf,2)*gamma)/(pow(mf,2)*(1.0 + gamma)) - sqrt(zeta));
  long double p_inf = p_0*((1.0 + pow(mf,2)*gamma*(1.0-1.0/rho_inf)));

  return p_inf/p_0/rho_inf;
}

long double RK4_High_Mach_Solver::Y1() {
 return (gamma-1)*0.5*mf*mf*(U2()*U2()-1)+Theta2()-1;
}

void RK4_High_Mach_Solver::setup_boundary_conditions() {
  Y[0] = Y1();
  U[0] = 1;
  Theta[0] =1;
  dYdx[0] = -1e-5;
}

long double RK4_High_Mach_Solver::dudx_(long double Theta_, long double U_, long double dYdx_, long double Y_) {
  return 0.75/(Pr*gamma*mf*mf)*(Theta_/U_-1+gamma*mf*mf*(U_-1));
}

long double RK4_High_Mach_Solver::dThetadx_(long double Theta_, long double U_, long double dYdx_, long double Y_) {
  return -1.0/Le*dYdx_+(gamma-1)/gamma*(U_-Theta_-gamma*mf*mf/2.0*(U_-1)*(U_-1))+Theta_-1+Y_-Y1();
}

long double RK4_High_Mach_Solver::ddyddx_(long double Theta_, long double U_, long double dYdx_, long double Y_) {
  if(fabs(Y_)<accuracy){
    return 0.0;
  }
  return Le*(dYdx_+lambda*Y_/U_*exp(-theta/Theta_));
}

long double RK4_High_Mach_Solver::dydx_(long double Theta_, long double U_, long double dYdx_, long double Y_) {
  if(fabs(Y_)<accuracy){
    return 0.0;
  }
  return dYdx_ + delta_x * ddyddx_(Theta_, U_, dYdx_, Y_);
}

///////////////////////////////////////////////////////////////////////////////
// Runge kutta 4
///////////////////////////////////////////////////////////////////////////////
void RK4_High_Mach_Solver::runge_kutta_4(int i) {

  long double F1dYdx =   delta_x*ddyddx_(  Theta[i], U[i], dYdx[i], Y[i]);
  long double F1Y =      delta_x*dydx_(    Theta[i], U[i], dYdx[i], Y[i]);
  long double F1U =      delta_x*dudx_(    Theta[i], U[i], dYdx[i], Y[i]);
  long double F1Theta =  delta_x*dThetadx_(Theta[i], U[i], dYdx[i], Y[i]);

  long double F2dYdx =   delta_x*ddyddx_(  Theta[i]+0.5*F1Theta, U[i]+0.5*F1U, dYdx[i]+0.5*F1dYdx, Y[i]+0.5*F1Y);
  long double F2Y =      delta_x*dydx_(    Theta[i]+0.5*F1Theta, U[i]+0.5*F1U, dYdx[i]+0.5*F1dYdx, Y[i]+0.5*F1Y);
  long double F2U =      delta_x*dudx_(    Theta[i]+0.5*F1Theta, U[i]+0.5*F1U, dYdx[i]+0.5*F1dYdx, Y[i]+0.5*F1Y);
  long double F2Theta =  delta_x*dThetadx_(Theta[i]+0.5*F1Theta, U[i]+0.5*F1U, dYdx[i]+0.5*F1dYdx, Y[i]+0.5*F1Y);

  long double F3dYdx =   delta_x*ddyddx_(  Theta[i]+0.5*F2Theta, U[i]+0.5*F2U, dYdx[i]+0.5*F2dYdx, Y[i]+0.5*F2Y);
  long double F3Y =      delta_x*dydx_(    Theta[i]+0.5*F2Theta, U[i]+0.5*F2U, dYdx[i]+0.5*F2dYdx, Y[i]+0.5*F2Y);
  long double F3U =      delta_x*dudx_(    Theta[i]+0.5*F2Theta, U[i]+0.5*F2U, dYdx[i]+0.5*F2dYdx, Y[i]+0.5*F2Y);
  long double F3Theta =  delta_x*dThetadx_(Theta[i]+0.5*F2Theta, U[i]+0.5*F2U, dYdx[i]+0.5*F2dYdx, Y[i]+0.5*F2Y);

  long double F4dYdx =   delta_x*ddyddx_(  Theta[i]+ F3Theta, U[i]+F3U, dYdx[i]+F3dYdx, Y[i]+ F3Y);
  long double F4Y =      delta_x*dydx_(    Theta[i]+ F3Theta, U[i]+F3U, dYdx[i]+F3dYdx, Y[i]+ F3Y);
  long double F4U =      delta_x*dudx_(    Theta[i]+ F3Theta, U[i]+F3U, dYdx[i]+F3dYdx, Y[i]+ F3Y);
  long double F4Theta =  delta_x*dThetadx_(Theta[i]+ F3Theta, U[i]+F3U, dYdx[i]+F3dYdx, Y[i]+ F3Y);


  dYdx[i+1] = dYdx[i] + 1.0/6.0*(F1dYdx + 2.0*F2dYdx + 2*F3dYdx + F4dYdx);
  if(fabs(Y[i])<accuracy){
    dYdx[i+1] = 0.0;
  }
  Y[i+1] = Y[i] + 1.0/6.0*(F1Y + 2.0*F2Y + 2*F3Y + F4Y);
  if(fabs(Y[i])<accuracy){
    Y[i+1] = 0.0;
  }
  U[i+1] = U[i] + 1.0/6.0*(F1U + 2.0*F2U + 2*F3U + F4U);
  Theta[i+1] = Theta[i] + 1.0/6.0*(F1Theta + 2.0*F2Theta + 2*F3Theta + F4Theta);
}


///////////////////////////////////////////////////////////////////////////////
// Bisection for flame speed
///////////////////////////////////////////////////////////////////////////////
inline void RK4_High_Mach_Solver::bisection_lambda(long double Y) {
  if (Y < 0.0){
    lambda_min = lambda;
    lambda = lambda + (lambda_max-lambda) / 2;
  } else if (Y > 9.1){
    lambda_max = lambda;
    lambda = lambda - (lambda - lambda_min)/2;
  }
  std::cout << "lambda: " << lambda << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// heavyside_T function
///////////////////////////////////////////////////////////////////////////////
// inline long double RK4_High_Mach_Solver::heavyside_T(long double x) {
//   return static_cast<long double>(x >= 0.0);
// }

#endif //#endif RK4_High_Mach_Solver_H
