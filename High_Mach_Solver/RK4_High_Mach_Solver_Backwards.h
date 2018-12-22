#ifndef RK4_HIGH_MACH_SOLVER_BACKWARDS_H
#define RK4_HIGH_MACH_SOLVER_BACKWARDS_H
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Header containg the definition of the class for a
///         RK4_High_Mach_Solver_Backwards
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Gnuplot_CS.h"
#include <math.h>

// using scalar_type = long double;
template <typename scalar_type>
class RK4_High_Mach_Solver_Backwards{

///////////////////////////////////////////////////////////////////////////////
  /// \brief  Header containg a solution to the time invarient Compressible_mach
  ///         Navier-Stokes Equations using Runge-Kutta 4. Based on "Travnikov,
  ///         O., Liberman, & Bychkov. (1997). Stability of a planar flame
  ///         front in a compressible flow. Physics of Fluids, 9(12), 3935-3937."
///////////////////////////////////////////////////////////////////////////////

public:

  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  RK4_High_Mach_Solver_Backwards() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  RK4_High_Mach_Solver_Backwards(const RK4_High_Mach_Solver_Backwards&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  RK4_High_Mach_Solver_Backwards(RK4_High_Mach_Solver_Backwards&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  RK4_High_Mach_Solver_Backwards& operator=(const RK4_High_Mach_Solver_Backwards&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  RK4_High_Mach_Solver_Backwards& operator=(RK4_High_Mach_Solver_Backwards&&) = default;

  RK4_High_Mach_Solver_Backwards(scalar_type le_in, scalar_type Q_in, scalar_type theta_in,
                                 scalar_type gamma_in, scalar_type mf_in, scalar_type Pr_in,
                                 scalar_type U_match, int number_of_nodes_in = 1e8,
                                 scalar_type domaine_length_in = 40000):
                                 number_of_nodes(number_of_nodes_in), domaine_length(domaine_length_in),
                                 Le(le_in), Q(Q_in), theta(theta_in), gamma(gamma_in), mf(mf_in), Pr(Pr_in) {
    delta_x = -domaine_length/static_cast<scalar_type>(number_of_nodes);
    if (Q < 4 || Q > 9) {
      std::cout << "Make sure that Q :" << Q << " is the correct value." << std::endl;
    }
    if (theta > 200) {
      std::cout << "Make sure that theta :" << theta << " is the correct value." << std::endl;
    }
    if (Le < 0.3 || Le > 1.8) {
      std::cout << "Make sure that Le: " << Le << " is the correct value" << std::endl;
    }
    make_reactive_solution();
    cut_data(U_match);
    export_data_to_dat();
  }

  void cut_data(scalar_type U_match){
    std::vector<scalar_type> Theta_cut;
    std::vector<scalar_type> U_cut;
    std::vector<scalar_type> Y_cut;
    std::vector<scalar_type> x_cut;
      for(size_t i = 0; i < Theta_vec.size()-1; ++i){
        if(U_vec[i+1] < U_match){
          break;
        }
        Theta_cut.push_back(Theta_vec[i]);
        U_cut.push_back(U_vec[i]);
        Y_cut.push_back(Y_vec[i]);
        x_cut.push_back(x[i]);
      }
      Theta_vec =  Theta_cut;
      U_vec =  U_cut;
      Y_vec =  Y_cut;
      x =  x_cut;
      double change_in_x = *std::min_element(std::begin(x), std::end(x));
      for (size_t i = 0; i < x.size(); ++i) {
        x[i] -= change_in_x;
      }
      std::vector<scalar_type> x_temp;
      std::vector<scalar_type> U_temp;
      std::vector<scalar_type> Y_temp;
      std::vector<scalar_type> Theta_temp;
      for (size_t i = x.size(); i > 0; --i) {
        x_temp.push_back(x[i-1]);
        U_temp.push_back(U_vec[i-1]);
        Y_temp.push_back(Y_vec[i-1]);
        Theta_temp.push_back(Theta_vec[i-1]);
      }
      x = x_temp;
      U_vec = U_temp;
      Y_vec = Y_temp;
      Theta_vec = Theta_temp;
  }

  scalar_type Power(scalar_type var, int itt){
    scalar_type sol = var;
    for(int i = 1; i < itt; ++i){
      sol *= var;
    }
    return sol;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Returns Y at a location.
  scalar_type get_Y(scalar_type location) {return 0.0;}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Returns Z at a location.
  // scalar_type get_Z(scalar_type location) {return interpolate(Z_x, location);}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Returns T at a location.
  scalar_type get_T(scalar_type location) {return interpolate(Theta_vec, location);}//*(gamma*mf*mf);}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Returns U at a location.
  scalar_type get_U(scalar_type location) {return interpolate(U_vec, location);}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Returns rho at a location.
  scalar_type get_rho(scalar_type location) {return 1 / get_U(location);}

  /////////////////////////////////////////////////////////////////////////
  /// \brief Returns lambda at a location.
  scalar_type get_lambda() {return lambda;}

  /////////////////////////////////////////////////////////////////////////
  /// \brief returns the length of the domaine.
  scalar_type length(){
    auto it_max = std::max_element(std::begin(x), std::end(x));
    auto it_min = std::min_element(std::begin(x), std::end(x));
    return *it_max - *it_min;
    // return fabs(delta_x * Y_vec.size());
  };

  /////////////////////////////////////////////////////////////////////////
  /// \brief Creates a .dat file containing solution to
  ///        the Compressible_mach solution.
  void export_data_to_dat();

  /////////////////////////////////////////////////////////////////////////
  // / \brief Gets the index from a location.
  int get_i(scalar_type location) {
    int i = 0;
    while(((x[i] - location) * (x[i+1] - location)) > 0){
      ++i;
    }
    return i;
    // return static_cast<int>(location/fabs(delta_x));
  }

  /////////////////////////////////////////////////////////////////////////
  // / \brief Gets a value between nodes.
  scalar_type interpolate(std::vector<scalar_type> T, scalar_type location);

  /////////////////////////////////////////////////////////////////////////
  /// \brief Build internal structure of the flame in temperature and
  ///        distance space.
  /// \param Le     Lewis Number.
  /// \param Q      Heat of reaction.
  /// \param theta  Activation Energy.
  /// \param T_ignition_in    Actication Temperature
  void make_reactive_solution();

  /////////////////////////////////////////////////////////////////////////
  /// \brief Sets up boundary conditions.
  void setup_boundary_conditions();

scalar_type dudx_(scalar_type Theta_, scalar_type U_, scalar_type dYdx_, scalar_type Y_);

scalar_type dThetadx_(scalar_type Theta_, scalar_type U_, scalar_type dYdx_, scalar_type Y_);

scalar_type ddyddx_(scalar_type Theta_, scalar_type U_, scalar_type dYdx_, scalar_type Y_);

scalar_type dydx_(scalar_type Theta_, scalar_type U_, scalar_type dYdx_, scalar_type Y_);

  /////////////////////////////////////////////////////////////////////////
  /// \brief One itteration of runge_kutta_4.
  bool runge_kutta_4(int i);

  /////////////////////////////////////////////////////////////////////////
  /// \brief  Modifies lambda, lambda_min and lambda_max using bisection to
  ///         find the actual lambda.
  void bisection_lambda(scalar_type Y);

  void add_lambda_gap();

  scalar_type U2();
  scalar_type Theta2();
  scalar_type Y1();
  // scalar_type eps();

  scalar_type lambda_min = 43151.2;
  scalar_type lambda = 43151.3;
  scalar_type lambda_max = 43151.4;
  int number_of_nodes;
  scalar_type domaine_length;
  scalar_type Le;
  scalar_type Q;
  scalar_type theta;
  scalar_type gamma;
  scalar_type mf;
  scalar_type Pr;
  std::vector<scalar_type> Y_vec;
  std::vector<scalar_type> U_vec;
  std::vector<scalar_type> Theta_vec;
  scalar_type Y;
  scalar_type U;
  scalar_type Theta;
  scalar_type dYdx;
  std::vector<scalar_type>  x;
  scalar_type delta_x;
  // scalar_type accuracy = 1e-10;
private:
  bool check = false;
  bool old_check1 = false;
  bool old_check2 = true;
  bool old_check3 = false;
};

///////////////////////////////////////////////////////////////////////////////
//Gets a value between nodes.
///////////////////////////////////////////////////////////////////////////////
template <typename scalar_type>
scalar_type RK4_High_Mach_Solver_Backwards<scalar_type>::interpolate(std::vector<scalar_type> T, scalar_type location) {
  int i = get_i(location);
  return (T[i+1] - T[i]) * (location - x[i]) / (x[i+1] - x[i]) + T[i];
}

///////////////////////////////////////////////////////////////////////////////
// Creates a .dat file containing solution to the Compressible_mach solution.
///////////////////////////////////////////////////////////////////////////////
template <typename scalar_type>
void RK4_High_Mach_Solver_Backwards<scalar_type>::export_data_to_dat() {
  std::vector<std::vector<scalar_type>> data;

  std::vector<std::vector<scalar_type>> data_usual;
  // re_position_positive();

  // std::vector<scalar_type> x(number_of_nodes);
  // for(int i = 0; i < number_of_nodes; ++i){
  //   x[i] = delta_x*i;
  // }
  std::vector<scalar_type> rho_usual(x.size());
  for(size_t i = 0; i < x.size(); ++i){
    rho_usual[i] = 1.0/U_vec[i];
  }
  std::vector<scalar_type> T_usual(x.size());
  for(size_t i = 0; i < x.size(); ++i){
    T_usual[i] = Theta_vec[i]/(gamma*mf*mf);
  }
  std::vector<scalar_type> p_usual(x.size());
  for(size_t i = 0; i < x.size(); ++i){
    p_usual[i] = T_usual[i]*rho_usual[i];
  }
  std::vector<scalar_type> Y_usual(x.size());
  for(size_t i = 0; i < x.size(); ++i){
    Y_usual[i] = 0;//Y_vec[i]/Y_vec[0];
  }

  data = {x, U_vec, Theta_vec, Y_vec};//, U, Theta, dYdx};
  data_usual = {x,rho_usual, U_vec, p_usual, T_usual, Y_usual};
  std::string variables = "# x Y";
  std::string variables_usual = "# x rho U p T Y";
  gnuplot_CS<scalar_type>("RK4_Burned", data, variables);
  gnuplot_CS<scalar_type>("RK4_Usual_Burned", data_usual, variables_usual);
}

///////////////////////////////////////////////////////////////////////////////
// Runs all required functions to solve the pde's
///////////////////////////////////////////////////////////////////////////////
template <typename scalar_type>
void RK4_High_Mach_Solver_Backwards<scalar_type>::make_reactive_solution() {
    x.clear();
    Y_vec.clear();
    U_vec.clear();
    Theta_vec.clear();
    setup_boundary_conditions();
    bool donecheck = 0;
    for(int i = 0; i < number_of_nodes; ++i) {
      donecheck = runge_kutta_4(i);
      if(donecheck == 1){
        if(!std::isfinite(Y)){
          bisection_lambda(100);
        }
        break;
      }
    }
    if(donecheck == 0){
      bisection_lambda(100.0);
    }
    export_data_to_dat();
}

///////////////////////////////////////////////////////////////////////////////
// setup all boundary conditions for the shooting method
///////////////////////////////////////////////////////////////////////////////
template <typename scalar_type>
scalar_type RK4_High_Mach_Solver_Backwards<scalar_type>::U2(){

  scalar_type zeta = (1.0 + mf*mf*mf*mf - 2.0 * mf*mf*(1.0+Q+gamma*Q)) / (pow(pow(mf, 2)*(1.0+gamma),2));
  scalar_type rho_inf = 1.0/((1.0 + pow(mf,2)*gamma)/(pow(mf,2)*(1.0 + gamma)) - sqrt(zeta));
  return 1.0/rho_inf;

}
template <typename scalar_type>
scalar_type RK4_High_Mach_Solver_Backwards<scalar_type>::Theta2() {
  scalar_type p_0 = 1/(gamma*mf*mf);
  scalar_type zeta = (1.0 + mf*mf*mf*mf - 2.0 * mf*mf*(1.0+Q+gamma*Q)) / (pow(pow(mf, 2)*(1.0+gamma),2));
  scalar_type rho_inf = 1.0/((1.0 + pow(mf,2)*gamma)/(pow(mf,2)*(1.0 + gamma)) - sqrt(zeta));
  scalar_type p_inf = p_0*((1.0 + pow(mf,2)*gamma*(1.0-1.0/rho_inf)));
  return p_inf/p_0/rho_inf;
}

template <typename scalar_type>
scalar_type RK4_High_Mach_Solver_Backwards<scalar_type>::Y1() {
 return (gamma-1)*0.5*mf*mf*(U2()*U2()-1)+Theta2()-1;
}

template <typename scalar_type>
void RK4_High_Mach_Solver_Backwards<scalar_type>::setup_boundary_conditions() {
  Y = 0.0;//Y1();
  U = 0.9999 * U2();//1.0;
  Theta = 1.0001 * Theta2();//1.0;
  dYdx =  0.0;//-1e-15;
}

template <typename scalar_type>
scalar_type RK4_High_Mach_Solver_Backwards<scalar_type>::dudx_(scalar_type Theta_, scalar_type U_, scalar_type, scalar_type) {
  return 0.75/(Pr*gamma*mf*mf)*(Theta_/U_-1+gamma*mf*mf*(U_-1));
}

template <typename scalar_type>
scalar_type RK4_High_Mach_Solver_Backwards<scalar_type>::dThetadx_(scalar_type Theta_, scalar_type U_, scalar_type dYdx_, scalar_type Y_) {
  return -1.0/Le*dYdx_+(gamma-1)/gamma*(U_-Theta_-gamma*mf*mf/2.0*(U_-1)*(U_-1))+Theta_-1+Y_-Y1();
}
template <typename scalar_type>
scalar_type RK4_High_Mach_Solver_Backwards<scalar_type>::ddyddx_(scalar_type Theta_, scalar_type U_, scalar_type dYdx_, scalar_type Y_) {
  // if(fabs(Y_)<accuracy){
  //   return 0.0;
  // }
  return Le*(dYdx_+lambda*Y_/U_*exp(-theta/Theta_));
}

template <typename scalar_type>
scalar_type RK4_High_Mach_Solver_Backwards<scalar_type>::dydx_(scalar_type Theta_, scalar_type U_, scalar_type dYdx_, scalar_type Y_) {
  // if(fabs(Y_)<accuracy){
  //   return 0.0;
  // }
  return dYdx_ + delta_x * ddyddx_(Theta_, U_, dYdx_, Y_);
}

///////////////////////////////////////////////////////////////////////////////
// Runge kutta 4
///////////////////////////////////////////////////////////////////////////////
template <typename scalar_type>
bool RK4_High_Mach_Solver_Backwards<scalar_type>::runge_kutta_4(int i) {

  scalar_type F1dYdx =   delta_x*ddyddx_(  Theta, U, dYdx, Y);
  scalar_type F1Y =      delta_x*dydx_(    Theta, U, dYdx, Y);
  scalar_type F1U =      delta_x*dudx_(    Theta, U, dYdx, Y);
  scalar_type F1Theta =  delta_x*dThetadx_(Theta, U, dYdx, Y);

  scalar_type F2dYdx =   delta_x*ddyddx_(  Theta+0.5*F1Theta, U+0.5*F1U, dYdx+0.5*F1dYdx, Y+0.5*F1Y);
  scalar_type F2Y =      delta_x*dydx_(    Theta+0.5*F1Theta, U+0.5*F1U, dYdx+0.5*F1dYdx, Y+0.5*F1Y);
  scalar_type F2U =      delta_x*dudx_(    Theta+0.5*F1Theta, U+0.5*F1U, dYdx+0.5*F1dYdx, Y+0.5*F1Y);
  scalar_type F2Theta =  delta_x*dThetadx_(Theta+0.5*F1Theta, U+0.5*F1U, dYdx+0.5*F1dYdx, Y+0.5*F1Y);

  scalar_type F3dYdx =   delta_x*ddyddx_(  Theta+0.5*F2Theta, U+0.5*F2U, dYdx+0.5*F2dYdx, Y+0.5*F2Y);
  scalar_type F3Y =      delta_x*dydx_(    Theta+0.5*F2Theta, U+0.5*F2U, dYdx+0.5*F2dYdx, Y+0.5*F2Y);
  scalar_type F3U =      delta_x*dudx_(    Theta+0.5*F2Theta, U+0.5*F2U, dYdx+0.5*F2dYdx, Y+0.5*F2Y);
  scalar_type F3Theta =  delta_x*dThetadx_(Theta+0.5*F2Theta, U+0.5*F2U, dYdx+0.5*F2dYdx, Y+0.5*F2Y);

  scalar_type F4dYdx =   delta_x*ddyddx_(  Theta+ F3Theta, U+F3U, dYdx+F3dYdx, Y+ F3Y);
  scalar_type F4Y =      delta_x*dydx_(    Theta+ F3Theta, U+F3U, dYdx+F3dYdx, Y+ F3Y);
  scalar_type F4U =      delta_x*dudx_(    Theta+ F3Theta, U+F3U, dYdx+F3dYdx, Y+ F3Y);
  scalar_type F4Theta =  delta_x*dThetadx_(Theta+ F3Theta, U+F3U, dYdx+F3dYdx, Y+ F3Y);


  dYdx = dYdx + 1.0/6.0*(F1dYdx + 2.0*F2dYdx + 2.0*F3dYdx + F4dYdx);
  Y = Y + 1.0/6.0*(F1Y + 2.0*F2Y + 2.0*F3Y + F4Y);
  U = U + 1.0/6.0*(F1U + 2.0*F2U + 2.0*F3U + F4U);
  Theta = Theta + 1.0/6.0*(F1Theta + 2.0*F2Theta + 2.0*F3Theta + F4Theta);

  if((i%(number_of_nodes/100000) == 0)){
    Y_vec.push_back(Y);
    Theta_vec.push_back(Theta);
    U_vec.push_back(U);
    x.push_back(i*delta_x);
  }
  return 0;
}


///////////////////////////////////////////////////////////////////////////////
// Bisection for flame speed
///////////////////////////////////////////////////////////////////////////////
template <typename scalar_type>
inline void RK4_High_Mach_Solver_Backwards<scalar_type>::bisection_lambda(scalar_type Y) {
  if (Y < 0.0){
    lambda_min = lambda;
    lambda = lambda + (lambda_max-lambda) / 2;
    check = 1;
  } else {
    lambda_max = lambda;
    lambda = lambda - (lambda - lambda_min)/2;
    check = 0;
  }
  add_lambda_gap();
}

template <typename scalar_type>
void RK4_High_Mach_Solver_Backwards<scalar_type>::add_lambda_gap() {
  if(check == old_check1 && check == old_check2 && check == old_check3){
    if (check == 0) {
      lambda_min -= (lambda_max-lambda_min);
    } else {
      lambda_max += (lambda_max-lambda_min);
    }
  }
  old_check3 = old_check2;
  old_check2 = old_check1;
  old_check1 = check;
}

#endif //#endif RK4_High_Mach_Solver_H
