#ifndef IMPLICIT_HLLE_H
#define IMPLICIT_HLLE_H

#include "HLLE.h"

template <typename grid_type, typename flow_type>
class Implicit_HLLE {
using global_solution_vector_type = typename grid_type::global_solution_vector_type;
using solution_vector_type = typename global_solution_vector_type::value_type;
using scalar_type = typename grid_type::scalar_type;
using size_type = typename grid_type::size_type;
using matrix_type = typename grid_type::matrix_type;
 public:
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  Implicit_HLLE() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  Implicit_HLLE(const Implicit_HLLE&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  Implicit_HLLE(Implicit_HLLE&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  Implicit_HLLE& operator=(const Implicit_HLLE&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  Implicit_HLLE& operator=(Implicit_HLLE&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Constructor setting up required inputs.
  Implicit_HLLE( const grid_type grid, const flow_type flow,
                            const global_solution_vector_type& DeltaUm, const scalar_type dt_in,
                            const scalar_type Theta_in, const scalar_type zeta_in,
                            const size_type cell_index) :
                            gamma(flow.gamma), Pr(flow.Pr), Le(flow.Le), Q(flow.Q()), lambda(flow.lambda),
                            theta(flow.theta()), dx(grid.dx()), dt(dt_in), zeta(zeta_in),
                            Theta(Theta_in), T_ignition(flow.T_ignition()) {

  solution_vector_mm  = grid.global_solution_vector[std::max(static_cast<int>(cell_index)-2 , 0)];
  solution_vector_m   = grid.global_solution_vector[cell_index-1];
  solution_vector     = grid.global_solution_vector[cell_index];
  solution_vector_p   = grid.global_solution_vector[cell_index+1];
  solution_vector_pp  = grid.global_solution_vector[std::min(cell_index+2, grid.number_of_cells()-1)];
  rhomm = solution_vector_mm[0];
  rhom = solution_vector_m[0];
  rho = solution_vector[0];
  rhop = solution_vector_p[0];
  rhopp = solution_vector_pp[0];
  umm = solution_vector_mm[1]/solution_vector_mm[0];
  um = solution_vector_m[1]/solution_vector_m[0];
  u = solution_vector[1]/solution_vector[0];
  up = solution_vector_p[1]/solution_vector_p[0];
  upp = solution_vector_pp[1]/solution_vector_pp[0];
  emm = solution_vector_mm[2];
  em = solution_vector_m[2];
  e = solution_vector[2];
  ep = solution_vector_p[2];
  epp = solution_vector_pp[2];
  Ymm = solution_vector_mm[3]/solution_vector_mm[0];
  Ym = solution_vector_m[3]/solution_vector_m[0];
  Y = solution_vector[3]/solution_vector[0];
  Yp = solution_vector_p[3]/solution_vector_p[0];
  Ypp = solution_vector_pp[3]/solution_vector_pp[0];
  pm = (solution_vector_m[2] - solution_vector_m[1]* solution_vector_m[1] /
       (2.0 * solution_vector_m[0])) * (gamma - 1);
  p =  (solution_vector[2] - solution_vector[1]* solution_vector[1] /
       (2.0 * solution_vector[0])) * (gamma - 1);
  pp = (solution_vector_p[2] - solution_vector_p[1]* solution_vector_p[1] /
       (2.0 * solution_vector_p[0])) * (gamma - 1);
  dUm1 = DeltaUm[cell_index-1][0];
  dUm2 = DeltaUm[cell_index-1][1];
  dUm3 = DeltaUm[cell_index-1][2];
  dUm4 = DeltaUm[cell_index-1][3];

  {
  scalar_type lambda_1 = u + sqrt(gamma*p/rho);
  scalar_type lambda_2 = u_rhoavg(rhom, rho, um, u) + sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  LlambdaR = std::max(lambda_1,lambda_2);
  }
  {
  scalar_type lambda_1 = um - sqrt(gamma*pm/rhom);
  scalar_type lambda_2 = u_rhoavg(rhom, rho, um, u) - sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  LlambdaL = std::min(lambda_1,lambda_2);
  }
  {
  scalar_type lambda_1 = up + sqrt(gamma*pp/rhop);
  scalar_type lambda_2 = u_rhoavg(rho, rhop, u, up) + sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  RlambdaR = std::max(lambda_1,lambda_2);
  }
  {
  scalar_type lambda_1 = u - sqrt(gamma*p/rho);
  scalar_type lambda_2 = u_rhoavg(rho, rhop, u, up) - sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  RlambdaL = std::min(lambda_1,lambda_2);
  }
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  matrix_type top_matrix();

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  matrix_type mid_matrix();

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  matrix_type bot_matrix();

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  solution_vector_type rhs_matrix();

 private:
  solution_vector_type solution_vector_mm;
  solution_vector_type solution_vector_m;
  solution_vector_type solution_vector;
  solution_vector_type solution_vector_p;
  solution_vector_type solution_vector_pp;
  scalar_type rhomm, rhom, rho, rhop, rhopp;
  scalar_type pmm, pm, p, pp, ppp;
  scalar_type umm, um, u, up, upp;
  scalar_type emm, em, e, ep, epp;
  scalar_type Ymm, Ym, Y, Yp, Ypp;
  char E;
  const scalar_type gamma, Pr, Le, Q;
  const scalar_type lambda, theta;
  const scalar_type dx, dt;
  const scalar_type zeta, Theta;
  scalar_type dUm1, dUm2, dUm3, dUm4;
  scalar_type LlambdaR;
  scalar_type LlambdaL;
  scalar_type RlambdaL;
  scalar_type RlambdaR;
  const scalar_type T_ignition;

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  template <typename T>
  scalar_type Power(const T& num, const int& expo) {
    return pow(num, expo);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  scalar_type Power(const char&, const scalar_type& expo) {
    return exp(expo);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  scalar_type u_rhoavg(scalar_type rho_left, scalar_type rho_right, scalar_type u_left, scalar_type u_right) {
    return (sqrt(rho_left) * u_left + sqrt(rho_right) * u_right) /
    (sqrt(rho_left) + sqrt(rho_right));
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  scalar_type rho_rhoavg(scalar_type rho_left, scalar_type rho_right) {
    return sqrt(rho_left * rho_right);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  scalar_type p_rhoavg(scalar_type rho_left, scalar_type rho_right, scalar_type u_left, scalar_type u_right,
                  scalar_type p_left, scalar_type p_right, scalar_type gamma) {

    scalar_type h_left = gamma * p_left / (rho_left * (gamma - 1.0)) + u_left * u_left * 0.5;
    scalar_type h_right = gamma * p_right / (rho_right * (gamma - 1.0)) + u_right * u_right * 0.5;

    scalar_type h_ravg = (sqrt(rho_left) * h_left + sqrt(rho_right) * h_right) /
    (sqrt(rho_left) + sqrt(rho_right));

    scalar_type u_ravg = u_rhoavg(rho_left, rho_right, u_left, u_right);

    scalar_type rho_ravg = rho_rhoavg(rho_left, rho_right);

    return (h_ravg - u_ravg*u_ravg * 0.5) * (gamma - 1.0) / gamma * rho_ravg;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  solution_vector_type limiter(const solution_vector_type &Ul,
                               const solution_vector_type &U,
                               const solution_vector_type &Ur,
                               scalar_type dx) {
  // solution_vector_type phi;
  // solution_vector_type a = (U - Ul) / dx;
  // solution_vector_type b = (Ur - U) / dx;
  // for(size_type i = 0; i < 4; ++i) {
  // phi[i] = sign(a[i])*std::max(0.0,static_cast<scalar_type>(std::min(fabs(a[i]), sign(a[i])*b[i])));
  // }
  //   return phi;
  solution_vector_type a = (U - Ul) / dx;
  solution_vector_type b = (Ur - U) / dx;
  double epsilon = 1.0e-6;
  solution_vector_type phi =  a.array() * b.array() * (a.array() + b.array()) / (a.array() * a.array() + b.array() * b.array() + epsilon);
  for (int i = 0; i < 4; ++i) {
    if (a[i] / b[i] <= 0.0 || b[i] == 0) {
      phi[i] = 0.0;
    }
  }
  return phi;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  template <typename T> int sign(T val) {return (T(0) < val) - (val < T(0));}

};

template <typename grid_type, typename flow_type>
typename grid_type::matrix_type Implicit_HLLE<grid_type, flow_type>::mid_matrix() {

  matrix_type b;
  matrix_type temp;
///////////////////////////////////////////////////////////////////////////////
// Needed
///////////////////////////////////////////////////////////////////////////////
  b << 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1;
///////////////////////////////////////////////////////////////////////////////
// Hyperbolic (Euler)
///////////////////////////////////////////////////////////////////////////////
#if defined(HYPERBOLIC)

  temp << -((dt*((LlambdaL*LlambdaR)/(-LlambdaL + LlambdaR) + (RlambdaL*RlambdaR)/(-RlambdaL + RlambdaR))*Theta)/(dx*(1 + zeta))),
   (dt*(-(LlambdaL*RlambdaL) + LlambdaR*RlambdaR)*Theta)/(dx*(LlambdaL - LlambdaR)*(RlambdaL - RlambdaR)*(1 + zeta)),0,0,
   -(dt*(-3 + gamma)*(LlambdaL*RlambdaL - LlambdaR*RlambdaR)*Theta*Power(u,2))/(2.*dx*(LlambdaL - LlambdaR)*(RlambdaL - RlambdaR)*(1 + zeta)),
   (dt*Theta*(-(LlambdaR*RlambdaR*(RlambdaL + (-3 + gamma)*u)) + LlambdaL*(LlambdaR*(RlambdaL - RlambdaR) + RlambdaL*(RlambdaR + (-3 + gamma)*u))))/
    (dx*(LlambdaL - LlambdaR)*(RlambdaL - RlambdaR)*(1 + zeta)),
   -((dt*(-1 + gamma)*(LlambdaL*RlambdaL - LlambdaR*RlambdaR)*Theta)/(dx*(LlambdaL - LlambdaR)*(RlambdaL - RlambdaR)*(1 + zeta))),0,
   (dt*(LlambdaL*RlambdaL - LlambdaR*RlambdaR)*Theta*(e*gamma*u - (-1 + gamma)*rho*Power(u,3)))/
    (dx*(LlambdaL - LlambdaR)*rho*(RlambdaL - RlambdaR)*(1 + zeta)),
   (dt*(LlambdaL*RlambdaL - LlambdaR*RlambdaR)*Theta*(-2*e*gamma + 3*(-1 + gamma)*rho*Power(u,2)))/
    (2.*dx*(LlambdaL - LlambdaR)*rho*(RlambdaL - RlambdaR)*(1 + zeta)),
   (dt*Theta*(LlambdaL*LlambdaR*(RlambdaL - RlambdaR) + LlambdaL*RlambdaL*(RlambdaR - gamma*u) + LlambdaR*RlambdaR*(-RlambdaL + gamma*u)))/
    (dx*(LlambdaL - LlambdaR)*(RlambdaL - RlambdaR)*(1 + zeta)),0,
   (dt*(LlambdaL*RlambdaL - LlambdaR*RlambdaR)*Theta*u*Y)/(dx*(LlambdaL - LlambdaR)*(RlambdaL - RlambdaR)*(1 + zeta)),
   (dt*(-(LlambdaL*RlambdaL) + LlambdaR*RlambdaR)*Theta*Y)/(dx*(LlambdaL - LlambdaR)*(RlambdaL - RlambdaR)*(1 + zeta)),0,
   (dt*Theta*(LlambdaL*LlambdaR*(RlambdaL - RlambdaR) + LlambdaL*RlambdaL*(RlambdaR - u) + LlambdaR*RlambdaR*(-RlambdaL + u)))/
    (dx*(LlambdaL - LlambdaR)*(RlambdaL - RlambdaR)*(1 + zeta));
  b += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Viscous (2nd order derivatives)
///////////////////////////////////////////////////////////////////////////////
#if defined(VISCOUS)

  temp << 0,0,0,0,(-8*dt*Pr*Theta*u)/(3.*Power(dx,2)*rho*(1 + zeta)),(8*dt*Pr*Theta)/(3.*Power(dx,2)*rho*(1 + zeta)),0,0,
   (-2*dt*Theta*(3*e*gamma + (-3*gamma + 4*Pr)*rho*Power(u,2)))/(3.*Power(dx,2)*Power(rho,2)*(1 + zeta)),
   (2*dt*(-3*gamma + 4*Pr)*Theta*u)/(3.*Power(dx,2)*rho*(1 + zeta)),(2*dt*gamma*Theta)/(Power(dx,2)*rho*(1 + zeta)),0,
   (-2*dt*Theta*Y)/(Power(dx,2)*Le*rho + Power(dx,2)*Le*rho*zeta),0,0,(2*dt*Theta)/(Power(dx,2)*Le*rho + Power(dx,2)*Le*rho*zeta);
  b += temp;


#endif
///////////////////////////////////////////////////////////////////////////////
// Source
///////////////////////////////////////////////////////////////////////////////
#if defined(SOURCE)

  temp << 0,0,0,0,0,0,0,0,(4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Q*rho*theta*Theta*(e - rho*Power(u,2))*Y)/
    ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),(4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Q*Power(rho,2)*
      theta*Theta*u*Y)/((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),
   (-4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Q*Power(rho,2)*theta*Theta*Y)/
    ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),-((dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Q*Theta)/
      (1 + zeta)),(4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*rho*theta*Theta*(-e + rho*Power(u,2))*Y)/
    ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),(-4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Power(rho,2)*
      theta*Theta*u*Y)/((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),
   (4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Power(rho,2)*theta*Theta*Y)/
    ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),(dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Theta)/(1 + zeta);

  if((gamma-1.0)/ rho * (e-0.5*rho*u*u) > T_ignition){
    b += temp;
  }
#endif

  return b;
}

template <typename grid_type, typename flow_type>
typename grid_type::matrix_type Implicit_HLLE<grid_type, flow_type>::top_matrix() {

  matrix_type c;
  matrix_type temp;
///////////////////////////////////////////////////////////////////////////////
// Needed
///////////////////////////////////////////////////////////////////////////////
  c << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
///////////////////////////////////////////////////////////////////////////////
// Hyperbolic (Euler)
///////////////////////////////////////////////////////////////////////////////
#if defined(HYPERBOLIC)

  temp << (dt*RlambdaL*RlambdaR*Theta)/(dx*(-RlambdaL + RlambdaR)*(1 + zeta)),(dt*RlambdaL*Theta)/(dx*(RlambdaL - RlambdaR)*(1 + zeta)),0,0,
   (dt*(-3 + gamma)*RlambdaL*Theta*Power(up,2))/(2.*dx*(RlambdaL - RlambdaR)*(1 + zeta)),
   -((dt*RlambdaL*Theta*(RlambdaR + (-3 + gamma)*up))/(dx*(RlambdaL - RlambdaR)*(1 + zeta))),
   (dt*(-1 + gamma)*RlambdaL*Theta)/(dx*(RlambdaL - RlambdaR)*(1 + zeta)),0,
   (dt*RlambdaL*Theta*up*(-(ep*gamma) + (-1 + gamma)*rhop*Power(up,2)))/(dx*rhop*(RlambdaL - RlambdaR)*(1 + zeta)),
   (dt*RlambdaL*Theta*(2*ep*gamma - 3*(-1 + gamma)*rhop*Power(up,2)))/(2.*dx*rhop*(RlambdaL - RlambdaR)*(1 + zeta)),
   (dt*RlambdaL*Theta*(-RlambdaR + gamma*up))/(dx*(RlambdaL - RlambdaR)*(1 + zeta)),0,(dt*RlambdaL*Theta*up*Yp)/(dx*(-RlambdaL + RlambdaR)*(1 + zeta)),
   (dt*RlambdaL*Theta*Yp)/(dx*(RlambdaL - RlambdaR)*(1 + zeta)),0,(dt*RlambdaL*Theta*(-RlambdaR + up))/(dx*(RlambdaL - RlambdaR)*(1 + zeta));
  c += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Viscous (2nd order derivatives)
///////////////////////////////////////////////////////////////////////////////
#if defined(VISCOUS)

  temp << 0,0,0,0,(dt*Pr*Theta*(4*rhom*rhop*u - rho*rhop*um + rho*rhom*up))/(3.*Power(dx,2)*rho*rhom*rhop*(1 + zeta)),-(dt*Pr*(rho*(rhom - rhop) + 4*rhom*rhop)*Theta)/(3.*Power(dx,2)*rho*rhom*rhop*(1 + zeta)),0,0,
   (dt*Theta*(-3*e*gamma*Power(rhom,2)*rhop*(Power(rho,2) - 4*Power(rhop,2)) + 3*ep*gamma*Power(rho,2)*Power(rhom,2)*(2*rho + rhop - 2*rhopp) +
        rho*rhop*(3*epp*gamma*rho*Power(rhom,2) + rhop*(-3*em*gamma*rho*rhop + (3*gamma - 4*Pr)*rhom*(rho*rhop*Power(um,2) - rhom*(4*rhop*Power(u,2) + rho*up*(-u + up + upp)))))))/
    (12.*Power(dx,2)*Power(rho,2)*Power(rhom,2)*Power(rhop,3)*(1 + zeta)),(dt*(3*gamma - 4*Pr)*Theta*(4*rhom*rhop*u - rho*rhop*um + rho*rhom*(-u + up + upp)))/(12.*Power(dx,2)*rho*rhom*rhop*(1 + zeta)),
   -(dt*gamma*(Power(rho,2)*rhom + rho*(rhom - rhop)*rhop + 4*rhom*Power(rhop,2) - rho*rhom*rhopp)*Theta)/(4.*Power(dx,2)*rho*rhom*Power(rhop,2)*(1 + zeta)),0,
   (4*dt*rhom*rhop*Theta*Y - dt*rho*rhop*Theta*Ym + dt*rho*rhom*Theta*Yp)/(4*Power(dx,2)*Le*rho*rhom*rhop + 4*Power(dx,2)*Le*rho*rhom*rhop*zeta),0,0,
   -((dt*rho*rhom*Theta - dt*rho*rhop*Theta + 4*dt*rhom*rhop*Theta)/(4*Power(dx,2)*Le*rho*rhom*rhop + 4*Power(dx,2)*Le*rho*rhom*rhop*zeta));
  c += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Source
///////////////////////////////////////////////////////////////////////////////
#if defined(SOURCE)

  temp << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
  c += temp;

#endif

  return c;
}

template <typename grid_type, typename flow_type>
typename grid_type::matrix_type Implicit_HLLE<grid_type, flow_type>::bot_matrix() {

  matrix_type a;
  matrix_type temp;
///////////////////////////////////////////////////////////////////////////////
// Needed
///////////////////////////////////////////////////////////////////////////////
  a << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
///////////////////////////////////////////////////////////////////////////////
// Hyperbolic (Euler)
///////////////////////////////////////////////////////////////////////////////
#if defined(HYPERBOLIC)

  temp << (dt*LlambdaL*LlambdaR*Theta)/(dx*(-LlambdaL + LlambdaR)*(1 + zeta)),(dt*LlambdaR*Theta)/(dx*(LlambdaL - LlambdaR)*(1 + zeta)),0,0,
   (dt*(-3 + gamma)*LlambdaR*Theta*Power(um,2))/(2.*dx*(LlambdaL - LlambdaR)*(1 + zeta)),
   -((dt*LlambdaR*Theta*(LlambdaL + (-3 + gamma)*um))/(dx*(LlambdaL - LlambdaR)*(1 + zeta))),
   (dt*(-1 + gamma)*LlambdaR*Theta)/(dx*(LlambdaL - LlambdaR)*(1 + zeta)),0,
   (dt*LlambdaR*Theta*um*(-(em*gamma) + (-1 + gamma)*rhom*Power(um,2)))/(dx*(LlambdaL - LlambdaR)*rhom*(1 + zeta)),
   (dt*LlambdaR*Theta*(2*em*gamma - 3*(-1 + gamma)*rhom*Power(um,2)))/(2.*dx*(LlambdaL - LlambdaR)*rhom*(1 + zeta)),
   (dt*LlambdaR*Theta*(-LlambdaL + gamma*um))/(dx*(LlambdaL - LlambdaR)*(1 + zeta)),0,(dt*LlambdaR*Theta*um*Ym)/(dx*(-LlambdaL + LlambdaR)*(1 + zeta)),
   (dt*LlambdaR*Theta*Ym)/(dx*(LlambdaL - LlambdaR)*(1 + zeta)),0,(dt*LlambdaR*Theta*(-LlambdaL + um))/(dx*(LlambdaL - LlambdaR)*(1 + zeta));
  a += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Viscous (2nd order derivatives)
///////////////////////////////////////////////////////////////////////////////
#if defined(VISCOUS)

  temp << 0,0,0,0,(dt*Pr*Theta*(4*rhom*rhop*u + rho*rhop*um - rho*rhom*up))/(3.*Power(dx,2)*rho*rhom*rhop*(1 + zeta)),
   (dt*Pr*(rho*rhom - (rho + 4*rhom)*rhop)*Theta)/(3.*Power(dx,2)*rho*rhom*rhop*(1 + zeta)),0,0,
   (dt*Theta*(-3*ep*gamma*Power(rho,2)*Power(rhom,3) + rhop*(-3*e*gamma*rhom*(Power(rho,2) - 4*Power(rhom,2))*rhop +
           3*em*gamma*Power(rho,2)*(2*rho + rhom - 2*rhomm)*rhop +
           rho*rhom*(3*emm*gamma*rho*rhop - (3*gamma - 4*Pr)*rhom*(4*rhom*rhop*Power(u,2) + rho*rhop*um*(-u + um + umm) - rho*rhom*Power(up,2))))))/
    (12.*Power(dx,2)*Power(rho,2)*Power(rhom,3)*Power(rhop,2)*(1 + zeta)),
   (dt*(3*gamma - 4*Pr)*Theta*(4*rhom*rhop*u + rho*rhop*(-u + um + umm) - rho*rhom*up))/(12.*Power(dx,2)*rho*rhom*rhop*(1 + zeta)),
   (dt*gamma*(rho*Power(rhom,2) - (Power(rho,2) + 4*Power(rhom,2) + rho*(rhom - rhomm))*rhop)*Theta)/(4.*Power(dx,2)*rho*Power(rhom,2)*rhop*(1 + zeta)),
   0,(4*dt*rhom*rhop*Theta*Y + dt*rho*rhop*Theta*Ym - dt*rho*rhom*Theta*Yp)/(4*Power(dx,2)*Le*rho*rhom*rhop + 4*Power(dx,2)*Le*rho*rhom*rhop*zeta),0,0,
   (dt*rho*rhom*Theta - dt*rho*rhop*Theta - 4*dt*rhom*rhop*Theta)/(4*Power(dx,2)*Le*rho*rhom*rhop + 4*Power(dx,2)*Le*rho*rhom*rhop*zeta);
  a += temp;


#endif
///////////////////////////////////////////////////////////////////////////////
// Source
///////////////////////////////////////////////////////////////////////////////
#if defined(SOURCE)

  temp << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
  a += temp;

#endif
  return a;
}

template <typename grid_type, typename flow_type>
typename grid_type::global_solution_vector_type::value_type
Implicit_HLLE<grid_type, flow_type>::rhs_matrix() {

  HLLE<global_solution_vector_type> hyperbolic_flux;
  solution_vector_type rhs;
  solution_vector_type temp;
///////////////////////////////////////////////////////////////////////////////
// Needed
///////////////////////////////////////////////////////////////////////////////
  rhs << (dUm1*zeta)/(1 + zeta),(dUm2*zeta)/(1 + zeta),(dUm3*zeta)/(1 + zeta),(dUm4*zeta)/(1 + zeta);
///////////////////////////////////////////////////////////////////////////////
// Hyperbolic (Euler)
///////////////////////////////////////////////////////////////////////////////
#if defined(HYPERBOLIC)

  // temp << (dt*((LlambdaL*LlambdaR*(rho - rhom))/(-LlambdaL + LlambdaR) - ((-rho + rhop)*RlambdaL*RlambdaR)/(-RlambdaL + RlambdaR) +
  //       (-(LlambdaL*rho*u) + LlambdaR*rhom*um)/(-LlambdaL + LlambdaR) - (rho*RlambdaR*u - rhop*RlambdaL*up)/(-RlambdaL + RlambdaR)))/(dx*(1 + zeta)),
  //  (dt*((LlambdaL*LlambdaR*(rho*u - rhom*um))/(-LlambdaL + LlambdaR) +
  //       (-(LlambdaL*(rho*Power(u,2) + (-1 + gamma)*(e - (rho*Power(u,2))/2.))) +
  //          LlambdaR*(rhom*Power(um,2) + (-1 + gamma)*(em - (rhom*Power(um,2))/2.)))/(-LlambdaL + LlambdaR) -
  //       (RlambdaL*RlambdaR*(-(rho*u) + rhop*up))/(-RlambdaL + RlambdaR) -
  //       (RlambdaR*(rho*Power(u,2) + (-1 + gamma)*(e - (rho*Power(u,2))/2.)) - RlambdaL*(rhop*Power(up,2) + (-1 + gamma)*(ep - (rhop*Power(up,2))/2.)))/
  //        (-RlambdaL + RlambdaR)))/(dx*(1 + zeta)),(dt*(((e - em)*LlambdaL*LlambdaR)/(-LlambdaL + LlambdaR) -
  //       ((-e + ep)*RlambdaL*RlambdaR)/(-RlambdaL + RlambdaR) +
  //       (-(LlambdaL*u*(e + (-1 + gamma)*(e - (rho*Power(u,2))/2.))) + LlambdaR*um*(em + (-1 + gamma)*(em - (rhom*Power(um,2))/2.)))/
  //        (-LlambdaL + LlambdaR) - (RlambdaR*u*(e + (-1 + gamma)*(e - (rho*Power(u,2))/2.)) -
  //          RlambdaL*up*(ep + (-1 + gamma)*(ep - (rhop*Power(up,2))/2.)))/(-RlambdaL + RlambdaR)))/(dx*(1 + zeta)),
  //  (dt*((LlambdaL*LlambdaR*(rho*Y - rhom*Ym))/(-LlambdaL + LlambdaR) + (-(LlambdaL*rho*u*Y) + LlambdaR*rhom*um*Ym)/(-LlambdaL + LlambdaR) -
  //       (RlambdaL*RlambdaR*(-(rho*Y) + rhop*Yp))/(-RlambdaL + RlambdaR) - (rho*RlambdaR*u*Y - rhop*RlambdaL*up*Yp)/(-RlambdaL + RlambdaR)))/
  //   (dx*(1 + zeta));
  // rhs += temp;


  auto var_vec_mm = Variable_Vector_Isolator<grid_type>(solution_vector_mm, gamma);
  auto var_vec_m = Variable_Vector_Isolator<grid_type>(solution_vector_m, gamma);
  auto var_vec = Variable_Vector_Isolator<grid_type>(solution_vector, gamma);
  auto var_vec_p = Variable_Vector_Isolator<grid_type>(solution_vector_p, gamma);
  auto var_vec_pp = Variable_Vector_Isolator<grid_type>(solution_vector_pp, gamma);
  solution_vector_type phi_m = limiter(var_vec_mm.w(),var_vec_m.w(), var_vec.w(), dx);
  solution_vector_type phi = limiter(var_vec_m.w(),var_vec.w(), var_vec_p.w(), dx);
  solution_vector_type phi_p = limiter(var_vec.w(),var_vec_p.w(), var_vec_pp.w(), dx);
  solution_vector_type LUl = var_vec_m.w().array() + phi_m.array() *dx / 2.0;
  solution_vector_type LUr =   var_vec.w().array() - phi.array() * dx / 2.0;
  solution_vector_type RUl =   var_vec.w().array() + phi.array() *dx / 2.0;
  solution_vector_type RUr = var_vec_p.w().array() - phi_p.array() * dx / 2.0;

   rhs += 1.0 / (1.0 + zeta) * (hyperbolic_flux.flux(LUl, LUr, gamma) - hyperbolic_flux.flux(RUl, RUr, gamma)) /dx*dt;

#endif
///////////////////////////////////////////////////////////////////////////////
// Viscous (2nd order derivatives)
///////////////////////////////////////////////////////////////////////////////
#if defined(VISCOUS)

  temp << 0,(4*dt*Pr*(-2*u + um + up))/(3.*Power(dx,2)*(1 + zeta)),
   (dt*((Pr*Power(-um + up,2))/(3.*Power(dx,2)) + (4*Pr*u*(-2*u + um + up))/(3.*Power(dx,2)) +
        (gamma*(((-1 + gamma)*Power(-rhom + rhop,2)*(e - (rho*Power(u,2))/2.))/(2.*Power(dx,2)*Power(rho,3)) -
             ((-1 + gamma)*(-2*rho + rhom + rhop)*(e - (rho*Power(u,2))/2.))/(Power(dx,2)*Power(rho,2)) -
             ((-1 + gamma)*(-rhom + rhop)*((-em + ep)/(2.*dx) - ((-rhom + rhop)*Power(u,2))/(4.*dx) - (rho*u*(-um + up))/(2.*dx)))/(dx*Power(rho,2)) +
             ((-1 + gamma)*((-2*e + em + ep)/Power(dx,2) - ((-2*rho + rhom + rhop)*Power(u,2))/(2.*Power(dx,2)) -
                  ((-rhom + rhop)*u*(-um + up))/(2.*Power(dx,2)) - (rho*Power(-um + up,2))/(4.*Power(dx,2)) - (rho*u*(-2*u + um + up))/Power(dx,2)))/rho
             ))/(-1 + gamma)))/(1 + zeta),(dt*(-2*Y + Ym + Yp))/(Power(dx,2)*Le*(1 + zeta));
  rhs += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Source
///////////////////////////////////////////////////////////////////////////////
#if defined(SOURCE)

  temp << 0.,0.,(dt*lambda*Q*rho*Y)/(Power(E,(rho*theta)/((-1 + gamma)*(e - (rho*Power(u,2))/2.)))*(1 + zeta)),
   -((dt*lambda*rho*Y)/(Power(E,(rho*theta)/((-1 + gamma)*(e - (rho*Power(u,2))/2.)))*(1 + zeta)));

  if((gamma-1.0)/ rho * (e-0.5*rho*u*u) > T_ignition){
    rhs += temp;
  }

#endif

  return rhs;
}

#endif //#ifndef IMPLICIT_HLLE_H
