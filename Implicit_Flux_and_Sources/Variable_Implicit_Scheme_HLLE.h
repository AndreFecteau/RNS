#ifndef CREATE_IMPLICIT_MATRIX_VECTOR_HLLE_H
#define CREATE_IMPLICIT_MATRIX_VECTOR_HLLE_H
//
#include "../Explicit_Flux_and_Sources/HLLE.h"
#include "../Explicit_Flux_and_Sources/Centered_Difference.h"
#include "../Explicit_Flux_and_Sources/Sources.h"
// #include "Solver/Explicit_Marching.h"

#define HYPERBOLIC
#define VISCOUS
#define SOURCE
typedef Eigen::Matrix<double, 4, 1> Vector_type;

template <typename T>
T Power(T num, int expo) {
  return pow(num, expo);
}

double Power(int num, double expo) {
  return exp(expo);
}

double u_rhoavg(double rho_left, double rho_right, double u_left, double u_right) {
  return (sqrt(rho_left) * u_left + sqrt(rho_right) * u_right) /
         (sqrt(rho_left) + sqrt(rho_right));
}

double rho_rhoavg(double rho_left, double rho_right) {
  return sqrt(rho_left * rho_right);
}

double p_rhoavg(double rho_left, double rho_right, double u_left, double u_right,
                double p_left, double p_right, double gamma) {

  double h_left = gamma * p_left / (rho_left * (gamma - 1.0)) + u_left * u_left * 0.5;
  double h_right = gamma * p_right / (rho_right * (gamma - 1.0)) + u_right * u_right * 0.5;

  double h_ravg = (sqrt(rho_left) * h_left + sqrt(rho_right) * h_right) /
                  (sqrt(rho_left) + sqrt(rho_right));

  double u_ravg = u_rhoavg(rho_left, rho_right, u_left, u_right);

  double rho_ravg = rho_rhoavg(rho_left, rho_right);

  return (h_ravg - u_ravg*u_ravg * 0.5) * (gamma - 1.0) / gamma * rho_ravg;
}

Vector_type limiter(const Vector_type &Ul,
                    const Vector_type &U,
                    const Vector_type &Ur,
                    double dx) {
  Vector_type a = (U - Ul) / dx;
  Vector_type b = (Ur - U) / dx;
  double epsilon = 1.0e-6;
  Vector_type phi =  a.array() * b.array() * (a.array() + b.array()) / (a.array() * a.array() + b.array() * b.array() + epsilon);
  for (int i = 0; i < 3; ++i) {
    if (a[i] / b[i] <= 0.0 || b[i] == 0) {
      phi[i] = 0.0;
    }
  }
  return phi;

}

template <typename solution_vector_type, typename Matrix_type>
Matrix_type create_mid_band_matrix(const solution_vector_type& solution_vector_mm,
                                       const solution_vector_type solution_vector_m,
                                       const solution_vector_type solution_vector,
                                       const solution_vector_type solution_vector_p,
                                       const solution_vector_type solution_vector_pp,
                                       const double& gamma, const double& Pr,
                                       const double& Le, const double& Q, const double& Lambda,
                                       const double& theta, const double& dx, const double& dt,
                                       const double& zeta, const double& Theta) {
  double rhomm = solution_vector_mm[0];
  double rhom = solution_vector_m[0];
  double rho = solution_vector[0];
  double rhop = solution_vector_p[0];
  double rhopp = solution_vector_pp[0];
  double umm = solution_vector_mm[1]/solution_vector_mm[0];
  double um = solution_vector_m[1]/solution_vector_m[0];
  double u = solution_vector[1]/solution_vector[0];
  double up = solution_vector_p[1]/solution_vector_p[0];
  double upp = solution_vector_pp[1]/solution_vector_pp[0];
  double emm = solution_vector_mm[2];
  double em = solution_vector_m[2];
  double e = solution_vector[2];
  double ep = solution_vector_p[2];
  double epp = solution_vector_pp[2];
  double Ymm = solution_vector_mm[3]/solution_vector_mm[0];
  double Ym = solution_vector_m[3]/solution_vector_m[0];
  double Y = solution_vector[3]/solution_vector[0];
  double Yp = solution_vector_p[3]/solution_vector_p[0];
  double Ypp = solution_vector_pp[3]/solution_vector_pp[0];
  double pm = (solution_vector_m[2] - solution_vector_m[1]* solution_vector_m[1] /
             (2.0 * solution_vector_m[0])) * (gamma - 1);
  double p = (solution_vector[2] - solution_vector[1]* solution_vector[1] /
             (2.0 * solution_vector[0])) * (gamma - 1);
  double pp = (solution_vector_p[2] - solution_vector_p[1]* solution_vector_p[1] /
             (2.0 * solution_vector_p[0])) * (gamma - 1);
  double LlambdaR;
  double LlambdaL;
  double RlambdaL;
  double RlambdaR;
  int E = 0;

  {
  double lambda_1 = u + sqrt(gamma*p/rho);
  double lambda_2 = u_rhoavg(rhom, rho, um, u) + sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  LlambdaR = std::max(lambda_1,lambda_2);
  }
  {
  double lambda_1 = um - sqrt(gamma*pm/rhom);
  double lambda_2 = u_rhoavg(rhom, rho, um, u) - sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  LlambdaL = std::min(lambda_1,lambda_2);
  }
  {
  double lambda_1 = up + sqrt(gamma*pp/rhop);
  double lambda_2 = u_rhoavg(rho, rhop, u, up) + sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  RlambdaR = std::max(lambda_1,lambda_2);
  }
  {
  double lambda_1 = u - sqrt(gamma*p/rho);
  double lambda_2 = u_rhoavg(rho, rhop, u, up) - sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  RlambdaL = std::min(lambda_1,lambda_2);
  }

  Matrix_type b;
  Matrix_type temp;
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

  temp << 0,0,0,0,0,0,0,0,(4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*Lambda*Q*rho*theta*Theta*(e - rho*Power(u,2))*Y)/
    ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),
   (4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*Lambda*Q*Power(rho,2)*theta*Theta*u*Y)/
    ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),
   (-4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*Lambda*Q*Power(rho,2)*theta*Theta*Y)/
    ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),
   -((dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*Lambda*Q*Theta)/(1 + zeta)),
   (4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*Lambda*rho*theta*Theta*(-e + rho*Power(u,2))*Y)/
    ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),
   (-4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*Lambda*Power(rho,2)*theta*Theta*u*Y)/
    ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),
   (4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*Lambda*Power(rho,2)*theta*Theta*Y)/
    ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),(dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*Lambda*Theta)/(1 + zeta);
  b += temp;

#endif

  return b;
}

template <typename solution_vector_type, typename Matrix_type>
Matrix_type create_top_band_matrix(const solution_vector_type& solution_vector_mm,
                                       const solution_vector_type solution_vector_m,
                                       const solution_vector_type solution_vector,
                                       const solution_vector_type solution_vector_p,
                                       const solution_vector_type solution_vector_pp,
                                       const double& gamma, const double& Pr,
                                       const double& Le, const double& Q, const double& Lambda,
                                       const double& theta, const double& dx, const double& dt,
                                       const double& zeta, const double& Theta) {
  double rhomm = solution_vector_mm[0];
  double rhom = solution_vector_m[0];
  double rho = solution_vector[0];
  double rhop = solution_vector_p[0];
  double rhopp = solution_vector_pp[0];
  double umm = solution_vector_mm[1]/solution_vector_mm[0];
  double um = solution_vector_m[1]/solution_vector_m[0];
  double u = solution_vector[1]/solution_vector[0];
  double up = solution_vector_p[1]/solution_vector_p[0];
  double upp = solution_vector_pp[1]/solution_vector_pp[0];
  double emm = solution_vector_mm[2];
  double em = solution_vector_m[2];
  double e = solution_vector[2];
  double ep = solution_vector_p[2];
  double epp = solution_vector_pp[2];
  double Ymm = solution_vector_mm[3]/solution_vector_mm[0];
  double Ym = solution_vector_m[3]/solution_vector_m[0];
  double Y = solution_vector[3]/solution_vector[0];
  double Yp = solution_vector_p[3]/solution_vector_p[0];
  double Ypp = solution_vector_pp[3]/solution_vector_pp[0];
  double pm = (solution_vector_m[2] - solution_vector_m[1]* solution_vector_m[1] /
             (2.0 * solution_vector_m[0])) * (gamma - 1);
  double p = (solution_vector[2] - solution_vector[1]* solution_vector[1] /
             (2.0 * solution_vector[0])) * (gamma - 1);
  double pp = (solution_vector_p[2] - solution_vector_p[1]* solution_vector_p[1] /
             (2.0 * solution_vector_p[0])) * (gamma - 1);
  double LlambdaR;
  double LlambdaL;
  double RlambdaL;
  double RlambdaR;
  int E = 0;

  {
  double lambda_1 = u + sqrt(gamma*p/rho);
  double lambda_2 = u_rhoavg(rhom, rho, um, u) + sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  LlambdaR = std::max(lambda_1,lambda_2);
  }
  {
  double lambda_1 = um - sqrt(gamma*pm/rhom);
  double lambda_2 = u_rhoavg(rhom, rho, um, u) - sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  LlambdaL = std::min(lambda_1,lambda_2);
  }
  {
  double lambda_1 = up + sqrt(gamma*pp/rhop);
  double lambda_2 = u_rhoavg(rho, rhop, u, up) + sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  RlambdaR = std::max(lambda_1,lambda_2);
  }
  {
  double lambda_1 = u - sqrt(gamma*p/rho);
  double lambda_2 = u_rhoavg(rho, rhop, u, up) - sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  RlambdaL = std::min(lambda_1,lambda_2);
  }
  Matrix_type c;
  Matrix_type temp;
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

template <typename solution_vector_type, typename Matrix_type>
Matrix_type create_bot_band_matrix(const solution_vector_type& solution_vector_mm,
                                       const solution_vector_type solution_vector_m,
                                       const solution_vector_type solution_vector,
                                       const solution_vector_type solution_vector_p,
                                       const solution_vector_type solution_vector_pp,
                                       const double& gamma, const double& Pr,
                                       const double& Le, const double& Q, const double& Lambda,
                                       const double& theta, const double& dx, const double& dt,
                                       const double& zeta, const double& Theta) {
  double rhomm = solution_vector_mm[0];
  double rhom = solution_vector_m[0];
  double rho = solution_vector[0];
  double rhop = solution_vector_p[0];
  double rhopp = solution_vector_pp[0];
  double umm = solution_vector_mm[1]/solution_vector_mm[0];
  double um = solution_vector_m[1]/solution_vector_m[0];
  double u = solution_vector[1]/solution_vector[0];
  double up = solution_vector_p[1]/solution_vector_p[0];
  double upp = solution_vector_pp[1]/solution_vector_pp[0];
  double emm = solution_vector_mm[2];
  double em = solution_vector_m[2];
  double e = solution_vector[2];
  double ep = solution_vector_p[2];
  double epp = solution_vector_pp[2];
  double Ymm = solution_vector_mm[3]/solution_vector_mm[0];
  double Ym = solution_vector_m[3]/solution_vector_m[0];
  double Y = solution_vector[3]/solution_vector[0];
  double Yp = solution_vector_p[3]/solution_vector_p[0];
  double Ypp = solution_vector_pp[3]/solution_vector_pp[0];
  double pm = (solution_vector_m[2] - solution_vector_m[1]* solution_vector_m[1] /
             (2.0 * solution_vector_m[0])) * (gamma - 1);
  double p = (solution_vector[2] - solution_vector[1]* solution_vector[1] /
             (2.0 * solution_vector[0])) * (gamma - 1);
  double pp = (solution_vector_p[2] - solution_vector_p[1]* solution_vector_p[1] /
             (2.0 * solution_vector_p[0])) * (gamma - 1);
  double LlambdaR;
  double LlambdaL;
  double RlambdaL;
  double RlambdaR;
  int E = 0;

  {
  double lambda_1 = u + sqrt(gamma*p/rho);
  double lambda_2 = u_rhoavg(rhom, rho, um, u) + sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  LlambdaR = std::max(lambda_1,lambda_2);
  }
  {
  double lambda_1 = um - sqrt(gamma*pm/rhom);
  double lambda_2 = u_rhoavg(rhom, rho, um, u) - sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  LlambdaL = std::min(lambda_1,lambda_2);
  }
  {
  double lambda_1 = up + sqrt(gamma*pp/rhop);
  double lambda_2 = u_rhoavg(rho, rhop, u, up) + sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  RlambdaR = std::max(lambda_1,lambda_2);
  }
  {
  double lambda_1 = u - sqrt(gamma*p/rho);
  double lambda_2 = u_rhoavg(rho, rhop, u, up) - sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  RlambdaL = std::min(lambda_1,lambda_2);
  }

  Matrix_type a;
  Matrix_type temp;
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

template <typename solution_vector_type, typename Matrix_type>
solution_vector_type create_rhs_vector(const solution_vector_type& solution_vector_mm,
                                       const solution_vector_type solution_vector_m,
                                       const solution_vector_type solution_vector,
                                       const solution_vector_type solution_vector_p,
                                       const solution_vector_type solution_vector_pp,
                                       const double& gamma, const double& Pr,
                                       const double& Le, const double& Q, const double& Lambda,
                                       const double& theta, const double& dx, const double& dt,
                                       const double& zeta, const double& Theta,
                                       const solution_vector_type& DeltaUm) {

  double rhomm = solution_vector_mm[0];
  double rhom = solution_vector_m[0];
  double rho = solution_vector[0];
  double rhop = solution_vector_p[0];
  double rhopp = solution_vector_pp[0];
  double umm = solution_vector_mm[1]/solution_vector_mm[0];
  double um = solution_vector_m[1]/solution_vector_m[0];
  double u = solution_vector[1]/solution_vector[0];
  double up = solution_vector_p[1]/solution_vector_p[0];
  double upp = solution_vector_pp[1]/solution_vector_pp[0];
  double emm = solution_vector_mm[2];
  double em = solution_vector_m[2];
  double e = solution_vector[2];
  double ep = solution_vector_p[2];
  double epp = solution_vector_pp[2];
  double Ymm = solution_vector_mm[3]/solution_vector_mm[0];
  double Ym = solution_vector_m[3]/solution_vector_m[0];
  double Y = solution_vector[3]/solution_vector[0];
  double Yp = solution_vector_p[3]/solution_vector_p[0];
  double Ypp = solution_vector_pp[3]/solution_vector_pp[0];
  double pm = (solution_vector_m[2] - solution_vector_m[1]* solution_vector_m[1] /
             (2.0 * solution_vector_m[0])) * (gamma - 1);
  double p = (solution_vector[2] - solution_vector[1]* solution_vector[1] /
             (2.0 * solution_vector[0])) * (gamma - 1);
  double pp = (solution_vector_p[2] - solution_vector_p[1]* solution_vector_p[1] /
             (2.0 * solution_vector_p[0])) * (gamma - 1);
  double dUm1 = DeltaUm[0];
  double dUm2 = DeltaUm[1];
  double dUm3 = DeltaUm[2];
  double dUm4 = DeltaUm[3];
  double LlambdaR;
  double LlambdaL;
  double RlambdaL;
  double RlambdaR;
  int E = 0;

  {
  double lambda_1 = u + sqrt(gamma*p/rho);
  double lambda_2 = u_rhoavg(rhom, rho, um, u) + sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  LlambdaR = std::max(lambda_1,lambda_2);
  }
  {
  double lambda_1 = um - sqrt(gamma*pm/rhom);
  double lambda_2 = u_rhoavg(rhom, rho, um, u) - sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  LlambdaL = std::min(lambda_1,lambda_2);
  }
  {
  double lambda_1 = up + sqrt(gamma*pp/rhop);
  double lambda_2 = u_rhoavg(rho, rhop, u, up) + sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  RlambdaR = std::max(lambda_1,lambda_2);
  }
  {
  double lambda_1 = u - sqrt(gamma*p/rho);
  double lambda_2 = u_rhoavg(rho, rhop, u, up) - sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  RlambdaL = std::min(lambda_1,lambda_2);
  }

  HLLE<std::vector<solution_vector_type>> hyperbolic_flux;
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


  auto var_vec_mm = Variable_Vector_Isolator<solution_vector_type>(solution_vector_mm, gamma);
  auto var_vec_m = Variable_Vector_Isolator<solution_vector_type>(solution_vector_m, gamma);
  auto var_vec = Variable_Vector_Isolator<solution_vector_type>(solution_vector, gamma);
  auto var_vec_p = Variable_Vector_Isolator<solution_vector_type>(solution_vector_p, gamma);
  auto var_vec_pp = Variable_Vector_Isolator<solution_vector_type>(solution_vector_pp, gamma);
  solution_vector_type phi_m = limiter(var_vec_mm.w(),var_vec_m.w(), var_vec.w(), dx);
  solution_vector_type phi = limiter(var_vec_m.w(),var_vec.w(), var_vec_p.w(), dx);
  solution_vector_type phi_p = limiter(var_vec.w(),var_vec_p.w(), var_vec_pp.w(), dx);
  solution_vector_type LUl = var_vec_m.w().array() + phi_m.array() *dx / 2.0;
  solution_vector_type LUr =   var_vec.w().array() - phi.array() * dx / 2.0;
  solution_vector_type RUl =   var_vec.w().array() + phi.array() *dx / 2.0;
  solution_vector_type RUr = var_vec_p.w().array() - phi_p.array() * dx / 2.0;
  // auto var_vec_l = Variable_Vector_Isolator<solution_vector_type>(solution_vector_m, gamma);
  // auto var_vec = Variable_Vector_Isolator<solution_vector_type>(solution_vector, gamma);
  // auto var_vec_r = Variable_Vector_Isolator<solution_vector_type>(solution_vector_p, gamma);

   rhs += (hyperbolic_flux.flux(LUl, LUr, gamma) - hyperbolic_flux.flux(RUl, RUr, gamma)) /dx*dt;

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

  temp << 0.,0.,(dt*Lambda*Q*rho*Y)/(Power(E,(rho*theta)/((-1 + gamma)*(e - (rho*Power(u,2))/2.)))*(1 + zeta)),
   -((dt*Lambda*rho*Y)/(Power(E,(rho*theta)/((-1 + gamma)*(e - (rho*Power(u,2))/2.)))*(1 + zeta)));

  rhs += temp;

#endif

  return rhs;
}

#endif //#ifndef CREATE_IMPLICIT_MATRIX_VECTOR_H
