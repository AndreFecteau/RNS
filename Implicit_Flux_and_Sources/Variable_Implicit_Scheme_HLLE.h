#ifndef CREATE_IMPLICIT_MATRIX_VECTOR_HLLE_H
#define CREATE_IMPLICIT_MATRIX_VECTOR_HLLE_H
//
#include "../Explicit_Flux_and_Sources/HLLE.h"
#include "../Explicit_Flux_and_Sources/Centered_Difference.h"
#include "../Explicit_Flux_and_Sources/Sources.h"

#define HYPERBOLIC
#define VISCOUS
#define SOURCE

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
  double Llambda;
  double Llambdam;
  double Rlambda;
  double Rlambdap;
  int E = 0;

  {
  double lambda_1 = u + sqrt(gamma*p/rho);
  double lambda_2 = u_rhoavg(rhom, rho, um, u) + sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  Llambda = std::max(lambda_1,lambda_2);
  }
  {
  double lambda_1 = um - sqrt(gamma*pm/rhom);
  double lambda_2 = u_rhoavg(rhom, rho, um, u) - sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  Llambdam = std::min(lambda_1,lambda_2);
  }
  {
  double lambda_1 = up + sqrt(gamma*pp/rhop);
  double lambda_2 = u_rhoavg(rho, rhop, u, up) + sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  Rlambdap = std::max(lambda_1,lambda_2);
  }
  {
  double lambda_1 = u - sqrt(gamma*p/rho);
  double lambda_2 = u_rhoavg(rho, rhop, u, up) - sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  Rlambda = std::min(lambda_1,lambda_2);
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

  temp << 0,(dt*(Llambdam/(Llambda - Llambdam) + Rlambdap/(-Rlambda + Rlambdap))*Theta)/(dx*(1 + zeta)),0,0,
   (dt*(-3 + gamma)*(Llambdam*Rlambda - Llambda*Rlambdap)*Theta*Power(u,2))/(2.*dx*(Llambda - Llambdam)*(Rlambda - Rlambdap)*(1 + zeta)),
   (dt*(-3 + gamma)*(-(Llambdam*Rlambda) + Llambda*Rlambdap)*Theta*u)/(dx*(Llambda - Llambdam)*(Rlambda - Rlambdap)*(1 + zeta)),
   (dt*(-1 + gamma)*(Llambdam/(Llambda - Llambdam) + Rlambdap/(-Rlambda + Rlambdap))*Theta)/(dx*(1 + zeta)),0,
   -((dt*(Llambdam*Rlambda - Llambda*Rlambdap)*Theta*(e*gamma*u - (-1 + gamma)*rho*Power(u,3)))/
      (dx*(Llambda - Llambdam)*rho*(Rlambda - Rlambdap)*(1 + zeta))),
   (dt*(Llambdam*Rlambda - Llambda*Rlambdap)*Theta*(2*e*gamma - 3*(-1 + gamma)*rho*Power(u,2)))/
    (2.*dx*(Llambda - Llambdam)*rho*(Rlambda - Rlambdap)*(1 + zeta)),
   (dt*gamma*(Llambdam/(Llambda - Llambdam) + Rlambdap/(-Rlambda + Rlambdap))*Theta*u)/(dx*(1 + zeta)),0,
   (dt*(-(Llambdam*Rlambda) + Llambda*Rlambdap)*Theta*u*Y)/(dx*(Llambda - Llambdam)*(Rlambda - Rlambdap)*(1 + zeta)),
   (dt*Theta*((Llambdam*Y)/(Llambda - Llambdam) + (Rlambdap*Y)/(-Rlambda + Rlambdap)))/(dx*(1 + zeta)),0,
   (dt*Theta*((Llambdam*u)/(Llambda - Llambdam) + (Rlambdap*u)/(-Rlambda + Rlambdap)))/(dx*(1 + zeta));
  b += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Viscous (2nd order derivatives)
///////////////////////////////////////////////////////////////////////////////
#if defined(VISCOUS)

  temp << 0,0,0,0,(2*dt*Pr*Theta*(Power(rhom - rhop,2)*u + 2*Power(rho,2)*(-2*u + um + up) - rho*(rhom*(2*u + um - up) + rhop*(2*u - um + up))))/
    (3.*Power(dx,2)*Power(rho,3)*(1 + zeta)),(2*dt*Pr*(-Power(rhom - rhop,2) + 2*rho*(rhom + rhop))*Theta)/(3.*Power(dx,2)*Power(rho,3)*(1 + zeta)),0,0,
   (dt*Theta*(-12*e*gamma*rho*rhom + 9*e*gamma*Power(rhom,2) + 6*ep*gamma*rho*(rho + rhom - rhop) - 12*e*gamma*rho*rhop - 18*e*gamma*rhom*rhop +
        9*e*gamma*Power(rhop,2) + 6*em*gamma*rho*(rho - rhom + rhop) + 24*gamma*Power(rho,3)*Power(u,2) - 32*Pr*Power(rho,3)*Power(u,2) +
        6*gamma*Power(rho,2)*rhom*Power(u,2) - 8*Pr*Power(rho,2)*rhom*Power(u,2) - 3*gamma*rho*Power(rhom,2)*Power(u,2) +
        4*Pr*rho*Power(rhom,2)*Power(u,2) + 6*gamma*Power(rho,2)*rhop*Power(u,2) - 8*Pr*Power(rho,2)*rhop*Power(u,2) +
        6*gamma*rho*rhom*rhop*Power(u,2) - 8*Pr*rho*rhom*rhop*Power(u,2) - 3*gamma*rho*Power(rhop,2)*Power(u,2) + 4*Pr*rho*Power(rhop,2)*Power(u,2) -
        12*gamma*Power(rho,3)*u*um + 16*Pr*Power(rho,3)*u*um + 6*gamma*Power(rho,2)*rhom*u*um - 8*Pr*Power(rho,2)*rhom*u*um -
        6*gamma*Power(rho,2)*rhop*u*um + 8*Pr*Power(rho,2)*rhop*u*um - 3*gamma*Power(rho,3)*Power(um,2) + 4*Pr*Power(rho,3)*Power(um,2) -
        2*(3*gamma - 4*Pr)*Power(rho,2)*((2*rho + rhom - rhop)*u - rho*um)*up + (-3*gamma + 4*Pr)*Power(rho,3)*Power(up,2)))/
    (6.*Power(dx,2)*Power(rho,4)*(1 + zeta)),-(dt*(3*gamma - 4*Pr)*Theta*
       (-(Power(rhom - rhop,2)*u) + rho*(rhom*(2*u + um - up) + rhop*(2*u - um + up)) + Power(rho,2)*(4*u - 2*(um + up))))/
    (6.*Power(dx,2)*Power(rho,3)*(1 + zeta)),(dt*gamma*(-Power(rhom - rhop,2) + 2*rho*(rhom + rhop))*Theta)/(2.*Power(dx,2)*Power(rho,3)*(1 + zeta)),0,
   (dt*Theta*(Power(rhom - rhop,2)*Y + 2*Power(rho,2)*(-2*Y + Ym + Yp) - rho*(rhom*(2*Y + Ym - Yp) + rhop*(2*Y - Ym + Yp))))/
    (2.*Power(dx,2)*Le*Power(rho,3)*(1 + zeta)),0,0,(dt*(-Power(rhom - rhop,2) + 2*rho*(rhom + rhop))*Theta)/(2.*Power(dx,2)*Le*Power(rho,3)*(1 + zeta));
  b += temp;


#endif
///////////////////////////////////////////////////////////////////////////////
// Source
///////////////////////////////////////////////////////////////////////////////
#if defined(SOURCE)

  temp << 0,0,0,0,0,0,0,0,(4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*Lambda*Q*rho*theta*Theta*(e - rho*Power(u,2))*Y)/
    ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),(4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*Lambda*Q*Power(rho,2)*
      theta*Theta*u*Y)/((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),
   (-4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*Lambda*Q*Power(rho,2)*theta*Theta*Y)/
    ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),-((dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*Lambda*Q*Theta)/
      (1 + zeta)),(4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*Lambda*rho*theta*Theta*(-e + rho*Power(u,2))*Y)/
    ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),(-4*dt*Power(E,(2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*Lambda*Power(rho,2)*
      theta*Theta*u*Y)/((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),
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
  double Llambda;
  double Llambdam;
  double Rlambda;
  double Rlambdap;
  int E = 0;

  {
  double lambda_1 = u + sqrt(gamma*p/rho);
  double lambda_2 = u_rhoavg(rhom, rho, um, u) + sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  Llambda = std::max(lambda_1,lambda_2);
  }
  {
  double lambda_1 = um - sqrt(gamma*pm/rhom);
  double lambda_2 = u_rhoavg(rhom, rho, um, u) - sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  Llambdam = std::min(lambda_1,lambda_2);
  }
  {
  double lambda_1 = up + sqrt(gamma*pp/rhop);
  double lambda_2 = u_rhoavg(rho, rhop, u, up) + sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  Rlambdap = std::max(lambda_1,lambda_2);
  }
  {
  double lambda_1 = u - sqrt(gamma*p/rho);
  double lambda_2 = u_rhoavg(rho, rhop, u, up) - sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  Rlambda = std::min(lambda_1,lambda_2);
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

  temp << 0,(dt*Rlambda*Theta)/(dx*(Rlambda - Rlambdap)*(1 + zeta)),0,0,(dt*(-3 + gamma)*Rlambda*Theta*Power(up,2))/(2.*dx*(Rlambda - Rlambdap)*(1 + zeta)),
   -((dt*(-3 + gamma)*Rlambda*Theta*up)/(dx*(Rlambda - Rlambdap)*(1 + zeta))),(dt*(-1 + gamma)*Rlambda*Theta)/(dx*(Rlambda - Rlambdap)*(1 + zeta)),0,
   (dt*Rlambda*Theta*up*(-(ep*gamma) + (-1 + gamma)*rhop*Power(up,2)))/(dx*rhop*(Rlambda - Rlambdap)*(1 + zeta)),
   (dt*Rlambda*Theta*(2*ep*gamma - 3*(-1 + gamma)*rhop*Power(up,2)))/(2.*dx*rhop*(Rlambda - Rlambdap)*(1 + zeta)),
   (dt*gamma*Rlambda*Theta*up)/(dx*(Rlambda - Rlambdap)*(1 + zeta)),0,(dt*Rlambda*Theta*up*Yp)/(dx*(-Rlambda + Rlambdap)*(1 + zeta)),
   (dt*Rlambda*Theta*Yp)/(dx*(Rlambda - Rlambdap)*(1 + zeta)),0,(dt*Rlambda*Theta*up)/(dx*(Rlambda - Rlambdap)*(1 + zeta));
  c += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Viscous (2nd order derivatives)
///////////////////////////////////////////////////////////////////////////////
#if defined(VISCOUS)

  temp << 0,0,0,0,(2*dt*Pr*Theta*((rhom - rhop)*u + rho*(2*u - um + up)))/(3.*Power(dx,2)*Power(rho,2)*(1 + zeta)),
   (-2*dt*Pr*(2*rho + rhom - rhop)*Theta)/(3.*Power(dx,2)*Power(rho,2)*(1 + zeta)),0,0,
   (dt*Theta*(6*e*gamma*(rho + rhom - rhop) + rho*(-3*em*gamma + 3*ep*gamma + (3*gamma - 4*Pr)*u*((-rhom + rhop)*u - 2*rho*(u - um + up)))))/
    (6.*Power(dx,2)*Power(rho,3)*(1 + zeta)),(dt*(3*gamma - 4*Pr)*Theta*((rhom - rhop)*u + rho*(2*u - um + up)))/
    (6.*Power(dx,2)*Power(rho,2)*(1 + zeta)),-(dt*gamma*(2*rho + rhom - rhop)*Theta)/(2.*Power(dx,2)*Power(rho,2)*(1 + zeta)),0,
   (dt*Theta*((rhom - rhop)*Y + rho*(2*Y - Ym + Yp)))/(2.*Power(dx,2)*Le*Power(rho,2)*(1 + zeta)),0,0,
   -(dt*(2*rho + rhom - rhop)*Theta)/(2.*Power(dx,2)*Le*Power(rho,2)*(1 + zeta));
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
  double Llambda;
  double Llambdam;
  double Rlambda;
  double Rlambdap;
  int E = 0;

  {
  double lambda_1 = u + sqrt(gamma*p/rho);
  double lambda_2 = u_rhoavg(rhom, rho, um, u) + sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  Llambda = std::max(lambda_1,lambda_2);
  }
  {
  double lambda_1 = um - sqrt(gamma*pm/rhom);
  double lambda_2 = u_rhoavg(rhom, rho, um, u) - sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  Llambdam = std::min(lambda_1,lambda_2);
  }
  {
  double lambda_1 = up + sqrt(gamma*pp/rhop);
  double lambda_2 = u_rhoavg(rho, rhop, u, up) + sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  Rlambdap = std::max(lambda_1,lambda_2);
  }
  {
  double lambda_1 = u - sqrt(gamma*p/rho);
  double lambda_2 = u_rhoavg(rho, rhop, u, up) - sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  Rlambda = std::min(lambda_1,lambda_2);
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

  temp << 0,-((dt*Llambda*Theta)/(dx*(Llambda - Llambdam)*(1 + zeta))),0,0,
   -(dt*(-3 + gamma)*Llambda*Theta*Power(um,2))/(2.*dx*(Llambda - Llambdam)*(1 + zeta)),
   (dt*(-3 + gamma)*Llambda*Theta*um)/(dx*(Llambda - Llambdam)*(1 + zeta)),-((dt*(-1 + gamma)*Llambda*Theta)/(dx*(Llambda - Llambdam)*(1 + zeta))),0,
   (dt*Llambda*Theta*(em*gamma*um - (-1 + gamma)*rhom*Power(um,3)))/(dx*(Llambda - Llambdam)*rhom*(1 + zeta)),
   (dt*Llambda*Theta*(-2*em*gamma + 3*(-1 + gamma)*rhom*Power(um,2)))/(2.*dx*(Llambda - Llambdam)*rhom*(1 + zeta)),
   -((dt*gamma*Llambda*Theta*um)/(dx*(Llambda - Llambdam)*(1 + zeta))),0,(dt*Llambda*Theta*um*Ym)/(dx*(Llambda - Llambdam)*(1 + zeta)),
   -((dt*Llambda*Theta*Ym)/(dx*(Llambda - Llambdam)*(1 + zeta))),0,-((dt*Llambda*Theta*um)/(dx*(Llambda - Llambdam)*(1 + zeta)));
  a += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Viscous (2nd order derivatives)
///////////////////////////////////////////////////////////////////////////////
#if defined(VISCOUS)

  temp << 0,0,0,0,(2*dt*Pr*Theta*((-rhom + rhop)*u + rho*(2*u + um - up)))/(3.*Power(dx,2)*Power(rho,2)*(1 + zeta)),
   (-2*dt*Pr*(2*rho - rhom + rhop)*Theta)/(3.*Power(dx,2)*Power(rho,2)*(1 + zeta)),0,0,
   (dt*Theta*(6*e*gamma*(rho - rhom + rhop) + rho*(3*em*gamma - 3*ep*gamma + (3*gamma - 4*Pr)*u*((rhom - rhop)*u - 2*rho*(u + um - up)))))/
    (6.*Power(dx,2)*Power(rho,3)*(1 + zeta)),(dt*(3*gamma - 4*Pr)*Theta*((-rhom + rhop)*u + rho*(2*u + um - up)))/
    (6.*Power(dx,2)*Power(rho,2)*(1 + zeta)),-(dt*gamma*(2*rho - rhom + rhop)*Theta)/(2.*Power(dx,2)*Power(rho,2)*(1 + zeta)),0,
   (dt*Theta*((-rhom + rhop)*Y + rho*(2*Y + Ym - Yp)))/(2.*Power(dx,2)*Le*Power(rho,2)*(1 + zeta)),0,0,
   -(dt*(2*rho - rhom + rhop)*Theta)/(2.*Power(dx,2)*Le*Power(rho,2)*(1 + zeta));
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
  double Llambda;
  double Llambdam;
  double Rlambda;
  double Rlambdap;
  int E = 0;

  {
  double lambda_1 = u + sqrt(gamma*p/rho);
  double lambda_2 = u_rhoavg(rhom, rho, um, u) + sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  Llambda = std::max(lambda_1,lambda_2);
  }
  {
  double lambda_1 = um - sqrt(gamma*pm/rhom);
  double lambda_2 = u_rhoavg(rhom, rho, um, u) - sqrt(gamma*p_rhoavg(rhom, rho, um, u, pm, p, gamma)/rho_rhoavg(rhom, rho));
  Llambdam = std::min(lambda_1,lambda_2);
  }
  {
  double lambda_1 = up + sqrt(gamma*pp/rhop);
  double lambda_2 = u_rhoavg(rho, rhop, u, up) + sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  Rlambdap = std::max(lambda_1,lambda_2);
  }
  {
  double lambda_1 = u - sqrt(gamma*p/rho);
  double lambda_2 = u_rhoavg(rho, rhop, u, up) - sqrt(gamma*p_rhoavg(rho, rhop, u, up, p, pp, gamma)/rho_rhoavg(rho, rhop));
  Rlambda = std::min(lambda_1,lambda_2);
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

  // temp <<(dt*((Llambda*Llambdam*(rho - rhom))/(Llambda - Llambdam) - ((-rho + rhop)*Rlambda*Rlambdap)/(-Rlambda + Rlambdap) +
  //       (-(Llambdam*rho*u) + Llambda*rhom*um)/(Llambda - Llambdam) - (rho*Rlambdap*u - rhop*Rlambda*up)/(-Rlambda + Rlambdap)))/(dx*(1 + zeta)),
  //  (dt*((Llambda*Llambdam*(rho*u - rhom*um))/(Llambda - Llambdam) +
  //       (-(Llambdam*(rho*Power(u,2) + (-1 + gamma)*(e - (rho*Power(u,2))/2.))) + Llambda*(rhom*Power(um,2) + (-1 + gamma)*(em - (rhom*Power(um,2))/2.)))/
  //        (Llambda - Llambdam) - (Rlambda*Rlambdap*(-(rho*u) + rhop*up))/(-Rlambda + Rlambdap) -
  //       (Rlambdap*(rho*Power(u,2) + (-1 + gamma)*(e - (rho*Power(u,2))/2.)) - Rlambda*(rhop*Power(up,2) + (-1 + gamma)*(ep - (rhop*Power(up,2))/2.)))/
  //        (-Rlambda + Rlambdap)))/(dx*(1 + zeta)),(dt*(((e - em)*Llambda*Llambdam)/(Llambda - Llambdam) -
  //       ((-e + ep)*Rlambda*Rlambdap)/(-Rlambda + Rlambdap) + (-(Llambdam*u*(e + (-1 + gamma)*(e - (rho*Power(u,2))/2.))) +
  //          Llambda*um*(em + (-1 + gamma)*(em - (rhom*Power(um,2))/2.)))/(Llambda - Llambdam) -
  //       (Rlambdap*u*(e + (-1 + gamma)*(e - (rho*Power(u,2))/2.)) - Rlambda*up*(ep + (-1 + gamma)*(ep - (rhop*Power(up,2))/2.)))/(-Rlambda + Rlambdap)))/
  //   (dx*(1 + zeta)),(dt*((Llambda*Llambdam*(rho*Y - rhom*Ym))/(Llambda - Llambdam) + (-(Llambdam*rho*u*Y) + Llambda*rhom*um*Ym)/(Llambda - Llambdam) -
  //       (Rlambda*Rlambdap*(-(rho*Y) + rhop*Yp))/(-Rlambda + Rlambdap) - (rho*Rlambdap*u*Y - rhop*Rlambda*up*Yp)/(-Rlambda + Rlambdap)))/(dx*(1 + zeta));
  // rhs += temp;
  auto var_vec_l = Variable_Vector_Isolator<solution_vector_type>(solution_vector_m, gamma);
  auto var_vec = Variable_Vector_Isolator<solution_vector_type>(solution_vector, gamma);
  auto var_vec_r = Variable_Vector_Isolator<solution_vector_type>(solution_vector_p, gamma);

   rhs += (hyperbolic_flux.flux(var_vec_l.w(), var_vec.w(), gamma) - hyperbolic_flux.flux(var_vec.w(), var_vec_r.w(), gamma)) /dx*dt;

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
                  ((-rhom + rhop)*u*(-um + up))/(2.*Power(dx,2)) - (rho*Power(-um + up,2))/(4.*Power(dx,2)) - (rho*u*(-2*u + um + up))/Power(dx,2)))/rho)
           )/(-1 + gamma)))/(1 + zeta),(dt*(-2*Y + Ym + Yp))/(Power(dx,2)*Le*(1 + zeta));
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
