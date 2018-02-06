#ifndef CREATE_IMPLICIT_MATRIX_VECTOR_3_H
#define CREATE_IMPLICIT_MATRIX_VECTOR_3_H

#include "../Explicit_Flux_and_Sources/HLLE.h"
#include "../Explicit_Flux_and_Sources/Centered_Difference.h"
#include "../Explicit_Flux_and_Sources/Sources.h"

#define HYPERBOLIC
#define VISCOUS
// #define SOURCE

template <typename T>
T Power(T num, int expo) {
  return pow(num, expo);
}

double Power(int num, double expo) {
  return exp(expo);
}

template <typename solution_vector_type, typename Matrix_type>
Matrix_type create_mid_band_matrix(const solution_vector_type& solution_vector_mm,
                                       const solution_vector_type solution_vector_m,
                                       const solution_vector_type solution_vector,
                                       const solution_vector_type solution_vector_p,
                                       const solution_vector_type solution_vector_pp,
                                       const double& gamma, const double& Pr,
                                       const double& Le, const double& Q, const double& lambda,
                                       const double& theta, const double& dx, const double& dt,
                                       const double& zeta, const double& Theta) {
  double rho = solution_vector[0];
  double rhom = solution_vector_m[0];
  double rhop = solution_vector_p[0];
  double u = solution_vector[1]/solution_vector[0];
  double um = solution_vector_m[1]/solution_vector_m[0];
  double up = solution_vector_p[1]/solution_vector_p[0];
  double e = solution_vector[2];
  double em = solution_vector_m[2];
  double ep = solution_vector_p[2];
  double Y = solution_vector[3]/solution_vector[0];
  double Ym = solution_vector_m[3]/solution_vector_m[0];
  double Yp = solution_vector_p[3]/solution_vector_p[0];
  int E = 0;

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

  temp << 0,0,0,0,-(dt*(-3 + gamma)*Theta*u*(um - up))/(2.*dx*(1 + zeta)),(dt*(-3 + gamma)*Theta*(um - up))/(2.*dx*(1 + zeta)),0,0,
   (dt*Theta*(em*gamma*rho*u - ep*gamma*rho*u - e*gamma*rhom*u + e*gamma*rhop*u + e*gamma*rho*um + 3*Power(rho,2)*Power(u,2)*um -
        3*gamma*Power(rho,2)*Power(u,2)*um + rho*(-(e*gamma) + 3*(-1 + gamma)*rho*Power(u,2))*up))/(2.*dx*Power(rho,2)*(1 + zeta)),
   (dt*Theta*(-(em*gamma*rho) + ep*gamma*rho + e*gamma*rhom - e*gamma*rhop - 3*Power(rho,2)*u*um + 3*gamma*Power(rho,2)*u*um -
        3*(-1 + gamma)*Power(rho,2)*u*up))/(2.*dx*Power(rho,2)*(1 + zeta)),(dt*gamma*Theta*(-um + up))/(2.*dx*(1 + zeta)),0,
   (dt*Theta*(um*Y - up*Y + u*(Ym - Yp)))/(2.*dx*(1 + zeta)),(dt*Theta*(-Ym + Yp))/(2.*dx*(1 + zeta)),0,(dt*Theta*(-um + up))/(2.*dx*(1 + zeta));
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
                                       const double& Le, const double& Q, const double& lambda,
                                       const double& theta, const double& dx, const double& dt,
                                       const double& zeta, const double& Theta) {
  double rho = solution_vector[0];
  double rhom = solution_vector_m[0];
  double rhop = solution_vector_p[0];
  double u = solution_vector[1]/solution_vector[0];
  double um = solution_vector_m[1]/solution_vector_m[0];
  double up = solution_vector_p[1]/solution_vector_p[0];
  double e = solution_vector[2];
  double em = solution_vector_m[2];
  double ep = solution_vector_p[2];
  double Y = solution_vector[3]/solution_vector[0];
  double Ym = solution_vector_m[3]/solution_vector_m[0];
  double Yp = solution_vector_p[3]/solution_vector_p[0];
  int E = 0;

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

  temp << 0,(dt*Theta)/(2*dx + 2*dx*zeta),0,0,(dt*(-3 + gamma)*Theta*Power(u,2))/(4.*dx*(1 + zeta)),-(dt*(-3 + gamma)*Theta*u)/(2.*dx*(1 + zeta)),
   (dt*(-1 + gamma)*Theta)/(2.*dx*(1 + zeta)),0,(dt*Theta*u*(-(e*gamma) + (-1 + gamma)*rho*Power(u,2)))/(2.*dx*rho*(1 + zeta)),
   (dt*Theta*(2*e*gamma - 3*(-1 + gamma)*rho*Power(u,2)))/(4.*dx*rho*(1 + zeta)),(dt*gamma*Theta*u)/(2*dx + 2*dx*zeta),0,
   -((dt*Theta*u*Y)/(2*dx + 2*dx*zeta)),(dt*Theta*Y)/(2*dx + 2*dx*zeta),0,(dt*Theta*u)/(2*dx + 2*dx*zeta);
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
                                       const double& Le, const double& Q, const double& lambda,
                                       const double& theta, const double& dx, const double& dt,
                                       const double& zeta, const double& Theta) {
  double rho = solution_vector[0];
  double rhom = solution_vector_m[0];
  double rhop = solution_vector_p[0];
  double u = solution_vector[1]/solution_vector[0];
  double um = solution_vector_m[1]/solution_vector_m[0];
  double up = solution_vector_p[1]/solution_vector_p[0];
  double e = solution_vector[2];
  double em = solution_vector_m[2];
  double ep = solution_vector_p[2];
  double Y = solution_vector[3]/solution_vector[0];
  double Ym = solution_vector_m[3]/solution_vector_m[0];
  double Yp = solution_vector_p[3]/solution_vector_p[0];
  int E = 0;

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

  temp << 0,-((dt*Theta)/(2*dx + 2*dx*zeta)),0,0,-(dt*(-3 + gamma)*Theta*Power(u,2))/(4.*dx*(1 + zeta)),(dt*(-3 + gamma)*Theta*u)/(2.*dx*(1 + zeta)),
   -(dt*(-1 + gamma)*Theta)/(2.*dx*(1 + zeta)),0,(dt*Theta*(e*gamma*u - (-1 + gamma)*rho*Power(u,3)))/(2.*dx*rho*(1 + zeta)),
   (dt*Theta*(-2*e*gamma + 3*(-1 + gamma)*rho*Power(u,2)))/(4.*dx*rho*(1 + zeta)),-((dt*gamma*Theta*u)/(2*dx + 2*dx*zeta)),0,
   (dt*Theta*u*Y)/(2*dx + 2*dx*zeta),-((dt*Theta*Y)/(2*dx + 2*dx*zeta)),0,-((dt*Theta*u)/(2*dx + 2*dx*zeta));
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
                                       const double& Le, const double& Q, const double& lambda,
                                       const double& theta, const double& dx, const double& dt,
                                       const double& zeta, const double& Theta,
                                       const solution_vector_type& DeltaUm) {


  double rho = solution_vector[0];
  double rhom = solution_vector_m[0];
  double rhop = solution_vector_p[0];
  double u = solution_vector[1]/solution_vector[0];
  double um = solution_vector_m[1]/solution_vector_m[0];
  double up = solution_vector_p[1]/solution_vector_p[0];
  double e = solution_vector[2];
  double em = solution_vector_m[2];
  double ep = solution_vector_p[2];
  double Y = solution_vector[3]/solution_vector[0];
  double Ym = solution_vector_m[3]/solution_vector_m[0];
  double Yp = solution_vector_p[3]/solution_vector_p[0];
  double dUm1 = DeltaUm[0];
  double dUm2 = DeltaUm[1];
  double dUm3 = DeltaUm[2];
  double dUm4 = DeltaUm[3];
  int E = 0;

  // HLLE<std::vector<solution_vector_type>> hyperbolic_flux;
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

  temp << (dt*(-((-rhom + rhop)*u)/(2.*dx) - (rho*(-um + up))/(2.*dx)))/(1 + zeta),
   (dt*(-((-rhom + rhop)*Power(u,2))/(2.*dx) - (rho*u*(-um + up))/dx -
        (-1 + gamma)*((-em + ep)/(2.*dx) - ((-rhom + rhop)*Power(u,2))/(4.*dx) - (rho*u*(-um + up))/(2.*dx))))/(1 + zeta),
   (dt*(-((e + (-1 + gamma)*(e - (rho*Power(u,2))/2.))*(-um + up))/(2.*dx) -
        u*((-em + ep)/(2.*dx) + (-1 + gamma)*((-em + ep)/(2.*dx) - ((-rhom + rhop)*Power(u,2))/(4.*dx) - (rho*u*(-um + up))/(2.*dx)))))/(1 + zeta),
   (dt*(-((-rhom + rhop)*u*Y)/(2.*dx) - (rho*(-um + up)*Y)/(2.*dx) - (rho*u*(-Ym + Yp))/(2.*dx)))/(1 + zeta);
  rhs += temp;

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

  temp << 0.,0.,(dt*lambda*Q*rho*Y)/(Power(E,(rho*theta)/((-1 + gamma)*(e - (rho*Power(u,2))/2.)))*(1 + zeta)),
   -((dt*lambda*rho*Y)/(Power(E,(rho*theta)/((-1 + gamma)*(e - (rho*Power(u,2))/2.)))*(1 + zeta)));
  rhs += temp;

#endif
  // auto var_vec_l = Variable_Vector_Isolator<solution_vector_type>(solution_vector_m, gamma);
  // auto var_vec = Variable_Vector_Isolator<solution_vector_type>(solution_vector, gamma);
  // auto var_vec_r = Variable_Vector_Isolator<solution_vector_type>(solution_vector_p, gamma);
  //
  //  rhs += (hyperbolic_flux.flux(var_vec_l.w(), var_vec.w(), gamma) - hyperbolic_flux.flux(var_vec.w(), var_vec_r.w(), gamma)) /dx*dt;

  return rhs;
}

#endif //#ifndef CREATE_IMPLICIT_MATRIX_VECTOR_H
