#ifndef CREATE_IMPLICIT_MATRIX_VECTOR_3_H
#define CREATE_IMPLICIT_MATRIX_VECTOR_3_H

#include "../Explicit_Flux_and_Sources/HLLE.h"
#include "../Explicit_Flux_and_Sources/Centered_Difference.h"
#include "../Explicit_Flux_and_Sources/Sources.h"

// #define HYPERBOLIC
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

  temp << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
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

  temp << 0,(dt*Theta)/(2*dx + 2*dx*zeta),0,0,(dt*(-3 + gamma)*Theta*Power(up,2))/(4.*dx*(1 + zeta)),-(dt*(-3 + gamma)*Theta*up)/(2.*dx*(1 + zeta)),
   (dt*(-1 + gamma)*Theta)/(2.*dx*(1 + zeta)),0,(dt*Theta*up*(-(ep*gamma) + (-1 + gamma)*rho*Power(up,2)))/(2.*dx*rho*(1 + zeta)),
   (dt*Theta*(2*ep*gamma - 3*(-1 + gamma)*rho*Power(up,2)))/(4.*dx*rho*(1 + zeta)),(dt*gamma*Theta*up)/(2*dx + 2*dx*zeta),0,
   -((dt*Theta*up*Yp)/(2*dx + 2*dx*zeta)),(dt*Theta*Yp)/(2*dx + 2*dx*zeta),0,(dt*Theta*up)/(2*dx + 2*dx*zeta);
  c += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Viscous (2nd order derivatives)
///////////////////////////////////////////////////////////////////////////////
#if defined(VISCOUS)

  temp << 0,0,0,0,(dt*Pr*Theta*(-(rho*um) + rhom*(4*u + up)))/(3.*Power(dx,2)*rho*rhom*(1 + zeta)),
   (dt*Pr*(rho - 5*rhom)*Theta)/(3.*Power(dx,2)*rho*rhom*(1 + zeta)),0,0,
   (dt*Theta*(-3*em*gamma*Power(rho,3) + rhom*(9*e*gamma*rho*rhom + 3*ep*gamma*rhom*(3*rho - 2*rhopp) +
           rho*(3*epp*gamma*rhom + (3*gamma - 4*Pr)*rho*(-4*rhom*Power(u,2) + rho*Power(um,2) + rhom*u*up - rhom*up*(up + upp))))))/
    (12.*Power(dx,2)*Power(rho,3)*Power(rhom,2)*(1 + zeta)),(dt*(3*gamma - 4*Pr)*Theta*(-(rho*um) + rhom*(3*u + up + upp)))/
    (12.*Power(dx,2)*rho*rhom*(1 + zeta)),(dt*gamma*(Power(rho,2) - 6*rho*rhom + rhom*rhopp)*Theta)/(4.*Power(dx,2)*Power(rho,2)*rhom*(1 + zeta)),0,
   (4*dt*rhom*Theta*Y - dt*rho*Theta*Ym + dt*rhom*Theta*Yp)/(4*Power(dx,2)*Le*rho*rhom + 4*Power(dx,2)*Le*rho*rhom*zeta),0,0,
   (dt*rho*Theta - 5*dt*rhom*Theta)/(4*Power(dx,2)*Le*rho*rhom + 4*Power(dx,2)*Le*rho*rhom*zeta);
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

  temp << 0,-((dt*Theta)/(2*dx + 2*dx*zeta)),0,0,-(dt*(-3 + gamma)*Theta*Power(um,2))/(4.*dx*(1 + zeta)),(dt*(-3 + gamma)*Theta*um)/(2.*dx*(1 + zeta)),
   -(dt*(-1 + gamma)*Theta)/(2.*dx*(1 + zeta)),0,(dt*Theta*(em*gamma*um - (-1 + gamma)*rhom*Power(um,3)))/(2.*dx*rhom*(1 + zeta)),
   (dt*Theta*(-2*em*gamma + 3*(-1 + gamma)*rhom*Power(um,2)))/(4.*dx*rhom*(1 + zeta)),-((dt*gamma*Theta*um)/(2*dx + 2*dx*zeta)),0,
   (dt*Theta*um*Ym)/(2*dx + 2*dx*zeta),-((dt*Theta*Ym)/(2*dx + 2*dx*zeta)),0,-((dt*Theta*um)/(2*dx + 2*dx*zeta));
  a += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Viscous (2nd order derivatives)
///////////////////////////////////////////////////////////////////////////////
#if defined(VISCOUS)

  temp << 0,0,0,0,(dt*Pr*Theta*(4*rhom*u + rho*um - rhom*up))/(3.*Power(dx,2)*rho*rhom*(1 + zeta)),
   -(dt*Pr*(rho + 3*rhom)*Theta)/(3.*Power(dx,2)*rho*rhom*(1 + zeta)),0,0,
   (dt*Theta*(3*em*gamma*Power(rho,3) + rhom*(15*e*gamma*rho*rhom + ep*gamma*(-9*rho*rhom + 6*rhom*rhopp) +
           rho*(-3*epp*gamma*rhom - (3*gamma - 4*Pr)*rho*(rho*Power(um,2) + rhom*(4*Power(u,2) + u*up - up*(up + upp)))))))/
    (12.*Power(dx,2)*Power(rho,3)*Power(rhom,2)*(1 + zeta)),(dt*(3*gamma - 4*Pr)*Theta*(rho*um + rhom*(5*u - up - upp)))/
    (12.*Power(dx,2)*rho*rhom*(1 + zeta)),-(dt*gamma*(Power(rho,2) + 2*rho*rhom + rhom*rhopp)*Theta)/(4.*Power(dx,2)*Power(rho,2)*rhom*(1 + zeta)),0,
   (4*dt*rhom*Theta*Y + dt*rho*Theta*Ym - dt*rhom*Theta*Yp)/(4*Power(dx,2)*Le*rho*rhom + 4*Power(dx,2)*Le*rho*rhom*zeta),0,0,
   -((dt*rho*Theta + 3*dt*rhom*Theta)/(4*Power(dx,2)*Le*rho*rhom + 4*Power(dx,2)*Le*rho*rhom*zeta));
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
  int E = 0;
  double dUm1 = DeltaUm[0];
  double dUm2 = DeltaUm[1];
  double dUm3 = DeltaUm[2];
  double dUm4 = DeltaUm[3];
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
