#ifndef CREATE_IMPLICIT_MATRIX_VECTOR_3_H
#define CREATE_IMPLICIT_MATRIX_VECTOR_3_H

#include "../Explicit_Flux_and_Sources/HLLE.h"
#include "../Explicit_Flux_and_Sources/Centered_Difference.h"
#include "../Explicit_Flux_and_Sources/Sources.h"



template <typename T>
T Power(T num, int expo) {
  return pow(num, expo);
}

template <typename solution_vector_type, typename Matrix_type>
Matrix_type create_mid_band_matrix(const solution_vector_type& solution_vector_m,
                                       const solution_vector_type solution_vector,
                                       const solution_vector_type solution_vector_p,
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

  Matrix_type b;

  b << 1,0,0,0,(dt*Theta*(-3*dx*(-3 + gamma)*Power(rho,3)*u*(um - up) +
        4*Pr*(Power(rhom - rhop,2)*u + 2*Power(rho,2)*(-2*u + um + up) - rho*(rhom*(2*u + um - up) + rhop*(2*u - um + up)))))/
    (6.*Power(dx,2)*Power(rho,3)*(1 + zeta)),1 + (dt*Theta*(-4*Pr*Power(rhom - rhop,2) + 8*Pr*rho*(rhom + rhop) +
         3*dx*(-3 + gamma)*Power(rho,3)*(um - up)))/(6.*Power(dx,2)*Power(rho,3)*(1 + zeta)),0,0,
   (dt*Theta*((6*gamma*(2*rho*(ep*(rho + rhom - rhop) + em*(rho - rhom + rhop)) + e*(4*Power(rho,2) + 3*Power(rhom - rhop,2) - 4*rho*(rhom + rhop))))/
         (Power(dx,2)*Power(rho,4)) - (2*(3*gamma - 4*Pr)*(4*Power(rho,2) + Power(rhom - rhop,2) - 2*rho*(rhom + rhop))*Power(u,2))/
         (Power(dx,2)*Power(rho,3)) + (3*(rhom - rhop)*u*(-2*e*gamma + 3*(-1 + gamma)*rho*Power(u,2)))/(dx*Power(rho,2)) -
        (8*(3*e*gamma + (-3*gamma + 4*Pr)*rho*Power(u,2)))/(Power(dx,2)*Power(rho,2)) +
        (4*(3*gamma - 4*Pr)*(rhom - rhop)*u*(um - up))/(Power(dx,2)*Power(rho,2)) + (9*(-1 + gamma)*Power(u,2)*(um - up))/dx +
        (3*(2*e*gamma - 3*(-1 + gamma)*rho*Power(u,2))*(um - up))/(dx*rho) +
        (2*(3*gamma - 4*Pr)*(8*Power(u,2) - Power(um - up,2) - 4*u*(um + up)))/(Power(dx,2)*rho) +
        (3*u*(2*em*gamma - 2*ep*gamma + 3*(-1 + gamma)*u*(-(rhom*u) + rhop*u - 2*rho*um + 2*rho*up)))/(dx*rho) +
        (48*e*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Q*rho*theta*Y)/((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)) -
        (48*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Q*Power(rho,2)*theta*Power(u,2)*Y)/
         ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2))))/(12.*(1 + zeta)),
   (dt*Theta*((6*rho*(-(em*gamma*rho) + ep*gamma*rho + e*gamma*rhom - e*gamma*rhop - 3*Power(rho,2)*u*um + 3*gamma*Power(rho,2)*u*um -
             3*(-1 + gamma)*Power(rho,2)*u*up))/dx + (2*(3*gamma - 4*Pr)*
           (Power(rhom - rhop,2)*u + 2*Power(rho,2)*(-2*u + um + up) - rho*(rhom*(2*u + um - up) + rhop*(2*u - um + up))))/Power(dx,2) +
        (48*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Q*Power(rho,5)*theta*u*Y)/((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)))
      )/(12.*Power(rho,3)*(1 + zeta)),1 + (dt*Theta*(-((gamma*(Power(rhom - rhop,2) - 2*rho*(rhom + rhop) + dx*Power(rho,3)*(um - up)))/
            (Power(dx,2)*Power(rho,3))) - (8*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Q*Power(rho,2)*theta*Y)/
          ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2))))/(2.*(1 + zeta)),
   -((dt*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Q*Theta)/(1 + zeta)),
   (dt*Theta*((8*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*rho*theta*(-e + rho*Power(u,2))*Y)/
         ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)) + (um*Y - up*Y + u*Ym - u*Yp)/dx +
        (Power(rhom - rhop,2)*Y + 2*Power(rho,2)*(-2*Y + Ym + Yp) - rho*(rhom*(2*Y + Ym - Yp) + rhop*(2*Y - Ym + Yp)))/(Power(dx,2)*Le*Power(rho,3))))/
    (2.*(1 + zeta)),(dt*Theta*((-4*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Power(rho,2)*theta*u*Y)/
         ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)) + (-Ym + Yp)/(2.*dx)))/(1 + zeta),
   (4*dt*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Power(rho,2)*theta*Theta*Y)/
    ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)*(1 + zeta)),1 +
    (dt*Theta*(2*Power(dx,2)*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Le*Power(rho,3) - Power(rhom - rhop,2) +
         2*rho*(rhom + rhop) + dx*Le*Power(rho,3)*(-um + up)))/(2.*Power(dx,2)*Le*Power(rho,3)*(1 + zeta));

  return b;
}

template <typename solution_vector_type, typename Matrix_type>
Matrix_type create_top_band_matrix(const solution_vector_type solution_vector_m,
                                       const solution_vector_type solution_vector,
                                       const solution_vector_type solution_vector_p,
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

  Matrix_type c;

  c << 0,(dt*Theta)/(2*dx + 2*dx*zeta),0,0,(dt*Theta*(3*dx*(-3 + gamma)*Power(rho,2)*Power(u,2) + 8*Pr*((rhom - rhop)*u + rho*(2*u - um + up))))/
    (12.*Power(dx,2)*Power(rho,2)*(1 + zeta)),-(dt*Theta*(4*Pr*(2*rho + rhom - rhop) + 3*dx*(-3 + gamma)*Power(rho,2)*u))/
    (6.*Power(dx,2)*Power(rho,2)*(1 + zeta)),(dt*(-1 + gamma)*Theta)/(2.*dx*(1 + zeta)),0,
   (dt*Theta*(6*e*gamma*(rho + rhom - rhop) - 3*dx*e*gamma*Power(rho,2)*u +
        rho*(-3*em*gamma + 3*ep*gamma + u*(4*Pr*(rhom - rhop)*u - 3*dx*Power(rho,2)*Power(u,2) + 8*Pr*rho*(u - um + up) +
              3*gamma*((-rhom + rhop)*u + dx*Power(rho,2)*Power(u,2) - 2*rho*(u - um + up))))))/(6.*Power(dx,2)*Power(rho,3)*(1 + zeta)),
   (dt*Theta*(3*dx*rho*(2*e*gamma - 3*(-1 + gamma)*rho*Power(u,2)) + 2*(3*gamma - 4*Pr)*((rhom - rhop)*u + rho*(2*u - um + up))))/
    (12.*Power(dx,2)*Power(rho,2)*(1 + zeta)),(dt*gamma*Theta*(-rhom + rhop + rho*(-2 + dx*rho*u)))/(2.*Power(dx,2)*Power(rho,2)*(1 + zeta)),0,
   (dt*Theta*((rhom - rhop)*Y - dx*Le*Power(rho,2)*u*Y + rho*(2*Y - Ym + Yp)))/(2.*Power(dx,2)*Le*Power(rho,2)*(1 + zeta)),
   (dt*Theta*Y)/(2*dx + 2*dx*zeta),0,(dt*Theta*(-rhom + rhop + rho*(-2 + dx*Le*rho*u)))/(2.*Power(dx,2)*Le*Power(rho,2)*(1 + zeta));

  return c;
}

template <typename solution_vector_type, typename Matrix_type>
Matrix_type create_bot_band_matrix(const solution_vector_type solution_vector_m,
                                       const solution_vector_type solution_vector,
                                       const solution_vector_type solution_vector_p,
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

  Matrix_type a;

    a << 0,-((dt*Theta)/(2*dx + 2*dx*zeta)),0,0,(dt*Theta*(-3*dx*(-3 + gamma)*Power(rho,2)*Power(u,2) + 8*Pr*((-rhom + rhop)*u + rho*(2*u + um - up))))/
    (12.*Power(dx,2)*Power(rho,2)*(1 + zeta)),(dt*Theta*(-4*Pr*(2*rho - rhom + rhop) + 3*dx*(-3 + gamma)*Power(rho,2)*u))/
    (6.*Power(dx,2)*Power(rho,2)*(1 + zeta)),-(dt*(-1 + gamma)*Theta)/(2.*dx*(1 + zeta)),0,
   (dt*Theta*(3*e*gamma*(2*(rho - rhom + rhop) + dx*Power(rho,2)*u) +
        rho*(3*em*gamma - 3*ep*gamma + u*(4*Pr*(-rhom + rhop)*u + 3*dx*Power(rho,2)*Power(u,2) -
              3*gamma*((-rhom + rhop)*u + dx*Power(rho,2)*Power(u,2) + 2*rho*(u + um - up)) + 8*Pr*rho*(u + um - up)))))/
    (6.*Power(dx,2)*Power(rho,3)*(1 + zeta)),(dt*Theta*(3*dx*rho*(-2*e*gamma + 3*(-1 + gamma)*rho*Power(u,2)) +
        2*(3*gamma - 4*Pr)*((-rhom + rhop)*u + rho*(2*u + um - up))))/(12.*Power(dx,2)*Power(rho,2)*(1 + zeta)),
   -(dt*gamma*Theta*(-rhom + rhop + rho*(2 + dx*rho*u)))/(2.*Power(dx,2)*Power(rho,2)*(1 + zeta)),0,
   (dt*Theta*((-rhom + rhop)*Y + dx*Le*Power(rho,2)*u*Y + rho*(2*Y + Ym - Yp)))/(2.*Power(dx,2)*Le*Power(rho,2)*(1 + zeta)),
   -((dt*Theta*Y)/(2*dx + 2*dx*zeta)),0,-(dt*Theta*(-rhom + rhop + rho*(2 + dx*Le*rho*u)))/(2.*Power(dx,2)*Le*Power(rho,2)*(1 + zeta));

  return a;
}

template <typename solution_vector_type, typename Matrix_type>
solution_vector_type create_rhs_vector(const solution_vector_type solution_vector_m,
                                       const solution_vector_type solution_vector,
                                       const solution_vector_type solution_vector_p,
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


  solution_vector_type rhs;

  rhs << (dt*(0. - ((-rhom + rhop)*u)/(2.*dx) - (rho*(-um + up))/(2.*dx)))/(1 + zeta) + (dUm1*zeta)/(1 + zeta),
   (dt*(0. - ((-rhom + rhop)*Power(u,2))/(2.*dx) - (rho*u*(-um + up))/dx + (4*Pr*(-2*u + um + up))/(3.*Power(dx,2)) -
         (-1 + gamma)*((-em + ep)/(2.*dx) - ((-rhom + rhop)*Power(u,2))/(4.*dx) - (rho*u*(-um + up))/(2.*dx))))/(1 + zeta) + (dUm2*zeta)/(1 + zeta),
   (dt*(-((e + (-1 + gamma)*(e - (rho*Power(u,2))/2.))*(-um + up))/(2.*dx) + (Pr*Power(-um + up,2))/(3.*Power(dx,2)) +
         (4*Pr*u*(-2*u + um + up))/(3.*Power(dx,2)) - u*((-em + ep)/(2.*dx) +
            (-1 + gamma)*((-em + ep)/(2.*dx) - ((-rhom + rhop)*Power(u,2))/(4.*dx) - (rho*u*(-um + up))/(2.*dx))) +
         (gamma*(((-1 + gamma)*Power(-rhom + rhop,2)*(e - (rho*Power(u,2))/2.))/(2.*Power(dx,2)*Power(rho,3)) -
              ((-1 + gamma)*(-2*rho + rhom + rhop)*(e - (rho*Power(u,2))/2.))/(Power(dx,2)*Power(rho,2)) -
              ((-1 + gamma)*(-rhom + rhop)*((-em + ep)/(2.*dx) - ((-rhom + rhop)*Power(u,2))/(4.*dx) - (rho*u*(-um + up))/(2.*dx)))/(dx*Power(rho,2)) +
              ((-1 + gamma)*((-2*e + em + ep)/Power(dx,2) - ((-2*rho + rhom + rhop)*Power(u,2))/(2.*Power(dx,2)) -
                   ((-rhom + rhop)*u*(-um + up))/(2.*Power(dx,2)) - (rho*Power(-um + up,2))/(4.*Power(dx,2)) - (rho*u*(-2*u + um + up))/Power(dx,2)))/rho
              ))/(-1 + gamma) + (lambda*Q*rho*Y)/exp((rho*theta)/((-1 + gamma)*(e - (rho*Power(u,2))/2.)))))/(1 + zeta) + (dUm3*zeta)/(1 + zeta),
   (dt*(-((lambda*rho*Y)/exp((rho*theta)/((-1 + gamma)*(e - (rho*Power(u,2))/2.)))) - ((-rhom + rhop)*u*Y)/(2.*dx) - (rho*(-um + up)*Y)/(2.*dx) -
         (rho*u*(-Ym + Yp))/(2.*dx) + (-2*Y + Ym + Yp)/(Power(dx,2)*Le)))/(1 + zeta) + (dUm4*zeta)/(1 + zeta);

  return rhs;
}

#endif //#ifndef CREATE_IMPLICIT_MATRIX_VECTOR_H
