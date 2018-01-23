#ifndef CREATE_IMPLICIT_MATRIX_VECTOR_H
#define CREATE_IMPLICIT_MATRIX_VECTOR_H

#include "../Explicit_Flux_and_Sources/HLLE.h"

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
                                       const double& theta, const double& dx, const double& dt) {
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

  b << 1,0,0,0,(dt*(-4*dx*(-3 + gamma)*Power(rho,3)*u*(um - up) +
        3*Pr*(Power(rhom - rhop,2)*u + 2*Power(rho,2)*(-2*u + um + up) - rho*(rhom*(2*u + um - up) + rhop*(2*u - um + up)))))/
    (8.*Power(dx,2)*Power(rho,3)),(8*Power(dx,2)*Power(rho,3) - 3*dt*Pr*Power(rhom - rhop,2) + 6*dt*Pr*rho*(rhom + rhop) +
      4*dt*dx*(-3 + gamma)*Power(rho,3)*(um - up))/(8.*Power(dx,2)*Power(rho,3)),0,0,
   (dt*((4*(em*gamma*rho*u - ep*gamma*rho*u - e*gamma*rhom*u + e*gamma*rhop*u + e*gamma*rho*um + 3*Power(rho,2)*Power(u,2)*um -
             3*gamma*Power(rho,2)*Power(u,2)*um + rho*(-(e*gamma) + 3*(-1 + gamma)*rho*Power(u,2))*up))/(dx*Power(rho,2)) +
        (-16*e*gamma*rho*rhom + 12*e*gamma*Power(rhom,2) + 8*ep*gamma*rho*(rho + rhom - rhop) - 16*e*gamma*rho*rhop - 24*e*gamma*rhom*rhop +
           12*e*gamma*Power(rhop,2) + 8*em*gamma*rho*(rho - rhom + rhop) + 32*gamma*Power(rho,3)*Power(u,2) - 24*Pr*Power(rho,3)*Power(u,2) +
           8*gamma*Power(rho,2)*rhom*Power(u,2) - 6*Pr*Power(rho,2)*rhom*Power(u,2) - 4*gamma*rho*Power(rhom,2)*Power(u,2) +
           3*Pr*rho*Power(rhom,2)*Power(u,2) + 8*gamma*Power(rho,2)*rhop*Power(u,2) - 6*Pr*Power(rho,2)*rhop*Power(u,2) +
           8*gamma*rho*rhom*rhop*Power(u,2) - 6*Pr*rho*rhom*rhop*Power(u,2) - 4*gamma*rho*Power(rhop,2)*Power(u,2) +
           3*Pr*rho*Power(rhop,2)*Power(u,2) - 16*gamma*Power(rho,3)*u*um + 12*Pr*Power(rho,3)*u*um + 8*gamma*Power(rho,2)*rhom*u*um -
           6*Pr*Power(rho,2)*rhom*u*um - 8*gamma*Power(rho,2)*rhop*u*um + 6*Pr*Power(rho,2)*rhop*u*um - 4*gamma*Power(rho,3)*Power(um,2) +
           3*Pr*Power(rho,3)*Power(um,2) - 2*(4*gamma - 3*Pr)*Power(rho,2)*((2*rho + rhom - rhop)*u - rho*um)*up +
           (-4*gamma + 3*Pr)*Power(rho,3)*Power(up,2))/(Power(dx,2)*Power(rho,4)) -
        (8*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Q*
           (4*e*(e*(-1 + gamma) - rho*theta) + 4*rho*(e - e*gamma + rho*theta)*Power(u,2) + (-1 + gamma)*Power(rho,2)*Power(u,4))*Y)/
         ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2))))/8.,(dt*
      ((4*rho*(-(em*gamma*rho) + ep*gamma*rho + e*gamma*rhom - e*gamma*rhop - 3*Power(rho,2)*u*um + 3*gamma*Power(rho,2)*u*um -
             3*(-1 + gamma)*Power(rho,2)*u*up))/dx + ((4*gamma - 3*Pr)*
           (Power(rhom - rhop,2)*u + 2*Power(rho,2)*(-2*u + um + up) - rho*(rhom*(2*u + um - up) + rhop*(2*u - um + up))))/Power(dx,2) +
        (32*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Q*Power(rho,5)*theta*u*Y)/((-1 + gamma)*Power(-2*e + rho*Power(u,2),2))
        ))/(8.*Power(rho,3)),1 - (dt*gamma*(Power(rhom - rhop,2) - 2*rho*(rhom + rhop) + dx*Power(rho,3)*(um - up)))/(2.*Power(dx,2)*Power(rho,3)) -
    (4*dt*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Q*Power(rho,2)*theta*Y)/((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)),
   -(dt*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Q*rho),
   (dt*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*
      (4*e*(e*(-1 + gamma) - rho*theta) + 4*rho*(e - e*gamma + rho*theta)*Power(u,2) + (-1 + gamma)*Power(rho,2)*Power(u,4))*Y)/
    ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)),dt*((-4*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Power(rho,2)*theta*u*Y)/
       ((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)) + (-Ym + Yp)/(2.*dx)),
   (4*dt*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda*Power(rho,2)*theta*Y)/((-1 + gamma)*Power(-2*e + rho*Power(u,2),2)),
   1 + (dt*(4 + dx*Le*((-rhom + rhop)*u + rho*(2*dx*exp((2*rho*theta)/((-1 + gamma)*(-2*e + rho*Power(u,2))))*lambda - um + up))))/
     (2.*Power(dx,2)*Le);

  return b;
}

template <typename solution_vector_type, typename Matrix_type>
Matrix_type create_top_band_matrix(const solution_vector_type solution_vector_m,
                                       const solution_vector_type solution_vector,
                                       const solution_vector_type solution_vector_p,
                                       const double& gamma, const double& Pr,
                                       const double& Le, const double& Q, const double& lambda,
                                       const double& theta, const double& dx, const double& dt) {
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

  c <<0,dt/(2.*dx),0,0,(dt*(2*dx*(-3 + gamma)*Power(rho,2)*Power(u,2) + 3*Pr*((rhom - rhop)*u + rho*(2*u - um + up))))/(8.*Power(dx,2)*Power(rho,2)),
   -(dt*(3*Pr*(2*rho + rhom - rhop) + 4*dx*(-3 + gamma)*Power(rho,2)*u))/(8.*Power(dx,2)*Power(rho,2)),(dt*(-1 + gamma))/(2.*dx),0,
   (dt*(8*e*gamma*(rho + rhom - rhop) - 4*dx*e*gamma*Power(rho,2)*u +
        rho*(-4*em*gamma + 4*ep*gamma + u*(3*Pr*(rhom - rhop)*u - 4*dx*Power(rho,2)*Power(u,2) + 6*Pr*rho*(u - um + up) +
              4*gamma*((-rhom + rhop)*u + dx*Power(rho,2)*Power(u,2) - 2*rho*(u - um + up))))))/(8.*Power(dx,2)*Power(rho,3)),
   (dt*(2*dx*rho*(2*e*gamma - 3*(-1 + gamma)*rho*Power(u,2)) + (4*gamma - 3*Pr)*((rhom - rhop)*u + rho*(2*u - um + up))))/(8.*Power(dx,2)*Power(rho,2)),
   (dt*gamma*(-rhom + rhop + rho*(-2 + dx*rho*u)))/(2.*Power(dx,2)*Power(rho,2)),0,0,(dt*Y)/(2.*dx),0,(dt*(-2 + dx*Le*rho*u))/(2.*Power(dx,2)*Le);

  return c;
}

template <typename solution_vector_type, typename Matrix_type>
Matrix_type create_bot_band_matrix(const solution_vector_type solution_vector_m,
                                       const solution_vector_type solution_vector,
                                       const solution_vector_type solution_vector_p,
                                       const double& gamma, const double& Pr,
                                       const double& Le, const double& Q, const double& lambda,
                                       const double& theta, const double& dx, const double& dt) {
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

    a << 0,-dt/(2.*dx),0,0,(dt*(-2*dx*(-3 + gamma)*Power(rho,2)*Power(u,2) + 3*Pr*((-rhom + rhop)*u + rho*(2*u + um - up))))/(8.*Power(dx,2)*Power(rho,2)),
   (dt*(-3*Pr*(2*rho - rhom + rhop) + 4*dx*(-3 + gamma)*Power(rho,2)*u))/(8.*Power(dx,2)*Power(rho,2)),-(dt*(-1 + gamma))/(2.*dx),0,
   (dt*(4*e*gamma*(2*(rho - rhom + rhop) + dx*Power(rho,2)*u) +
        rho*(4*em*gamma - 4*ep*gamma + u*(3*Pr*(-rhom + rhop)*u + 4*dx*Power(rho,2)*Power(u,2) -
              4*gamma*((-rhom + rhop)*u + dx*Power(rho,2)*Power(u,2) + 2*rho*(u + um - up)) + 6*Pr*rho*(u + um - up)))))/(8.*Power(dx,2)*Power(rho,3)),
   (dt*(2*dx*rho*(-2*e*gamma + 3*(-1 + gamma)*rho*Power(u,2)) + (4*gamma - 3*Pr)*((-rhom + rhop)*u + rho*(2*u + um - up))))/
    (8.*Power(dx,2)*Power(rho,2)),-(dt*gamma*(-rhom + rhop + rho*(2 + dx*rho*u)))/(2.*Power(dx,2)*Power(rho,2)),0,0,-(dt*Y)/(2.*dx),0,
   -(dt*(2 + dx*Le*rho*u))/(2.*Power(dx,2)*Le);

  return a;
}

template <typename solution_vector_type, typename Matrix_type>
solution_vector_type create_rhs_vector(const solution_vector_type solution_vector_mm,
                                       const solution_vector_type solution_vector_m,
                                       const solution_vector_type solution_vector,
                                       const solution_vector_type solution_vector_p,
                                       const solution_vector_type solution_vector_pp,
                                       const double& gamma, const double& Pr,
                                       const double& Le, const double& Q, const double& lambda,
                                       const double& theta, const double& dx, const double& dt) {

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

  solution_vector_type rhs;
  solution_vector_type rhs2;
  // using global_solution_vector_type = std::vector<solution_vector_type>;
  // HLLE<global_solution_vector_type> hyperbolic_flux;

  rhs << dt*(0. - (-rhom/(2.*dx) + rhop/(2.*dx))*u - rho*(-um/(2.*dx) + up/(2.*dx))),
   dt*(0. - (-rhom/(2.*dx) + rhop/(2.*dx))*Power(u,2) + (3*Pr*((-2*u)/Power(dx,2) + um/Power(dx,2) + up/Power(dx,2)))/4. -
      2*rho*u*(-um/(2.*dx) + up/(2.*dx)) - (-1 + gamma)*(-em/(2.*dx) + ep/(2.*dx) - ((-rhom/(2.*dx) + rhop/(2.*dx))*Power(u,2))/2. -
         rho*u*(-um/(2.*dx) + up/(2.*dx)))),dt*((3*Pr*u*((-2*u)/Power(dx,2) + um/Power(dx,2) + up/Power(dx,2)))/4. -
      (e + (-1 + gamma)*(e - (rho*Power(u,2))/2.))*(-um/(2.*dx) + up/(2.*dx)) + (3*Pr*Power(-um/(2.*dx) + up/(2.*dx),2))/4. -
      u*(-em/(2.*dx) + ep/(2.*dx) + (-1 + gamma)*(-em/(2.*dx) + ep/(2.*dx) - ((-rhom/(2.*dx) + rhop/(2.*dx))*Power(u,2))/2. -
            rho*u*(-um/(2.*dx) + up/(2.*dx)))) + (gamma*(-(((-1 + gamma)*((-2*rho)/Power(dx,2) + rhom/Power(dx,2) + rhop/Power(dx,2))*
                (e - (rho*Power(u,2))/2.))/Power(rho,2)) + (2*(-1 + gamma)*Power(-rhom/(2.*dx) + rhop/(2.*dx),2)*(e - (rho*Power(u,2))/2.))/
            Power(rho,3) - (2*(-1 + gamma)*(-rhom/(2.*dx) + rhop/(2.*dx))*
              (-em/(2.*dx) + ep/(2.*dx) - ((-rhom/(2.*dx) + rhop/(2.*dx))*Power(u,2))/2. - rho*u*(-um/(2.*dx) + up/(2.*dx))))/Power(rho,2) +
           ((-1 + gamma)*((-2*e)/Power(dx,2) + em/Power(dx,2) + ep/Power(dx,2) -
                (((-2*rho)/Power(dx,2) + rhom/Power(dx,2) + rhop/Power(dx,2))*Power(u,2))/2. -
                rho*u*((-2*u)/Power(dx,2) + um/Power(dx,2) + up/Power(dx,2)) - 2*(-rhom/(2.*dx) + rhop/(2.*dx))*u*(-um/(2.*dx) + up/(2.*dx)) -
                rho*Power(-um/(2.*dx) + up/(2.*dx),2)))/rho))/(-1 + gamma) +
      (lambda*Q*rho*Y)/exp((rho*theta)/((-1 + gamma)*(e - (rho*Power(u,2))/2.)))),
   dt*(-((lambda*rho*Y)/exp((rho*theta)/((-1 + gamma)*(e - (rho*Power(u,2))/2.)))) - (-rhom/(2.*dx) + rhop/(2.*dx))*u*Y -
      rho*(-um/(2.*dx) + up/(2.*dx))*Y + ((-2*Y)/Power(dx,2) + Ym/Power(dx,2) + Yp/Power(dx,2))/Le - rho*u*(-Ym/(2.*dx) + Yp/(2.*dx)));

  rhs2 << 0.2/8.0*(solution_vector_pp-4.0*solution_vector_p + 6.0*solution_vector - 4.0*solution_vector_m+solution_vector_mm);
  rhs -= rhs2;

      // std::cout << "rhs2" << rhs2 << std::endl;

  return rhs;
}

#endif //#ifndef CREATE_IMPLICIT_MATRIX_VECTOR_H
