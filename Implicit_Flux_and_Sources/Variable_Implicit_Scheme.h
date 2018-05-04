#ifndef VARIABLE_IMPLICIT_SCHEME_H
#define VARIABLE_IMPLICIT_SCHEME_H

template <typename grid_type, typename flow_type>
class Variable_Implicit_Scheme {
using global_solution_vector_type = typename grid_type::global_solution_vector_type;
using solution_vector_type = typename global_solution_vector_type::value_type;
using scalar_type = typename grid_type::scalar_type;
using size_type = typename grid_type::size_type;
using matrix_type = typename grid_type::matrix_type;

 public:
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  Variable_Implicit_Scheme() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  Variable_Implicit_Scheme(const Variable_Implicit_Scheme&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  Variable_Implicit_Scheme(Variable_Implicit_Scheme&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  Variable_Implicit_Scheme& operator=(const Variable_Implicit_Scheme&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  Variable_Implicit_Scheme& operator=(Variable_Implicit_Scheme&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Constructor setting up required inputs.
  Variable_Implicit_Scheme( const grid_type grid, const flow_type flow,
                            const global_solution_vector_type& DeltaUm, const scalar_type dt_in,
                            const scalar_type Theta_in, const scalar_type zeta_in,
                            const size_type cell_index) :
                            gamma(flow.gamma), Pr(flow.Pr), Le(flow.Le), Q(flow.Q()), lambda(flow.lambda),
                            theta(flow.theta()), dx(grid.dx()), dt(dt_in), zeta(zeta_in),
                            Theta(Theta_in), T_ignition(flow.T_ignition()) {
  solution_vector_type solution_vector_m   = grid.global_solution_vector[cell_index-1];
  solution_vector_type solution_vector     = grid.global_solution_vector[cell_index];
  solution_vector_type solution_vector_p   = grid.global_solution_vector[cell_index+1];
  rho = solution_vector[0];
  rhom = solution_vector_m[0];
  rhop = solution_vector_p[0];
  u = solution_vector[1]/solution_vector[0];
  um = solution_vector_m[1]/solution_vector_m[0];
  up = solution_vector_p[1]/solution_vector_p[0];
  e = solution_vector[2];
  em = solution_vector_m[2];
  ep = solution_vector_p[2];
  Y = solution_vector[3]/solution_vector[0];
  Ym = solution_vector_m[3]/solution_vector_m[0];
  Yp = solution_vector_p[3]/solution_vector_p[0];
  dUm1 = DeltaUm[cell_index-1][0];
  dUm2 = DeltaUm[cell_index-1][1];
  dUm3 = DeltaUm[cell_index-1][2];
  dUm4 = DeltaUm[cell_index-1][3];
  }

  matrix_type top_matrix();
  matrix_type mid_matrix();
  matrix_type bot_matrix();
  solution_vector_type rhs_matrix();

 private:
  scalar_type rho, rhom, rhop;
  scalar_type u, um, up;
  scalar_type e, em, ep;
  scalar_type Y, Ym, Yp;
  char E;
  scalar_type gamma, Pr, Le, Q;
  scalar_type lambda, theta;
  scalar_type dx, dt;
  scalar_type zeta, Theta;
  scalar_type dUm1, dUm2, dUm3, dUm4;
  scalar_type mf, T_ignition;

  template <typename T>
  scalar_type Power(const T num, const int expo) {
    return pow(num, expo);
  }

  scalar_type Power(const char, const scalar_type expo) {
    return exp(expo);
  }
};


template <typename grid_type, typename flow_type>
typename grid_type::matrix_type Variable_Implicit_Scheme<grid_type, flow_type>::mid_matrix() {

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

   if((gamma-1.0)/ rho * (e-0.5*rho*u*u) > T_ignition){
    b += temp;
  }

#endif

  return b;
}

template <typename grid_type, typename flow_type>
typename grid_type::matrix_type Variable_Implicit_Scheme<grid_type, flow_type>::top_matrix() {

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

template <typename grid_type, typename flow_type>
typename grid_type::matrix_type Variable_Implicit_Scheme<grid_type, flow_type>::bot_matrix() {

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

template <typename grid_type, typename flow_type>
typename grid_type::global_solution_vector_type::value_type
Variable_Implicit_Scheme<grid_type, flow_type>::rhs_matrix() {

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

   if((gamma-1.0)/ rho * (e-0.5*rho*u*u) > T_ignition){
     rhs += temp;
   }

#endif
  // auto var_vec_l = Variable_Vector_Isolator<grid_type>(solution_vector_m, gamma);
  // auto var_vec = Variable_Vector_Isolator<grid_type>(solution_vector, gamma);
  // auto var_vec_r = Variable_Vector_Isolator<grid_type>(solution_vector_p, gamma);
  //
  //  rhs += (hyperbolic_flux.flux(var_vec_l.w(), var_vec.w(), gamma) - hyperbolic_flux.flux(var_vec.w(), var_vec_r.w(), gamma)) /dx*dt;

  return rhs;
}

#endif //#ifndef VARIABLE_IMPLICIT_SCHEME_H
