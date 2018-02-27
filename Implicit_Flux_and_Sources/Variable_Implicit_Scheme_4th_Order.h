#ifndef IMPLICIT_MATRIX_ENTRIES
#define IMPLICIT_MATRIX_ENTRIES

template <typename global_solution_vector_type, typename matrix_type>
class Implicit_Matrix_Entries {
using solution_vector_type = typename global_solution_vector_type::value_type;
 public:
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  Implicit_Matrix_Entries() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  Implicit_Matrix_Entries(const Implicit_Matrix_Entries&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  Implicit_Matrix_Entries(Implicit_Matrix_Entries&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  Implicit_Matrix_Entries& operator=(const Implicit_Matrix_Entries&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  Implicit_Matrix_Entries& operator=(Implicit_Matrix_Entries&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Constructor setting up required inputs.
  Implicit_Matrix_Entries(const solution_vector_type& solution_vector_mm,
                             const solution_vector_type& solution_vector_m,
                             const solution_vector_type& solution_vector,
                             const solution_vector_type& solution_vector_p,
                             const solution_vector_type& solution_vector_pp,
                             const solution_vector_type& DeltaUm,
                             const double& gamma_in, const double& Pr_in,
                             const double& Le_in, const double& Q_in, const double& lambda_in,
                             const double& theta_in, const double& dx_in, const double& dt_in,
                             const double& zeta_in, const double& Theta_in) :
                             gamma(gamma_in), Pr(Pr_in), Le(Le_in), Q(Q_in), lambda(lambda_in),
                             theta(theta_in), dx(dx_in), dt(dt_in), zeta(zeta_in),
                             Theta(Theta_in) {
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
  // pm = (solution_vector_m[2] - solution_vector_m[1]* solution_vector_m[1] /
  //      (2.0 * solution_vector_m[0])) * (gamma - 1);
  // p =  (solution_vector[2] - solution_vector[1]* solution_vector[1] /
  //      (2.0 * solution_vector[0])) * (gamma - 1);
  // pp = (solution_vector_p[2] - solution_vector_p[1]* solution_vector_p[1] /
  //      (2.0 * solution_vector_p[0])) * (gamma - 1);
  dUm1 = DeltaUm[0];
  dUm2 = DeltaUm[1];
  dUm3 = DeltaUm[2];
  dUm4 = DeltaUm[3];
  }

  matrix_type top2_matrix();
  matrix_type top_matrix();
  matrix_type mid_matrix();
  matrix_type bot_matrix();
  matrix_type bot2_matrix();
  solution_vector_type rhs_matrix();

 private:
  double rhomm, rhom, rho, rhop, rhopp;
  // double pmm, pm, p, pp, ppp;
  double umm, um, u, up, upp;
  double emm, em, e, ep, epp;
  double Ymm, Ym, Y, Yp, Ypp;
  char E;
  double gamma, Pr, Le, Q;
  double lambda, theta;
  double dx, dt;
  double zeta, Theta;
  double dUm1, dUm2, dUm3, dUm4;

  template <typename T>
  double Power(const T num, const int expo) {
    return pow(num, expo);
  }

  double Power(const char num, const double expo) {
    return exp(expo);
  }
};

template <typename global_solution_vector_type, typename matrix_type>
matrix_type Implicit_Matrix_Entries<global_solution_vector_type, matrix_type>::top2_matrix() {

  matrix_type top2;
  matrix_type temp;
///////////////////////////////////////////////////////////////////////////////
// Needed
///////////////////////////////////////////////////////////////////////////////
  top2 << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
///////////////////////////////////////////////////////////////////////////////
// Hyperbolic (Euler)
///////////////////////////////////////////////////////////////////////////////
#if defined(HYPERBOLIC)

  temp << 0,-(dt*Theta)/(12.*dx*(1 + zeta)),0,0,-(dt*(-3 + gamma)*Theta*Power(u,2))/(24.*dx*(1 + zeta)),(dt*(-3 + gamma)*Theta*u)/(12.*dx*(1 + zeta)),
   -(dt*(-1 + gamma)*Theta)/(12.*dx*(1 + zeta)),0,(dt*Theta*(e*gamma*u - (-1 + gamma)*rho*Power(u,3)))/(12.*dx*rho*(1 + zeta)),
   (dt*Theta*(-2*e*gamma + 3*(-1 + gamma)*rho*Power(u,2)))/(24.*dx*rho*(1 + zeta)),-(dt*gamma*Theta*u)/(12.*dx*(1 + zeta)),0,
   (dt*Theta*u*Y)/(12*dx + 12*dx*zeta),-(dt*Theta*Y)/(12.*dx*(1 + zeta)),0,-(dt*Theta*u)/(12.*dx*(1 + zeta));
  top2 += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Viscous (2nd order derivatives)
///////////////////////////////////////////////////////////////////////////////
#if defined(VISCOUS)

  temp << 0,0,0,0,-(dt*Pr*Theta*((8*rhom - rhomm - 8*rhop + rhopp)*u + rho*(6*u - 8*um + umm + 8*up - upp)))/(54.*Power(dx,2)*Power(rho,2)*(1 + zeta)),
   (dt*Pr*(6*rho + 8*rhom - rhomm - 8*rhop + rhopp)*Theta)/(54.*Power(dx,2)*Power(rho,2)*(1 + zeta)),0,0,
   (dt*Theta*(-6*e*gamma*(3*rho + 8*rhom - rhomm - 8*rhop + rhopp) +
        rho*(24*em*gamma - 3*emm*gamma - 24*ep*gamma + 3*epp*gamma + 18*gamma*rho*Power(u,2) - 24*Pr*rho*Power(u,2) + 24*gamma*rhom*Power(u,2) -
           32*Pr*rhom*Power(u,2) - 3*gamma*rhomm*Power(u,2) + 4*Pr*rhomm*Power(u,2) - 24*gamma*rhop*Power(u,2) + 32*Pr*rhop*Power(u,2) +
           3*gamma*rhopp*Power(u,2) - 4*Pr*rhopp*Power(u,2) - 48*gamma*rho*u*um + 64*Pr*rho*u*um + 6*gamma*rho*u*umm - 8*Pr*rho*u*umm +
           48*gamma*rho*u*up - 64*Pr*rho*u*up - 6*gamma*rho*u*upp + 8*Pr*rho*u*upp)))/(216.*Power(dx,2)*Power(rho,3)*(1 + zeta)),
   -(dt*(3*gamma - 4*Pr)*Theta*((8*rhom - rhomm - 8*rhop + rhopp)*u + rho*(6*u - 8*um + umm + 8*up - upp)))/(216.*Power(dx,2)*Power(rho,2)*(1 + zeta)),
   (dt*gamma*(6*rho + 8*rhom - rhomm - 8*rhop + rhopp)*Theta)/(72.*Power(dx,2)*Power(rho,2)*(1 + zeta)),0,
   -(dt*Theta*((8*rhom - rhomm - 8*rhop + rhopp)*Y + rho*(6*Y - 8*Ym + Ymm + 8*Yp - Ypp)))/(72.*Power(dx,2)*Le*Power(rho,2)*(1 + zeta)),0,0,
   (dt*(6*rho + 8*rhom - rhomm - 8*rhop + rhopp)*Theta)/(72.*Power(dx,2)*Le*Power(rho,2)*(1 + zeta));
  top2 += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Source
///////////////////////////////////////////////////////////////////////////////
#if defined(SOURCE)

  temp << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
  top2 += temp;

#endif

  return top2;
}


template <typename global_solution_vector_type, typename matrix_type>
matrix_type Implicit_Matrix_Entries<global_solution_vector_type, matrix_type>::top_matrix() {

  matrix_type top;
  matrix_type temp;
///////////////////////////////////////////////////////////////////////////////
// Needed
///////////////////////////////////////////////////////////////////////////////
  top << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
///////////////////////////////////////////////////////////////////////////////
// Hyperbolic (Euler)
///////////////////////////////////////////////////////////////////////////////
#if defined(HYPERBOLIC)

  temp << 0,(2*dt*Theta)/(3*dx + 3*dx*zeta),0,0,(dt*(-3 + gamma)*Theta*Power(u,2))/(3.*dx*(1 + zeta)),(-2*dt*(-3 + gamma)*Theta*u)/(3.*dx*(1 + zeta)),
   (2*dt*(-1 + gamma)*Theta)/(3.*dx*(1 + zeta)),0,(2*dt*Theta*u*(-(e*gamma) + (-1 + gamma)*rho*Power(u,2)))/(3.*dx*rho*(1 + zeta)),
   (dt*Theta*(2*e*gamma - 3*(-1 + gamma)*rho*Power(u,2)))/(3.*dx*rho*(1 + zeta)),(2*dt*gamma*Theta*u)/(3*dx + 3*dx*zeta),0,
   (-2*dt*Theta*u*Y)/(3*dx + 3*dx*zeta),(2*dt*Theta*Y)/(3*dx + 3*dx*zeta),0,(2*dt*Theta*u)/(3*dx + 3*dx*zeta);
  top += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Viscous (2nd order derivatives)
///////////////////////////////////////////////////////////////////////////////
#if defined(VISCOUS)

  temp << 0,0,0,0,(4*dt*Pr*Theta*((8*rhom - rhomm - 8*rhop + rhopp)*u + rho*(12*u - 8*um + umm + 8*up - upp)))/(27.*Power(dx,2)*Power(rho,2)*(1 + zeta)),
   (-4*dt*Pr*(12*rho + 8*rhom - rhomm - 8*rhop + rhopp)*Theta)/(27.*Power(dx,2)*Power(rho,2)*(1 + zeta)),0,0,
   (dt*Theta*(6*e*gamma*(6*rho + 8*rhom - rhomm - 8*rhop + rhopp) +
        rho*(-24*em*gamma + 3*emm*gamma + 24*ep*gamma - 3*epp*gamma - 36*gamma*rho*Power(u,2) + 48*Pr*rho*Power(u,2) - 24*gamma*rhom*Power(u,2) +
           32*Pr*rhom*Power(u,2) + 3*gamma*rhomm*Power(u,2) - 4*Pr*rhomm*Power(u,2) + 24*gamma*rhop*Power(u,2) - 32*Pr*rhop*Power(u,2) -
           3*gamma*rhopp*Power(u,2) + 4*Pr*rhopp*Power(u,2) + 48*gamma*rho*u*um - 64*Pr*rho*u*um - 6*gamma*rho*u*umm + 8*Pr*rho*u*umm -
           48*gamma*rho*u*up + 64*Pr*rho*u*up + 6*gamma*rho*u*upp - 8*Pr*rho*u*upp)))/(27.*Power(dx,2)*Power(rho,3)*(1 + zeta)),
   (dt*(3*gamma - 4*Pr)*Theta*((8*rhom - rhomm - 8*rhop + rhopp)*u + rho*(12*u - 8*um + umm + 8*up - upp)))/(27.*Power(dx,2)*Power(rho,2)*(1 + zeta)),
   -(dt*gamma*(12*rho + 8*rhom - rhomm - 8*rhop + rhopp)*Theta)/(9.*Power(dx,2)*Power(rho,2)*(1 + zeta)),0,
   (dt*Theta*((8*rhom - rhomm - 8*rhop + rhopp)*Y + rho*(12*Y - 8*Ym + Ymm + 8*Yp - Ypp)))/(9.*Power(dx,2)*Le*Power(rho,2)*(1 + zeta)),0,0,
   -(dt*(12*rho + 8*rhom - rhomm - 8*rhop + rhopp)*Theta)/(9.*Power(dx,2)*Le*Power(rho,2)*(1 + zeta));
  top += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Source
///////////////////////////////////////////////////////////////////////////////
#if defined(SOURCE)

  temp << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
  top += temp;

#endif

  return top;
}

template <typename global_solution_vector_type, typename matrix_type>
matrix_type Implicit_Matrix_Entries<global_solution_vector_type, matrix_type>::mid_matrix() {

  matrix_type mid;
  matrix_type temp;
///////////////////////////////////////////////////////////////////////////////
// Needed
///////////////////////////////////////////////////////////////////////////////
  mid << 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1;
///////////////////////////////////////////////////////////////////////////////
// Hyperbolic (Euler)
///////////////////////////////////////////////////////////////////////////////
#if defined(HYPERBOLIC)

  temp << 0,0,0,0,-(dt*(-3 + gamma)*Theta*u*(8*um - umm - 8*up + upp))/(12.*dx*(1 + zeta)),
   (dt*(-3 + gamma)*Theta*(8*um - umm - 8*up + upp))/(12.*dx*(1 + zeta)),0,0,
   (dt*Theta*(8*em*gamma*rho*u - emm*gamma*rho*u - 8*ep*gamma*rho*u + epp*gamma*rho*u - 8*e*gamma*rhom*u + e*gamma*rhomm*u + 8*e*gamma*rhop*u -
        e*gamma*rhopp*u + 8*e*gamma*rho*um + 24*Power(rho,2)*Power(u,2)*um - 24*gamma*Power(rho,2)*Power(u,2)*um - e*gamma*rho*umm -
        3*Power(rho,2)*Power(u,2)*umm + 3*gamma*Power(rho,2)*Power(u,2)*umm - 8*e*gamma*rho*up - 24*Power(rho,2)*Power(u,2)*up +
        24*gamma*Power(rho,2)*Power(u,2)*up + e*gamma*rho*upp + 3*Power(rho,2)*Power(u,2)*upp - 3*gamma*Power(rho,2)*Power(u,2)*upp))/
    (12.*dx*Power(rho,2)*(1 + zeta)),(dt*Theta*(-8*em*gamma*rho + emm*gamma*rho + 8*ep*gamma*rho - epp*gamma*rho + 8*e*gamma*rhom - e*gamma*rhomm -
        8*e*gamma*rhop + e*gamma*rhopp - 24*Power(rho,2)*u*um + 24*gamma*Power(rho,2)*u*um + 3*Power(rho,2)*u*umm - 3*gamma*Power(rho,2)*u*umm +
        24*Power(rho,2)*u*up - 24*gamma*Power(rho,2)*u*up + 3*(-1 + gamma)*Power(rho,2)*u*upp))/(12.*dx*Power(rho,2)*(1 + zeta)),
   (dt*gamma*Theta*(-8*um + umm + 8*up - upp))/(12.*dx*(1 + zeta)),0,
   (dt*Theta*((8*um - umm - 8*up + upp)*Y + u*(8*Ym - Ymm - 8*Yp + Ypp)))/(12.*dx*(1 + zeta)),(dt*Theta*(-8*Ym + Ymm + 8*Yp - Ypp))/(12.*dx*(1 + zeta)),
   0,(dt*Theta*(-8*um + umm + 8*up - upp))/(12.*dx*(1 + zeta));
  mid += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Viscous (2nd order derivatives)
///////////////////////////////////////////////////////////////////////////////
#if defined(VISCOUS)

  temp << 0,0,0,0,(dt*Pr*Theta*(Power(8*rhom - rhomm - 8*rhop + rhopp,2)*u - 6*Power(rho,2)*(30*u - 16*um + umm - 16*up + upp) +
        rho*(-96*rhop*u + 6*rhopp*u + 64*rhop*um - 8*rhopp*um - 8*rhop*umm + rhopp*umm - 64*rhop*up + 8*rhopp*up + 8*rhop*upp - rhopp*upp +
           rhomm*(6*u + 8*um - umm - 8*up + upp) - 8*rhom*(12*u + 8*um - umm - 8*up + upp))))/(54.*Power(dx,2)*Power(rho,3)*(1 + zeta)),
   (dt*Pr*(96*rho*rhom - 6*rho*(rhomm - 16*rhop + rhopp) - Power(8*rhom - rhomm - 8*rhop + rhopp,2))*Theta)/(54.*Power(dx,2)*Power(rho,3)*(1 + zeta)),0,
   0,(dt*Theta*(288*ep*gamma*Power(rho,2) - 18*epp*gamma*Power(rho,2) - 576*e*gamma*rho*rhom + 384*ep*gamma*rho*rhom - 48*epp*gamma*rho*rhom +
        576*e*gamma*Power(rhom,2) + 36*e*gamma*rho*rhomm - 48*ep*gamma*rho*rhomm + 6*epp*gamma*rho*rhomm - 144*e*gamma*rhom*rhomm +
        9*e*gamma*Power(rhomm,2) - 576*e*gamma*rho*rhop - 384*ep*gamma*rho*rhop + 48*epp*gamma*rho*rhop - 1152*e*gamma*rhom*rhop +
        144*e*gamma*rhomm*rhop + 576*e*gamma*Power(rhop,2) - 6*emm*gamma*rho*(3*rho - 8*rhom + rhomm + 8*rhop - rhopp) +
        48*em*gamma*rho*(6*rho - 8*rhom + rhomm + 8*rhop - rhopp) + 36*e*gamma*rho*rhopp + 48*ep*gamma*rho*rhopp - 6*epp*gamma*rho*rhopp +
        144*e*gamma*rhom*rhopp - 18*e*gamma*rhomm*rhopp - 144*e*gamma*rhop*rhopp + 9*e*gamma*Power(rhopp,2) + 1080*gamma*Power(rho,3)*Power(u,2) -
        1440*Pr*Power(rho,3)*Power(u,2) + 288*gamma*Power(rho,2)*rhom*Power(u,2) - 384*Pr*Power(rho,2)*rhom*Power(u,2) -
        192*gamma*rho*Power(rhom,2)*Power(u,2) + 256*Pr*rho*Power(rhom,2)*Power(u,2) - 18*gamma*Power(rho,2)*rhomm*Power(u,2) +
        24*Pr*Power(rho,2)*rhomm*Power(u,2) + 48*gamma*rho*rhom*rhomm*Power(u,2) - 64*Pr*rho*rhom*rhomm*Power(u,2) -
        3*gamma*rho*Power(rhomm,2)*Power(u,2) + 4*Pr*rho*Power(rhomm,2)*Power(u,2) + 288*gamma*Power(rho,2)*rhop*Power(u,2) -
        384*Pr*Power(rho,2)*rhop*Power(u,2) + 384*gamma*rho*rhom*rhop*Power(u,2) - 512*Pr*rho*rhom*rhop*Power(u,2) -
        48*gamma*rho*rhomm*rhop*Power(u,2) + 64*Pr*rho*rhomm*rhop*Power(u,2) - 192*gamma*rho*Power(rhop,2)*Power(u,2) +
        256*Pr*rho*Power(rhop,2)*Power(u,2) - 18*gamma*Power(rho,2)*rhopp*Power(u,2) + 24*Pr*Power(rho,2)*rhopp*Power(u,2) -
        48*gamma*rho*rhom*rhopp*Power(u,2) + 64*Pr*rho*rhom*rhopp*Power(u,2) + 6*gamma*rho*rhomm*rhopp*Power(u,2) - 8*Pr*rho*rhomm*rhopp*Power(u,2) +
        48*gamma*rho*rhop*rhopp*Power(u,2) - 64*Pr*rho*rhop*rhopp*Power(u,2) - 3*gamma*rho*Power(rhopp,2)*Power(u,2) +
        4*Pr*rho*Power(rhopp,2)*Power(u,2) - 576*gamma*Power(rho,3)*u*um + 768*Pr*Power(rho,3)*u*um + 384*gamma*Power(rho,2)*rhom*u*um -
        512*Pr*Power(rho,2)*rhom*u*um - 48*gamma*Power(rho,2)*rhomm*u*um + 64*Pr*Power(rho,2)*rhomm*u*um - 384*gamma*Power(rho,2)*rhop*u*um +
        512*Pr*Power(rho,2)*rhop*u*um + 48*gamma*Power(rho,2)*rhopp*u*um - 64*Pr*Power(rho,2)*rhopp*u*um - 192*gamma*Power(rho,3)*Power(um,2) +
        256*Pr*Power(rho,3)*Power(um,2) + 36*gamma*Power(rho,3)*u*umm - 48*Pr*Power(rho,3)*u*umm - 48*gamma*Power(rho,2)*rhom*u*umm +
        64*Pr*Power(rho,2)*rhom*u*umm + 6*gamma*Power(rho,2)*rhomm*u*umm - 8*Pr*Power(rho,2)*rhomm*u*umm + 48*gamma*Power(rho,2)*rhop*u*umm -
        64*Pr*Power(rho,2)*rhop*u*umm - 6*gamma*Power(rho,2)*rhopp*u*umm + 8*Pr*Power(rho,2)*rhopp*u*umm + 48*gamma*Power(rho,3)*um*umm -
        64*Pr*Power(rho,3)*um*umm - 3*gamma*Power(rho,3)*Power(umm,2) + 4*Pr*Power(rho,3)*Power(umm,2) - 576*gamma*Power(rho,3)*u*up +
        768*Pr*Power(rho,3)*u*up - 384*gamma*Power(rho,2)*rhom*u*up + 512*Pr*Power(rho,2)*rhom*u*up + 48*gamma*Power(rho,2)*rhomm*u*up -
        64*Pr*Power(rho,2)*rhomm*u*up + 384*gamma*Power(rho,2)*rhop*u*up - 512*Pr*Power(rho,2)*rhop*u*up - 48*gamma*Power(rho,2)*rhopp*u*up +
        64*Pr*Power(rho,2)*rhopp*u*up + 384*gamma*Power(rho,3)*um*up - 512*Pr*Power(rho,3)*um*up - 48*gamma*Power(rho,3)*umm*up +
        64*Pr*Power(rho,3)*umm*up - 192*gamma*Power(rho,3)*Power(up,2) + 256*Pr*Power(rho,3)*Power(up,2) +
        2*(3*gamma - 4*Pr)*Power(rho,2)*((8*rhom - rhomm - 8*rhop + rhopp)*u + rho*(6*u - 8*um + umm + 8*up))*upp +
        (-3*gamma + 4*Pr)*Power(rho,3)*Power(upp,2)))/(216.*Power(dx,2)*Power(rho,4)*(1 + zeta)),
   -(dt*(3*gamma - 4*Pr)*Theta*(-(Power(8*rhom - rhomm - 8*rhop + rhopp,2)*u) + 6*Power(rho,2)*(30*u - 16*um + umm - 16*up + upp) +
         rho*(96*rhop*u - 6*rhopp*u - 64*rhop*um + 8*rhopp*um + 8*rhop*umm - rhopp*umm + 64*rhop*up - 8*rhopp*up +
            rhomm*(-6*u - 8*um + umm + 8*up - upp) - 8*rhop*upp + rhopp*upp + 8*rhom*(12*u + 8*um - umm - 8*up + upp))))/
    (216.*Power(dx,2)*Power(rho,3)*(1 + zeta)),(dt*gamma*(96*rho*rhom - 6*rho*(rhomm - 16*rhop + rhopp) - Power(8*rhom - rhomm - 8*rhop + rhopp,2))*
      Theta)/(72.*Power(dx,2)*Power(rho,3)*(1 + zeta)),0,(dt*Theta*
      (Power(8*rhom - rhomm - 8*rhop + rhopp,2)*Y - 6*Power(rho,2)*(30*Y - 16*Ym + Ymm - 16*Yp + Ypp) +
        rho*(-96*rhop*Y + 6*rhopp*Y + 64*rhop*Ym - 8*rhopp*Ym - 8*rhop*Ymm + rhopp*Ymm - 64*rhop*Yp + 8*rhopp*Yp + 8*rhop*Ypp - rhopp*Ypp +
           rhomm*(6*Y + 8*Ym - Ymm - 8*Yp + Ypp) - 8*rhom*(12*Y + 8*Ym - Ymm - 8*Yp + Ypp))))/(72.*Power(dx,2)*Le*Power(rho,3)*(1 + zeta)),0,0,
   (dt*(96*rho*rhom - 6*rho*(rhomm - 16*rhop + rhopp) - Power(8*rhom - rhomm - 8*rhop + rhopp,2))*Theta)/(72.*Power(dx,2)*Le*Power(rho,3)*(1 + zeta));
  mid += temp;

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

   if((gamma-1.0)/ rho * (e-0.5*rho*u*u) > 1.1/(gamma*0.005*0.005)){
    mid += temp;
  }

#endif

  return mid;
}

template <typename global_solution_vector_type, typename matrix_type>
matrix_type Implicit_Matrix_Entries<global_solution_vector_type, matrix_type>::bot_matrix() {

  matrix_type bot;
  matrix_type temp;
///////////////////////////////////////////////////////////////////////////////
// Needed
///////////////////////////////////////////////////////////////////////////////
  bot << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
///////////////////////////////////////////////////////////////////////////////
// Hyperbolic (Euler)
///////////////////////////////////////////////////////////////////////////////
#if defined(HYPERBOLIC)

  temp << 0,(-2*dt*Theta)/(3*dx + 3*dx*zeta),0,0,-(dt*(-3 + gamma)*Theta*Power(u,2))/(3.*dx*(1 + zeta)),(2*dt*(-3 + gamma)*Theta*u)/(3.*dx*(1 + zeta)),
   (-2*dt*(-1 + gamma)*Theta)/(3.*dx*(1 + zeta)),0,(2*dt*Theta*(e*gamma*u - (-1 + gamma)*rho*Power(u,3)))/(3.*dx*rho*(1 + zeta)),
   (dt*Theta*(-2*e*gamma + 3*(-1 + gamma)*rho*Power(u,2)))/(3.*dx*rho*(1 + zeta)),(-2*dt*gamma*Theta*u)/(3*dx + 3*dx*zeta),0,
   (2*dt*Theta*u*Y)/(3*dx + 3*dx*zeta),(-2*dt*Theta*Y)/(3*dx + 3*dx*zeta),0,(-2*dt*Theta*u)/(3*dx + 3*dx*zeta);
  bot += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Viscous (2nd order derivatives)
///////////////////////////////////////////////////////////////////////////////
#if defined(VISCOUS)

  temp << 0,0,0,0,(4*dt*Pr*Theta*((-8*rhom + rhomm + 8*rhop - rhopp)*u + rho*(12*u + 8*um - umm - 8*up + upp)))/(27.*Power(dx,2)*Power(rho,2)*(1 + zeta)),
   (-4*dt*Pr*(12*rho - 8*rhom + rhomm + 8*rhop - rhopp)*Theta)/(27.*Power(dx,2)*Power(rho,2)*(1 + zeta)),0,0,
   (dt*Theta*(6*e*gamma*(6*rho - 8*rhom + rhomm + 8*rhop - rhopp) +
        rho*(24*em*gamma - 3*emm*gamma - 24*ep*gamma + 3*epp*gamma - 36*gamma*rho*Power(u,2) + 48*Pr*rho*Power(u,2) + 24*gamma*rhom*Power(u,2) -
           32*Pr*rhom*Power(u,2) - 3*gamma*rhomm*Power(u,2) + 4*Pr*rhomm*Power(u,2) - 24*gamma*rhop*Power(u,2) + 32*Pr*rhop*Power(u,2) +
           3*gamma*rhopp*Power(u,2) - 4*Pr*rhopp*Power(u,2) - 48*gamma*rho*u*um + 64*Pr*rho*u*um + 6*gamma*rho*u*umm - 8*Pr*rho*u*umm +
           48*gamma*rho*u*up - 64*Pr*rho*u*up - 6*gamma*rho*u*upp + 8*Pr*rho*u*upp)))/(27.*Power(dx,2)*Power(rho,3)*(1 + zeta)),
   (dt*(3*gamma - 4*Pr)*Theta*((-8*rhom + rhomm + 8*rhop - rhopp)*u + rho*(12*u + 8*um - umm - 8*up + upp)))/(27.*Power(dx,2)*Power(rho,2)*(1 + zeta)),
   -(dt*gamma*(12*rho - 8*rhom + rhomm + 8*rhop - rhopp)*Theta)/(9.*Power(dx,2)*Power(rho,2)*(1 + zeta)),0,
   (dt*Theta*((-8*rhom + rhomm + 8*rhop - rhopp)*Y + rho*(12*Y + 8*Ym - Ymm - 8*Yp + Ypp)))/(9.*Power(dx,2)*Le*Power(rho,2)*(1 + zeta)),0,0,
   -(dt*(12*rho - 8*rhom + rhomm + 8*rhop - rhopp)*Theta)/(9.*Power(dx,2)*Le*Power(rho,2)*(1 + zeta));
  bot += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Source
///////////////////////////////////////////////////////////////////////////////
#if defined(SOURCE)

  temp << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
  bot += temp;

#endif
  return bot;
}

template <typename global_solution_vector_type, typename matrix_type>
matrix_type Implicit_Matrix_Entries<global_solution_vector_type, matrix_type>::bot2_matrix() {

  matrix_type bot2;
  matrix_type temp;
///////////////////////////////////////////////////////////////////////////////
// Needed
///////////////////////////////////////////////////////////////////////////////
  bot2 << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
///////////////////////////////////////////////////////////////////////////////
// Hyperbolic (Euler)
///////////////////////////////////////////////////////////////////////////////
#if defined(HYPERBOLIC)

  temp << 0,(dt*Theta)/(12*dx + 12*dx*zeta),0,0,(dt*(-3 + gamma)*Theta*Power(u,2))/(24.*dx*(1 + zeta)),-(dt*(-3 + gamma)*Theta*u)/(12.*dx*(1 + zeta)),
   (dt*(-1 + gamma)*Theta)/(12.*dx*(1 + zeta)),0,(dt*Theta*u*(-(e*gamma) + (-1 + gamma)*rho*Power(u,2)))/(12.*dx*rho*(1 + zeta)),
   (dt*Theta*(2*e*gamma - 3*(-1 + gamma)*rho*Power(u,2)))/(24.*dx*rho*(1 + zeta)),(dt*gamma*Theta*u)/(12*dx + 12*dx*zeta),0,
   -(dt*Theta*u*Y)/(12.*dx*(1 + zeta)),(dt*Theta*Y)/(12*dx + 12*dx*zeta),0,(dt*Theta*u)/(12*dx + 12*dx*zeta);
  bot2 += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Viscous (2nd order derivatives)
///////////////////////////////////////////////////////////////////////////////
#if defined(VISCOUS)

  temp << 0,0,0,0,-(dt*Pr*Theta*((-8*rhom + rhomm + 8*rhop - rhopp)*u + rho*(6*u + 8*um - umm - 8*up + upp)))/(54.*Power(dx,2)*Power(rho,2)*(1 + zeta)),
   (dt*Pr*(6*rho - 8*rhom + rhomm + 8*rhop - rhopp)*Theta)/(54.*Power(dx,2)*Power(rho,2)*(1 + zeta)),0,0,
   (dt*Theta*(-6*e*gamma*(3*rho - 8*rhom + rhomm + 8*rhop - rhopp) +
        rho*(-24*em*gamma + 3*emm*gamma + 24*ep*gamma - 3*epp*gamma + 18*gamma*rho*Power(u,2) - 24*Pr*rho*Power(u,2) - 24*gamma*rhom*Power(u,2) +
           32*Pr*rhom*Power(u,2) + 3*gamma*rhomm*Power(u,2) - 4*Pr*rhomm*Power(u,2) + 24*gamma*rhop*Power(u,2) - 32*Pr*rhop*Power(u,2) -
           3*gamma*rhopp*Power(u,2) + 4*Pr*rhopp*Power(u,2) + 48*gamma*rho*u*um - 64*Pr*rho*u*um - 6*gamma*rho*u*umm + 8*Pr*rho*u*umm -
           48*gamma*rho*u*up + 64*Pr*rho*u*up + 6*gamma*rho*u*upp - 8*Pr*rho*u*upp)))/(216.*Power(dx,2)*Power(rho,3)*(1 + zeta)),
   -(dt*(3*gamma - 4*Pr)*Theta*((-8*rhom + rhomm + 8*rhop - rhopp)*u + rho*(6*u + 8*um - umm - 8*up + upp)))/(216.*Power(dx,2)*Power(rho,2)*(1 + zeta)),
   (dt*gamma*(6*rho - 8*rhom + rhomm + 8*rhop - rhopp)*Theta)/(72.*Power(dx,2)*Power(rho,2)*(1 + zeta)),0,
   -(dt*Theta*((-8*rhom + rhomm + 8*rhop - rhopp)*Y + rho*(6*Y + 8*Ym - Ymm - 8*Yp + Ypp)))/(72.*Power(dx,2)*Le*Power(rho,2)*(1 + zeta)),0,0,
   (dt*(6*rho - 8*rhom + rhomm + 8*rhop - rhopp)*Theta)/(72.*Power(dx,2)*Le*Power(rho,2)*(1 + zeta));
  bot2 += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Source
///////////////////////////////////////////////////////////////////////////////
#if defined(SOURCE)

  temp << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
  bot2 += temp;

#endif
  return bot2;
}


template <typename global_solution_vector_type, typename matrix_type>
typename global_solution_vector_type::value_type
Implicit_Matrix_Entries<global_solution_vector_type, matrix_type>::rhs_matrix() {

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

  temp << (dt*(-((((-2*rhom)/3. + rhomm/12. + (2*rhop)/3. - rhopp/12.)*u)/dx) - (rho*((-2*um)/3. + umm/12. + (2*up)/3. - upp/12.))/dx))/(1 + zeta),
   (dt*(-((((-2*rhom)/3. + rhomm/12. + (2*rhop)/3. - rhopp/12.)*Power(u,2))/dx) -
        (-1 + gamma)*(((-2*em)/3. + emm/12. + (2*ep)/3. - epp/12.)/dx - (((-2*rhom)/3. + rhomm/12. + (2*rhop)/3. - rhopp/12.)*Power(u,2))/(2.*dx) -
           (rho*u*((-2*um)/3. + umm/12. + (2*up)/3. - upp/12.))/dx) - (2*rho*u*((-2*um)/3. + umm/12. + (2*up)/3. - upp/12.))/dx))/(1 + zeta),
   (dt*(-(u*(((-2*em)/3. + emm/12. + (2*ep)/3. - epp/12.)/dx +
             (-1 + gamma)*(((-2*em)/3. + emm/12. + (2*ep)/3. - epp/12.)/dx - (((-2*rhom)/3. + rhomm/12. + (2*rhop)/3. - rhopp/12.)*Power(u,2))/(2.*dx) -
                (rho*u*((-2*um)/3. + umm/12. + (2*up)/3. - upp/12.))/dx))) -
        ((e + (-1 + gamma)*(e - (rho*Power(u,2))/2.))*((-2*um)/3. + umm/12. + (2*up)/3. - upp/12.))/dx))/(1 + zeta),
   (dt*(-((((-2*rhom)/3. + rhomm/12. + (2*rhop)/3. - rhopp/12.)*u*Y)/dx) - (rho*((-2*um)/3. + umm/12. + (2*up)/3. - upp/12.)*Y)/dx -
        (rho*u*((-2*Ym)/3. + Ymm/12. + (2*Yp)/3. - Ypp/12.))/dx))/(1 + zeta);
  rhs += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Viscous (2nd order derivatives)
///////////////////////////////////////////////////////////////////////////////
#if defined(VISCOUS)

  temp << 0,(4*dt*Pr*((-5*u)/2. + (4*um)/3. - umm/12. + (4*up)/3. - upp/12.))/(3.*Power(dx,2)*(1 + zeta)),
   (dt*((gamma*((2*(-1 + gamma)*Power((-2*rhom)/3. + rhomm/12. + (2*rhop)/3. - rhopp/12.,2)*(e - (rho*Power(u,2))/2.))/(Power(dx,2)*Power(rho,3)) -
             ((-1 + gamma)*((-5*rho)/2. + (4*rhom)/3. - rhomm/12. + (4*rhop)/3. - rhopp/12.)*(e - (rho*Power(u,2))/2.))/(Power(dx,2)*Power(rho,2)) -
             (2*(-1 + gamma)*((-2*rhom)/3. + rhomm/12. + (2*rhop)/3. - rhopp/12.)*
                (((-2*em)/3. + emm/12. + (2*ep)/3. - epp/12.)/dx - (((-2*rhom)/3. + rhomm/12. + (2*rhop)/3. - rhopp/12.)*Power(u,2))/(2.*dx) -
                  (rho*u*((-2*um)/3. + umm/12. + (2*up)/3. - upp/12.))/dx))/(dx*Power(rho,2)) +
             ((-1 + gamma)*(((-5*e)/2. + (4*em)/3. - emm/12. + (4*ep)/3. - epp/12.)/Power(dx,2) -
                  (((-5*rho)/2. + (4*rhom)/3. - rhomm/12. + (4*rhop)/3. - rhopp/12.)*Power(u,2))/(2.*Power(dx,2)) -
                  (2*((-2*rhom)/3. + rhomm/12. + (2*rhop)/3. - rhopp/12.)*u*((-2*um)/3. + umm/12. + (2*up)/3. - upp/12.))/Power(dx,2) -
                  (rho*Power((-2*um)/3. + umm/12. + (2*up)/3. - upp/12.,2))/Power(dx,2) -
                  (rho*u*((-5*u)/2. + (4*um)/3. - umm/12. + (4*up)/3. - upp/12.))/Power(dx,2)))/rho))/(-1 + gamma) +
        (4*Pr*Power((-2*um)/3. + umm/12. + (2*up)/3. - upp/12.,2))/(3.*Power(dx,2)) +
        (4*Pr*u*((-5*u)/2. + (4*um)/3. - umm/12. + (4*up)/3. - upp/12.))/(3.*Power(dx,2))))/(1 + zeta),
   (dt*((-5*Y)/2. + (4*Ym)/3. - Ymm/12. + (4*Yp)/3. - Ypp/12.))/(Power(dx,2)*Le*(1 + zeta));
  rhs += temp;

#endif
///////////////////////////////////////////////////////////////////////////////
// Source
///////////////////////////////////////////////////////////////////////////////
#if defined(SOURCE)

  temp << 0.,0.,(dt*lambda*Q*rho*Y)/(Power(E,(rho*theta)/((-1 + gamma)*(e - (rho*Power(u,2))/2.)))*(1 + zeta)),
   -((dt*lambda*rho*Y)/(Power(E,(rho*theta)/((-1 + gamma)*(e - (rho*Power(u,2))/2.)))*(1 + zeta)));

   if((gamma-1.0)/ rho * (e-0.5*rho*u*u) > 1.1/(gamma*0.005*0.005)){
     rhs += temp;
   }

#endif
  // auto var_vec_l = Variable_Vector_Isolator<solution_vector_type>(solution_vector_m, gamma);
  // auto var_vec = Variable_Vector_Isolator<solution_vector_type>(solution_vector, gamma);
  // auto var_vec_r = Variable_Vector_Isolator<solution_vector_type>(solution_vector_p, gamma);
  //
  //  rhs += (hyperbolic_flux.flux(var_vec_l.w(), var_vec.w(), gamma) - hyperbolic_flux.flux(var_vec.w(), var_vec_r.w(), gamma)) /dx*dt;

  return rhs;
}

#endif //#ifndef IMPLICIT_MATRIX_ENTRIES_CD
