#ifndef NON_DIMENSIONAL_NAVIER_STOKES_H
#define NON_DIMENSIONAL_NAVIER_STOKES_H

template <typename scalar_type_in>
struct Non_Dimensional_Navier_Stokes {

 public:

  using scalar_type = scalar_type_in;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  Non_Dimensional_Navier_Stokes() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  Non_Dimensional_Navier_Stokes(const Non_Dimensional_Navier_Stokes&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  Non_Dimensional_Navier_Stokes(Non_Dimensional_Navier_Stokes&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  Non_Dimensional_Navier_Stokes& operator=(const Non_Dimensional_Navier_Stokes&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  Non_Dimensional_Navier_Stokes& operator=(Non_Dimensional_Navier_Stokes&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Constructor setting up required inputs.
  Non_Dimensional_Navier_Stokes(scalar_type Pr_in, scalar_type Le_in, scalar_type Q_in,
                                scalar_type theta_in, scalar_type mf_in, scalar_type gamma_in,
                                scalar_type lambda_in, scalar_type T_ignition_in) :
                                Pr(Pr_in), Le(Le_in), Q(Q_in), theta(theta_in), mf(mf_in),
                                gamma(gamma_in), lambda(lambda_in), T_ignition(T_ignition_in) {}

  scalar_type Pr;
  scalar_type Le;
  scalar_type Q;
  scalar_type theta;
  scalar_type mf;
  scalar_type gamma;
  scalar_type lambda;
  scalar_type T_ignition;

  scalar_type Q_low_mach() {return Q*mf*mf*(gamma - 1);}
  scalar_type theta_low_mach() {return theta*mf*mf*gamma;}

  template<typename Archive>
  void serialize(Archive& archive) {
    archive(Pr, Le, Q, theta, mf, gamma);
  }

};


#endif //#ifndef NON_DIMENSIONAL_NAVIER_STOKES_H
