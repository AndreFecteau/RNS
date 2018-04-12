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
  Non_Dimensional_Navier_Stokes(scalar_type Pr_in, scalar_type Le_in, scalar_type Q_low_mach_in,
                                scalar_type theta_low_mach_in, scalar_type mf_in, scalar_type gamma_in,
                                scalar_type lambda_in, scalar_type T_ignition_scalar_in) :
                                Pr(Pr_in), Le(Le_in), Q_low_mach(Q_low_mach_in),
                                theta_low_mach(theta_low_mach_in), mf(mf_in),
                                gamma(gamma_in), lambda(lambda_in), T_ignition_scalar(T_ignition_scalar_in) {}

  scalar_type Pr;
  scalar_type Le;
  scalar_type Q_low_mach;
  scalar_type theta_low_mach;
  scalar_type mf;
  scalar_type gamma;
  scalar_type lambda;
  scalar_type T_ignition_scalar;


  scalar_type T_ignition() {return T_ignition_scalar/(gamma*mf*mf);}
  scalar_type Q() {return Q_low_mach/(mf*mf*(gamma - 1));}
  scalar_type theta() {return theta_low_mach/(mf*mf*gamma);}

  template<typename Archive>
  void serialize(Archive& archive) {
    archive(Pr, Le, Q_low_mach, theta_low_mach, mf, gamma, lambda, T_ignition_scalar);
  }

};


#endif //#ifndef NON_DIMENSIONAL_NAVIER_STOKES_H
