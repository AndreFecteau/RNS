#ifndef SOURCES_H
#define SOURCES_H
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Header containg the definition of the class for a
///         Sources.h
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
  /// \brief  Header containg the Sources solver for first order terms.
///////////////////////////////////////////////////////////////////////////////
#include <math.h>
#include <algorithm>
#include <iostream>
#include "../Usefull_Headers/Variable_Vector_Isolator.h"

template <typename global_solution_vector_type>
class Sources {
 public:
  using solution_vector_type = typename global_solution_vector_type::value_type;
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  Sources() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  Sources(const Sources&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  Sources(Sources&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  Sources& operator=(const Sources&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  Sources& operator=(Sources&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the flux of the given cells.
  /// \return Flux between cells.
  solution_vector_type flux(solution_vector_type solution_vector,
                            double gamma_in, double Q_in, double Lambda_in,
                            double theta_in) {
    gamma = gamma_in;
    Q = Q_in;
    Lambda = Lambda_in;
    theta = theta_in;
    solution_vector_var = Variable_Vector_Isolator<solution_vector_type>(solution_vector, gamma);
    solution_vector_type sources_vector;
    sources_vector << 0.0, 0.0, energy_source(), species_source();
    return sources_vector;
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  double energy_source() {
    return Q*Lambda*solution_vector_var.rho()*solution_vector_var.Y()*exp(-theta/solution_vector_var.T());
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  double species_source() {
    return -Lambda*solution_vector_var.rho()*solution_vector_var.Y()*exp(-theta/solution_vector_var.T());
  }

 private:
  Variable_Vector_Isolator<solution_vector_type> solution_vector_var;
  double gamma;
  double Q;
  double Lambda;
  double theta;
};

#endif //#ifndef SOURCES_H
