#ifndef CENTERED_DIFFERENCE_H
#define CENTERED_DIFFERENCE_H
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Header containg the definition of the class for a
///         Centered_Difference.h
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
  /// \brief  Header containg the Centered_Difference solver for first
  ///         order terms.
///////////////////////////////////////////////////////////////////////////////
#include <math.h>
#include <algorithm>
#include <iostream>
// #include "Eigen/Dense"
#include "../Usefull_Headers/Variable_Vector_Isolator.h"

template <typename global_solution_vector_type>
class Centered_Difference {
 public:
   using solution_vector_type = typename global_solution_vector_type::value_type;
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  Centered_Difference() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  Centered_Difference(const Centered_Difference&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  Centered_Difference(Centered_Difference&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  Centered_Difference& operator=(const Centered_Difference&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  Centered_Difference& operator=(Centered_Difference&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the flux of the given cells.
  /// \return Flux of a cell.
  solution_vector_type flux(solution_vector_type solution_vector_min1,
                            solution_vector_type solution_vector_centered,
                            solution_vector_type solution_vector_plus1,
                            double gamma_in, double Le_in, double Pr_in, double dx_in) {
    dx = dx_in;
    Pr = Pr_in;
    Le = Le_in;
    gamma = gamma_in;
    min1_var      = Variable_Vector_Isolator<solution_vector_type>(solution_vector_min1, gamma);
    centered_var  = Variable_Vector_Isolator<solution_vector_type>(solution_vector_centered, gamma);
    plus1_var     = Variable_Vector_Isolator<solution_vector_type>(solution_vector_plus1, gamma);

    solution_vector_type flux_vector;
    flux_vector[0] = 0.0;
    flux_vector[1] = momentum_cons_centered_difference();
    flux_vector[2] = energy_centered_difference();
    flux_vector[3] = species_cons_centered_difference();
    return flux_vector;
  }

private:

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the the mass conservation higher
  ///        order terms
  double momentum_cons_centered_difference() {
    return Pr * 4.0 / 3.0 * (plus1_var.u() - 2.0 * centered_var.u() +
           min1_var.u()) / (dx * dx);
  }

  double energy_centered_difference() {
    return gamma / (gamma - 1.0) * (plus1_var.T() - 2.0 * centered_var.T() +
           min1_var.T()) / (dx * dx) +
           4.0 / 3.0 * Pr * (plus1_var.u() * plus1_var.u() -
           2.0 * centered_var.u() * centered_var.u() +
           min1_var.u() * min1_var.u()) * 0.5 / (dx * dx);
  }


  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to calculate the the species conservation higher
  ///        order terms
  double species_cons_centered_difference() {
    return 1.0 / Le * (plus1_var.Y() - 2.0 * centered_var.Y() + min1_var.Y()) / (dx * dx);
  }

  Variable_Vector_Isolator<solution_vector_type> min1_var;
  Variable_Vector_Isolator<solution_vector_type> centered_var;
  Variable_Vector_Isolator<solution_vector_type> plus1_var;
  double gamma;
  double Le;
  double Pr;
  double dx;
};

#endif //#ifndef CENTERED_DIFFERENCE_H
