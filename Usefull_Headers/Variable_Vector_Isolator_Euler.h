#ifndef VARIABLE_VECTOR_ISOLATOR_EULER_H
#define VARIABLE_VECTOR_ISOLATOR_EULER_H
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Header containg the definition of the class for a
///         Variable_Vector_Isolator_Euler.h
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
  /// \brief  Header containg an object to simplify returning isolated variables.
///////////////////////////////////////////////////////////////////////////////
template <typename solution_vector_type>
class Variable_Vector_Isolator_Euler{
public:
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  Variable_Vector_Isolator_Euler() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  Variable_Vector_Isolator_Euler(const Variable_Vector_Isolator_Euler&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  Variable_Vector_Isolator_Euler(Variable_Vector_Isolator_Euler&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  Variable_Vector_Isolator_Euler& operator=(const Variable_Vector_Isolator_Euler&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  Variable_Vector_Isolator_Euler& operator=(Variable_Vector_Isolator_Euler&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Constructor isolating variables.
  Variable_Vector_Isolator_Euler(solution_vector_type solution_vector, double gamma_in) : gamma(gamma_in){
    u_var = isolate_u(solution_vector);
    rho_var = isolate_rho(solution_vector);
    p_var = isolate_p(solution_vector);
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to return u.
  double u(){return u_var;};

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to return rho.
  double rho(){return rho_var;};

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to return p.
  double p(){return p_var;};

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  solution_vector_type w(){
    solution_vector_type temp; temp << rho(), u(), p();
    return temp;
  }

private:

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to isolate u.
  double isolate_u(solution_vector_type solution_vector) {
    return solution_vector[1] / solution_vector[0];
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to isolate rho.
  double isolate_rho(solution_vector_type solution_vector) {
    return solution_vector[0];
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to isolate T.
  double isolate_p(solution_vector_type solution_vector) {
    return (solution_vector[2] - solution_vector[1] *solution_vector[1] / solution_vector[0] / 2.0 )*(gamma-1.0);
  }

  double u_var;
  double rho_var;
  double p_var;
  double gamma;
};

#endif //#endif VARIABLE_VECTOR_ISOLATOR_EULER_H
