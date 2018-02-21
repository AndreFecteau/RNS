#ifndef VARIABLE_VECTOR_ISOLATOR_H
#define VARIABLE_VECTOR_ISOLATOR_H
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Header containg the definition of the class for a
///         Variable_Vector_Isolator.h
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
  /// \brief  Header containg an object to simplify returning isolated variables.
///////////////////////////////////////////////////////////////////////////////
template <typename solution_vector_type>
class Variable_Vector_Isolator{
public:
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  Variable_Vector_Isolator() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  Variable_Vector_Isolator(const Variable_Vector_Isolator&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  Variable_Vector_Isolator(Variable_Vector_Isolator&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  Variable_Vector_Isolator& operator=(const Variable_Vector_Isolator&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  Variable_Vector_Isolator& operator=(Variable_Vector_Isolator&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Constructor isolating variables.
  Variable_Vector_Isolator(solution_vector_type solution_vector, double gamma_in) : gamma(gamma_in){
    u_var = isolate_u(solution_vector);
    T_var = isolate_T(solution_vector);
    rho_var = isolate_rho(solution_vector);
    Y_var = isolate_Y(solution_vector);
  }
  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to return T.
  double T(){return T_var;};

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to return u.
  double u(){return u_var;};

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to return rho.
  double rho(){return rho_var;};

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to return Y.
  double Y(){return Y_var;};

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to return p.
  double p(){return T_var * rho_var;};

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to return primitive variable vector.
  solution_vector_type w(){
    solution_vector_type temp; temp << rho(), u(), p(), Y();
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
  double isolate_T(solution_vector_type solution_vector) {
    return (solution_vector[2] - solution_vector[1]* solution_vector[1] /
            (2.0 * solution_vector[0])) * (gamma - 1) / solution_vector[0];
  }

  /////////////////////////////////////////////////////////////////////////
  /// \brief Function to isolate Y.
  double isolate_Y(solution_vector_type solution_vector) {
    return solution_vector[3]/solution_vector[0];
  }

  double T_var;
  double u_var;
  double rho_var;
  double Y_var;
  double gamma;
};

#endif //#endif VARIABLE_VECTOR_ISOLATOR_H
