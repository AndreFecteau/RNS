#ifndef BLOCK_TRIAGONAL_MATRIX_INVERSE_H
#define BLOCK_TRIAGONAL_MATRIX_INVERSE_H
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief  Header containg the definition of the class for a
///         block_triagonal_matrix_inverse.h
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
  /// \brief  Header containg the Triag_Inverse solver for first
  ///         order terms.
///////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <iostream>
#include "Eigen/Core"
#include "Eigen/Dense"

template <typename matrix_type, typename vector_type>
std::vector<vector_type> block_triagonal_matrix_inverse(std::vector<matrix_type>& diag,
                                                        std::vector<matrix_type>& diag_top,
                                                        std::vector<matrix_type>& diag_bot,
                                                        std::vector<vector_type> global_solution_vector) {

  std::vector<vector_type> delta_global_solution_vector(diag.size());
  diag_top[0] = diag[0].inverse() * diag_top[0];
  global_solution_vector[0] = diag[0].inverse() * global_solution_vector[0];

  for(size_t i = 1; i < diag.size(); ++i){
    diag_top[i] = (diag[i]-diag_bot[i]*diag_top[i-1]).inverse() * diag_top[i];
    global_solution_vector[i] = (diag[i] - diag_bot[i]*diag_top[i-1]).inverse()*
                                (global_solution_vector[i] - diag_bot[i]*global_solution_vector[i-1]);
  }

  delta_global_solution_vector[diag.size()-1] = global_solution_vector[diag.size()-1];
  for(int i = diag.size() - 2; i >= 0; --i) {
    delta_global_solution_vector[i] = global_solution_vector[i] -
                                      diag_top[i] * delta_global_solution_vector[i+1];
  }

return delta_global_solution_vector;
}

#endif //#ifndef BLOCK_TRIAGONAL_MATRIX_INVERSE_H
