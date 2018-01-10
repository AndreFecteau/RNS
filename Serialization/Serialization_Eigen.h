#ifndef SERIALIZATION_EIGEN_H
#define SERIALIZATION_EIGEN_H
////////////////////////////////////////////////////////////////////////////
// This file is part of the BLawB software for the numerical solution
// of balance laws.  It is licensed under the MIT licence.  A copy of
// this license, in a file named LICENSE.md, should have been
// distributed with this file.  A copy of this license is also
// currently available at "http://opensource.org/licenses/MIT".
//
// Unless explicitly stated, all contributions intentionally submitted
// to this project shall also be under the terms and conditions of this
// license, without any additional terms or conditions.
////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Header containing serialization of Eigen matrices.
///////////////////////////////////////////////////////////////////////////
#include<stdexcept>
#include"Eigen/Dense"

namespace Eigen {

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <class Archive, class Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
void save(Archive & ar,
          const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>& m,
          const unsigned int) {
  using index_type = typename Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>::Index;
  index_type rows = m.rows();
  index_type cols = m.cols();
  ar << rows;
  ar << cols;
  for(index_type i = 0; i < rows; ++i) {
    for(index_type j = 0; j < cols; ++j) {
      ar << m(i,j);
    }
  }
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <class Archive, class Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
void load(Archive & ar,
          Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> & m,
          const unsigned int) {
  using index_type = typename Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>::Index;
  index_type rows = m.rows();
  index_type cols = m.cols();
  ar >> rows;
  ar >> cols;

  if(m.RowsAtCompileTime == Eigen::Dynamic && m.ColsAtCompileTime == Eigen::Dynamic) {
    m.resize(rows, cols);
  } else if (m.RowsAtCompileTime == Eigen::Dynamic && m.ColsAtCompileTime == cols) {
    m.resize(rows, cols);
  } else if (m.RowsAtCompileTime == rows && m.ColsAtCompileTime == Eigen::Dynamic) {
    m.resize(rows, cols);
  } else if (m.RowsAtCompileTime != rows || m.ColsAtCompileTime != cols) {
    throw std::out_of_range("Attempted to deserialize static matrix into wrong size.");
  }
  for(index_type i = 0; i < rows; ++i) {
    for(index_type j = 0; j < cols; ++j) {
      ar >> m(i,j);
    }
  }
}

#ifdef BLAWB_SERIALIZATION_BOOST
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <class Archive, class Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
void serialize(Archive & ar,
               Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> & m,
               const unsigned int version) {
  boost::serialization::split_free(ar, m, version);
}
#endif //#ifdef BLAWB_SERIALIZATION_BOOST


} //namespace Eigen

#endif //#ifndef SERIALIZATION_EIGEN_H
