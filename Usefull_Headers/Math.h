#ifndef MATH_H
#define MATH_H

namespace math{

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  template <typename T>
  T Power(const T& num, const double& expo) {return pow(num, expo);}

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  template <typename T> int sign(T val) {return (T(0) < val) - (val < T(0));}

  /////////////////////////////////////////////////////////////////////////
  /// \brief
  /// \param
  double Power(const char, const double& expo) {
    return exp(expo);
  }
} // namespace math
#endif //#ifndef MATH_H
