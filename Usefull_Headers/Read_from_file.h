#ifndef READ_FROM_FILE_H
#define READ_FROM_FILE_H

// #include <omp.h>
// #include <math.h>
// #include <limits>
#include <vector>
#include "Eigen/Core"
#include "Eigen/Dense"

class Read_from_file {

  /////////////////////////////////////////////////////////////////////////
  /// \brief type for individual cell solution vector.
  typedef Eigen::Matrix<double, 4, 1> solution_vector_type;
  using global_solution_vector_type = std::vector<solution_vector_type>;

 public:
  /////////////////////////////////////////////////////////////////////////
  /// \brief Default constructor.
  Read_from_file() = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy constructor.
  Read_from_file(const Read_from_file&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move constructor.
  Read_from_file(Read_from_file&&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Copy assignment operator.
  Read_from_file& operator=(const Read_from_file&) = default;

  /////////////////////////////////////////////////////////////////////////
  /// \brief Move assignment operator.
  Read_from_file& operator=(Read_from_file&&) = default;

  Read_from_file(std::string filename) {
        solution_vector_type temp;
        std::ifstream fin;
        std::string line;
        fin.open(filename.c_str());
        if (!fin) {
          throw std::runtime_error("Could not load input file " + filename);
        }
          std::getline(fin, line);
          double x, rho, u, T, Y, garb;
          fin >> x >> rho >> u >> garb >> T >> Y >> garb;
          x_vector.push_back(x);
          temp << rho, u, T, Y;
          primitive_variable.push_back(temp);
        while(fin >> x >> rho >> u >> garb >> T >> Y >> garb){
          temp << rho, u, T, Y;
          primitive_variable.push_back(temp);
        }
  }

  Eigen::Matrix<double, 4, 1> interpolate(double location);

 private:
  global_solution_vector_type primitive_variable;
  std::vector<double> x_vector;

  int get_i(double location);

};

///////////////////////////////////////////////////////////////////////////////
//Gets a value between nodes.
///////////////////////////////////////////////////////////////////////////////
Eigen::Matrix<double, 4, 1> Read_from_file::interpolate(double location) {
  int i = get_i(location);
  return (primitive_variable[i+1] - primitive_variable[i]) * (location - x_vector[i]) / (x_vector[i+1] - x_vector[i]) + primitive_variable[i];
}

int Read_from_file::get_i(double location) {
  int i = 0;
  while(((x_vector[i] - location) * (x_vector[i+1] - location)) > 0){
    ++i;
  }
  return i;
}
#endif //#ifndef READ_FROM_FILE_H
