// #include "Block_Triagonal_Matrix_Inverse.h"
#include "Matrix_Inverse/Gaussian_Block_Penagonal_Matrix_Inverse.h"
#include "Eigen/Core"
#include "Eigen/Dense"

int main() {

std::cout << "//////////////////////////////////////////////" << std::endl;
std::cout << "1X1 matrices" << std::endl;
std::cout << "//////////////////////////////////////////////" << std::endl;
  {
  using matrix_type  = Eigen::Matrix<double, 1, 1>;
  using vector_type  = Eigen::Matrix<double, 1, 1>;
  matrix_type a; a << 1;
  matrix_type b; b << -2;
  matrix_type c; c << 1;
  matrix_type d; d << 0.5;
  vector_type sol; sol << 1;
  std::vector<vector_type>  global_solution_vector = {sol, sol, sol, sol};
  std::vector<matrix_type> diag_2top = {d, d, d, d};
  std::vector<matrix_type> diag_top = {c, c, c, c};
  std::vector<matrix_type> diag = {b, b, b, b};
  std::vector<matrix_type> diag_bot = {a, a, a, a};
  std::vector<matrix_type> diag_2bot = {d, d, d, d};
  std::vector<vector_type> matrix = block_triagonal_matrix_inverse<matrix_type>(diag_2top, diag_top, diag, diag_bot, diag_2bot, global_solution_vector);
  std::cout << matrix[0](0,0) << "," << " to be -2" << std::endl;
  std::cout << matrix[1](0,0) << "," << " to be -3" << std::endl;
  std::cout << matrix[2](0,0) << "," << " to be -3" << std::endl;
  std::cout << matrix[3](0,0) << "," << " to be -2" << std::endl;
  }
//
//
// std::cout << "//////////////////////////////////////////////" << std::endl;
// std::cout << "1X1 matrices" << std::endl;
// std::cout << "//////////////////////////////////////////////" << std::endl;
//   {
//   using matrix_type  = Eigen::Matrix<double, 1, 1>;
//   using vector_type  = Eigen::Matrix<double, 1, 1>;
//   matrix_type a; a << 1;
//   matrix_type b; b << -4;
//   matrix_type c; c << -1;
//   vector_type sol; sol << 1;
//   std::vector<vector_type>  global_solution_vector = {sol, sol, sol, sol};
//   std::vector<matrix_type> diag = {a, b, c, a};
//   std::vector<matrix_type> diag_top = {b, a, b, c};
//   std::vector<matrix_type> diag_bot = {c, b, c, a};
//   std::vector<vector_type> matrix = block_triagonal_matrix_inverse<matrix_type>(diag, diag_top, diag_bot, global_solution_vector);
//   std::cout << matrix[0](0,0) << "," << " to be 0.322034" << std::endl;
//   std::cout << matrix[1](0,0) << "," << " to be -0.169492" << std::endl;
//   std::cout << matrix[2](0,0) << "," << " to be 1.61017" << std::endl;
//   std::cout << matrix[3](0,0) << "," << " to be -0.610169" << std::endl;
//   }
//
//
std::cout << "//////////////////////////////////////////////" << std::endl;
std::cout << "2X2 matrices" << std::endl;
std::cout << "//////////////////////////////////////////////" << std::endl;
 {
  using matrix_type = Eigen::Matrix<double, 2, 2>;
  using vector_type = Eigen::Matrix<double, 2, 1>;
  matrix_type a; a << 1, 2,
                      3, 4;
  matrix_type b; b << 2, 6,
                      9, 1;
  matrix_type c; c << 5, 1,
                      0, 4;
  matrix_type d; d << 1, 2,
                      2, 1;
  vector_type sol; sol << 1, 1;
  std::vector<vector_type>  global_solution_vector = {sol, sol, sol, sol};
  std::vector<matrix_type> diag_2top = {d, d, d, d};
  std::vector<matrix_type> diag_top = {c, c, c, c};
  std::vector<matrix_type> diag = {b, b, b, b};
  std::vector<matrix_type> diag_bot = {a, a, a, a};
  std::vector<matrix_type> diag_2bot = {d, d, d, d};
  std::cout << a.eigenvalues() << std::endl;
  std::cout << b.eigenvalues() << std::endl;
  std::cout << c.eigenvalues() << std::endl;
  std::vector<vector_type> matrix = block_triagonal_matrix_inverse<matrix_type>(diag_2top, diag_top, diag, diag_bot, diag_2bot, global_solution_vector);
  std::cout << matrix[0](0,0) << ", to be 0.0503721" <<  std::endl;
  std::cout << matrix[0](1,0) << ", to be 0.127831" <<  std::endl;
  std::cout << matrix[1](0,0) << ", to be 0.0055125" <<  std::endl;
  std::cout << matrix[1](1,0) << ", to be 0.104705" <<  std::endl;
  std::cout << matrix[2](0,0) << ", to be 0.00178003" <<  std::endl;
  std::cout << matrix[2](1,0) << ", to be 0.0458102" <<  std::endl;
  std::cout << matrix[3](0,0) << ", to be 0.0761907" <<  std::endl;
  std::cout << matrix[3](1,0) << ", to be 0.125703" <<  std::endl;
  }
// std::cout << "//////////////////////////////////////////////" << std::endl;
// std::cout << "3X3 matrices" << std::endl;
// std::cout << "//////////////////////////////////////////////" << std::endl;
//  {
//   using matrix_type = Eigen::Matrix<double, 3, 3>;
//   using vector_type = Eigen::Matrix<double, 3, 1>;
//   matrix_type a; a << 2, 3, 2,
//                       2, 3, 2,
//                       1, 3, 6;
//   matrix_type b; b << 2, 3, 5,
//                       1, 6, 4,
//                       3, 9, 6;
//   matrix_type c; c << 2, 3, 1,
//                       2, 3, 2,
//                       2, 3, 6;
//   vector_type sol; sol << 1, 1, 1;
//   std::vector<vector_type>  global_solution_vector = {sol, sol, sol, sol};
//   std::vector<matrix_type> diag = {b, b, b, b};
//   std::vector<matrix_type> diag_top = {c, c, c, c};
//   std::vector<matrix_type> diag_bot = {a, a, a, a};
//   std::vector<vector_type> matrix = block_triagonal_matrix_inverse<matrix_type>(diag, diag_top, diag_bot, global_solution_vector);
//   std::cout << matrix[0](0,0) << ", to be 0.214286" <<  std::endl;
//   std::cout << matrix[0](1,0) << ", to be 0.142857" <<  std::endl;
//   std::cout << matrix[0](2,0) << ", to be 0.0" <<  std::endl;
//   std::cout << matrix[1](0,0) << ", to be 0.214286" <<  std::endl;
//   std::cout << matrix[1](1,0) << ", to be -0.0238095" <<  std::endl;
//   std::cout << matrix[1](2,0) << ", to be -0.214286" <<  std::endl;
//   std::cout << matrix[2](0,0) << ", to be 0.25" <<  std::endl;
//   std::cout << matrix[2](1,0) << ", to be 0.0952381" <<  std::endl;
//   std::cout << matrix[2](2,0) << ", to be 0.0714286" <<  std::endl;
//   std::cout << matrix[3](0,0) << ", to be -0.047619" <<  std::endl;
//   std::cout << matrix[3](1,0) << ", to be -0.00396825" <<  std::endl;
//   std::cout << matrix[3](2,0) << ", to be 0.0357143" <<  std::endl;
//   }
}
