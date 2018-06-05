#include<vector>
#include<Eigen/Core>
#include <chaiscript/chaiscript.hpp>
#include <chaiscript/chaiscript_stdlib.hpp>
#include "../Grid/Grid1D.h"
// #include "Solver/Solver.h"

using size_type = size_t;
using scalar_type = double;
using matrix_type = Eigen::Matrix<scalar_type, 4,4>;
using solution_vector_type = Eigen::Matrix<scalar_type,4,1>;
using global_solution_vector_type = std::vector<solution_vector_type>;

using Grid1D_type = Grid1D<scalar_type, size_type, global_solution_vector_type, matrix_type>;

std::string helloWorld(const std::string &t_name) {
  return "Hello " + t_name + "!";
}

int main()
{
  chaiscript::ChaiScript chai;
  chai.add(chaiscript::bootstrap::standard_library::vector_type<std::vector<solution_vector_type> >("global_solution_vector_type"));
////////////////////////////////////////////////////////////////////////////////////////////////////
//----------------------------------------------Type----------------------------------------------//
  chai.add(chaiscript::user_type<Grid1D_type>(),                          "Grid1D");
//------------------------------------------Constructors------------------------------------------//
  chai.add(chaiscript::constructor<Grid1D_type ()>(),                     "Grid1D");
  chai.add(chaiscript::constructor<Grid1D_type (const Grid1D_type &)>(),  "Grid1D");
  chai.add(chaiscript::constructor<Grid1D_type (scalar_type, scalar_type, global_solution_vector_type)>(),  "Grid1D");
  chai.add(chaiscript::constructor<Grid1D_type (scalar_type, scalar_type, scalar_type)>(),  "Grid1D");
//----------------------------------------Member Functions----------------------------------------//
  chai.add(chaiscript::fun(&Grid1D_type::number_of_cells),                  "number_of_cells");
  chai.add(chaiscript::fun(&Grid1D_type::domaine_length),                   "domaine_length");
  chai.add(chaiscript::fun(&Grid1D_type::dx),                               "dx");
  chai.add(chaiscript::fun(&Grid1D_type::per_FL),                           "per_FL");
//----------------------------------------Member Variable-----------------------------------------//
  chai.add(chaiscript::fun(&Grid1D_type::x_min),                            "x_min");
  chai.add(chaiscript::fun(&Grid1D_type::x_max),                            "x_max");
  chai.add(chaiscript::fun(&Grid1D_type::global_solution_vector),           "global_solution_vector");

////////////////////////////////////////////////////////////////////////////////////////////////////


  chai.add(chaiscript::fun(&helloWorld), "helloWorld");

try {
  chai.eval_file("test.chai");
} catch (const chaiscript::exception::eval_error &ee) {
  std::cout << ee.pretty_print();
  if ( !ee.call_stack.empty() ) {
  std::cout << "during evaluation at (" << ee.call_stack[0].start().line << ", " << ee.call_stack[0].start().column << ")";
  }
std::cout << '\n';
}
}
