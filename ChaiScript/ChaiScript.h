#define HYPERBOLIC
#define VISCOUS
#define SOURCE
#define MANUFACTURED
// #define RECENTER_FLAME

#include<vector>
#include <Eigen/Core>
#include <chaiscript/chaiscript.hpp>
#include <chaiscript/chaiscript_stdlib.hpp>
#include <chaiscript/extras/math.hpp>

#include "../Solver/Implicit_Marching.h"
#include "../Implicit_Flux_and_Sources/Implicit_HLLE.h"
#include "../Implicit_Flux_and_Sources/Implicit_Centered_Difference_2nd_Order.h"
#include "../Grid/Grid1D.h"
#include "../Physical_Property/Non_Dimensional_Navier_Stokes.h"
#include "../Solver/Solver.h"
#include "../Initial_Condition/Initial_Conditions.h"
#include "../Serialization/Serialize.h"
#include "../Usefull_Headers/Handle_Itterative_Lambda.h"

using size_type = size_t;
using scalar_type = long double;
using matrix_type = Eigen::Matrix<scalar_type, 4,4>;
using solution_vector_type = Eigen::Matrix<scalar_type,4,1>;
using global_solution_vector_type = std::vector<solution_vector_type>;
using grid1D_type = Grid1D<scalar_type, size_type, global_solution_vector_type, matrix_type>;
using flow_type = Non_Dimensional_Navier_Stokes<scalar_type>;
using time_stepping_type = Implicit_Marching<grid1D_type, flow_type>;
using HLLE_flux_type = Implicit_HLLE<grid1D_type, flow_type>;
using CD2_flux_type = Implicit_Centered_Difference_2nd_Order<grid1D_type, flow_type>;
using HLLE_solver_type = Solver<flow_type, grid1D_type, HLLE_flux_type, time_stepping_type>;
using CD2_solver_type = Solver<flow_type, grid1D_type, CD2_flux_type, time_stepping_type>;

void register_RNS_Chaiscript(chaiscript::ChaiScript& chai) {

  chai.add(chaiscript::bootstrap::standard_library::vector_type<std::vector<solution_vector_type> >("global_solution_vector_type"));
  chai.add(chaiscript::user_type<HLLE_flux_type>(), "HLLE");
  chai.add(chaiscript::user_type<CD2_flux_type>(), "Centered_Difference");
  chai.add(chaiscript::user_type<time_stepping_type>(), "Stepping");
  auto mathlib = chaiscript::extras::math::bootstrap();
  chai.add(mathlib);
///////////////////////////////////////////////Eigen////////////////////////////////////////////////
//----------------------------------------------Type----------------------------------------------//
  chai.add(chaiscript::user_type<solution_vector_type>(),               "Solution_vector_type");
//------------------------------------------Constructors------------------------------------------//
  chai.add(chaiscript::constructor<solution_vector_type ()>(),          "Solution_vector_type");
  chai.add(chaiscript::constructor<solution_vector_type(const solution_vector_type&)>(),   "Solution_vector_type");
//----------------------------------------Member Functions----------------------------------------//
  chai.add(chaiscript::fun([](const solution_vector_type& v){return v.size();}),"size");
  auto subscript0 = [](solution_vector_type& v, typename solution_vector_type::Index i) -> typename solution_vector_type::Scalar& {
    if(i < 0) {
      throw std::length_error("Negative index given to vector subscript.");
    } else if(i >= v.size()) {
      throw std::length_error("Index given to vector subscript is out of range.");
    }
    return v[i];
  };
  const auto assignment = [](solution_vector_type& a, const solution_vector_type& b) {
    return a=b;
  };
  const auto addition = [](const solution_vector_type& a, const solution_vector_type& b) {
    if(a.size()!=b.size()) {
      throw std::length_error("Error, tried to add vectors of different sizes");
    }
    return (a+b).eval();
  };
  const auto subtraction = [](const solution_vector_type& a, const solution_vector_type& b) {
    if(a.size()!=b.size()) {
      throw std::length_error("Error, tried to subtract vectors of different sizes");
    }
    return (a-b).eval();
  };
  const auto multiplication1 = [](const solution_vector_type& a, const scalar_type& b) {
    return (a*b).eval();
  };
  const auto multiplication2 = [](const scalar_type& a, const solution_vector_type& b) {
    return (a*b).eval();
  };
  const auto division = [](const solution_vector_type& a, const scalar_type& b) {
    return (a/b).eval();
  };
  const auto vec_to_string = [](const solution_vector_type& v) {
    if(v.size()==0) return std::string("[]");
    std::stringstream oss;
    oss.precision(16);
    oss << "[" << v[0];
    for(int i = 1; i< v.size(); ++i) {
      oss << ", " << v[i];
    }
    oss << "]";
    return oss.str();
  };
  chai.add(chaiscript::fun(subscript0),"[]");
  chai.add(chaiscript::fun(assignment),"=");
  chai.add(chaiscript::fun(addition),"+");
  chai.add(chaiscript::fun(subtraction),"-");
  chai.add(chaiscript::fun(multiplication1),"*");
  chai.add(chaiscript::fun(multiplication2),"*");
  chai.add(chaiscript::fun(division),"/");
  chai.add(chaiscript::fun([](solution_vector_type& v, const scalar_type d) {v.fill(d);}),"fill");
  chai.add(chaiscript::fun(vec_to_string),"to_string");

/////////////////////////////////////////Initial_Condition//////////////////////////////////////////
//---------------------------------------General Functions----------------------------------------//
chai.add(chaiscript::fun(&set_initial_solution<grid1D_type>),             "set_initial_solution");
chai.add(chaiscript::fun(&set_initial_solution<grid1D_type, flow_type>),  "set_initial_solution");
chai.add(chaiscript::fun(&deflagration_CJ_point<grid1D_type, flow_type>), "deflagration_CJ_point");
chai.add(chaiscript::fun(&detonation_CJ_point<grid1D_type, flow_type>),   "detonation_CJ_point");
chai.add(chaiscript::fun(&RK4_mf_point<grid1D_type, flow_type>),          "deflagration_mf_point");
chai.add(chaiscript::fun(&Incompressible_mf_point<grid1D_type, flow_type>),"Incompressible_mf_point");
chai.add(chaiscript::fun(&manufactured_solution<grid1D_type, flow_type>), "manufactured_solution");
chai.add(chaiscript::fun(&unserialize_to_file<HLLE_solver_type>),         "unserialize_to_file");
chai.add(chaiscript::fun(&unserialize_to_file<CD2_solver_type>),          "unserialize_to_file");

/////////////////////////////////////Handle_Itterative_Lambda///////////////////////////////////////
//---------------------------------------General Functions----------------------------------------//
chai.add(chaiscript::fun(&bisection_lambda),             "bisection_lambda");
chai.add(chaiscript::fun(&add_lambda_gap),             "add_lambda_gap");

///////////////////////////////////////////////Grid1D///////////////////////////////////////////////
//----------------------------------------------Type----------------------------------------------//
  chai.add(chaiscript::user_type<grid1D_type>(),                          "Grid1D");
//------------------------------------------Constructors------------------------------------------//
  chai.add(chaiscript::constructor<grid1D_type ()>(),                     "Grid1D");
  chai.add(chaiscript::constructor<grid1D_type (const grid1D_type &)>(),  "Grid1D");
  chai.add(chaiscript::constructor<grid1D_type (scalar_type, scalar_type, global_solution_vector_type)>(),  "Grid1D");
  chai.add(chaiscript::constructor<grid1D_type (scalar_type, scalar_type, scalar_type)>(),  "Grid1D");
//----------------------------------------Member Functions----------------------------------------//
  chai.add(chaiscript::fun(&grid1D_type::number_of_cells),                  "number_of_cells");
  chai.add(chaiscript::fun(&grid1D_type::domaine_length),                   "domaine_length");
  chai.add(chaiscript::fun(&grid1D_type::dx),                               "dx");
  chai.add(chaiscript::fun(&grid1D_type::per_FL),                           "per_FL");
  chai.add(chaiscript::fun(&grid1D_type::add_space_in_back),                "add_space_in_back");
  chai.add(chaiscript::fun(&grid1D_type::add_space_in_front),               "add_space_in_front");
  chai.add(chaiscript::fun(&grid1D_type::refine),                           "refine");

//----------------------------------------Member Variable-----------------------------------------//
  chai.add(chaiscript::fun(&grid1D_type::x_min),                            "x_min");
  chai.add(chaiscript::fun(&grid1D_type::x_max),                            "x_max");
  chai.add(chaiscript::fun(&grid1D_type::global_solution_vector),           "global_solution_vector");

//////////////////////////////////////////////flow_type/////////////////////////////////////////////
//----------------------------------------------Type----------------------------------------------//
  chai.add(chaiscript::user_type<flow_type>(),                          "Flow");
//------------------------------------------Constructors------------------------------------------//
  chai.add(chaiscript::constructor<flow_type ()>(),                     "Flow");
  chai.add(chaiscript::constructor<flow_type (const flow_type &)>(),    "Flow");
  chai.add(chaiscript::constructor<flow_type (scalar_type, scalar_type, scalar_type,
                                              scalar_type, scalar_type, scalar_type,
                                              scalar_type, scalar_type)>(),  "Flow");
//----------------------------------------Member Functions----------------------------------------//
  chai.add(chaiscript::fun(&flow_type::T_ignition),                   "T_ignition");
  chai.add(chaiscript::fun(&flow_type::Q),                            "Q");
  chai.add(chaiscript::fun(&flow_type::theta),                        "theta");
//----------------------------------------Member Variable-----------------------------------------//
  chai.add(chaiscript::fun(&flow_type::Pr),                            "Pr");
  chai.add(chaiscript::fun(&flow_type::Le),                            "Le");
  chai.add(chaiscript::fun(&flow_type::Q_low_mach),                    "Q_low_mach");
  chai.add(chaiscript::fun(&flow_type::theta_low_mach),                "theta_low_mach");
  chai.add(chaiscript::fun(&flow_type::mf),                            "mf");
  chai.add(chaiscript::fun(&flow_type::gamma),                         "gamma");
  chai.add(chaiscript::fun(&flow_type::lambda),                        "lambda");
  chai.add(chaiscript::fun(&flow_type::T_ignition_scalar),             "T_ignition_scalar");

/////////////////////////////////////////////Solver_HLLE////////////////////////////////////////////
//----------------------------------------------Type----------------------------------------------//
  chai.add(chaiscript::user_type<HLLE_solver_type>(),                          "HLLE_solver");
//------------------------------------------Constructors------------------------------------------//
  chai.add(chaiscript::constructor<HLLE_solver_type ()>(),                     "HLLE_solver");
  chai.add(chaiscript::constructor<HLLE_solver_type (const HLLE_solver_type &)>(),  "HLLE_solver");
  chai.add(chaiscript::constructor<HLLE_solver_type (flow_type, grid1D_type, scalar_type,
                                              scalar_type, scalar_type, scalar_type,
                                              std::string, scalar_type, scalar_type)>(),  "HLLE_solver");
//----------------------------------------Member Functions----------------------------------------//
  chai.add(chaiscript::fun(&HLLE_solver_type::solve),                        "solve");
  chai.add(chaiscript::fun(&HLLE_solver_type::get_flow),                     "get_flow");
  chai.add(chaiscript::fun(&HLLE_solver_type::change_flow),                  "change_flow");
  chai.add(chaiscript::fun(&HLLE_solver_type::get_grid),                     "get_grid");
  chai.add(chaiscript::fun(&HLLE_solver_type::get_lambda),                   "get_lambda");
  chai.add(chaiscript::fun(&HLLE_solver_type::change_lambda),                "change_lambda");
  chai.add(chaiscript::fun(&HLLE_solver_type::recenter_solution),            "recenter_solution");
  // chai.add(chaiscript::fun(&HLLE_solver_type::add_space_in_back),            "add_space_in_back");
  // chai.add(chaiscript::fun(&HLLE_solver_type::add_space_in_front),           "add_space_in_front");
  chai.add(chaiscript::fun(&HLLE_solver_type::print_stats),                  "print_stats");
  chai.add(chaiscript::fun(&HLLE_solver_type::change_filename),              "change_filename");
  chai.add(chaiscript::fun(&HLLE_solver_type::reset_frame_number),           "reset_frame_number");
  // chai.add(chaiscript::fun(&HLLE_solver_type::refine),                       "refine");
  chai.add(chaiscript::fun(&HLLE_solver_type::change_frame_time),            "change_frame_time");
  chai.add(chaiscript::fun(&HLLE_solver_type::change_CFL),                   "change_CFL");
  chai.add(chaiscript::fun(&HLLE_solver_type::change_flame_location),        "change_flame_location");
  chai.add(chaiscript::fun(&HLLE_solver_type::plot_limiter),                 "plot_limiter");
  chai.add(chaiscript::fun(&HLLE_solver_type::grid),                         "grid");
//----------------------------------------Member Variable-----------------------------------------//
  chai.add(chaiscript::fun(&HLLE_solver_type::plot_global_solution_vector),  "plot_global_solution_vector");

/////////////////////////////////////////////Solver_HLLE////////////////////////////////////////////
//----------------------------------------------Type----------------------------------------------//
  chai.add(chaiscript::user_type<CD2_solver_type>(),                          "CD2_solver");
//------------------------------------------Constructors------------------------------------------//
  chai.add(chaiscript::constructor<CD2_solver_type ()>(),                     "CD2_solver");
  chai.add(chaiscript::constructor<CD2_solver_type (const CD2_solver_type &)>(),  "CD2_solver");
  chai.add(chaiscript::constructor<CD2_solver_type (flow_type, grid1D_type, scalar_type,
                                              scalar_type, scalar_type, scalar_type,
                                              std::string, scalar_type, scalar_type)>(),  "CD2_solver");
//----------------------------------------Member Functions----------------------------------------//
  chai.add(chaiscript::fun(&CD2_solver_type::solve),                        "solve");
  chai.add(chaiscript::fun(&CD2_solver_type::get_flow),                     "get_flow");
  chai.add(chaiscript::fun(&CD2_solver_type::change_flow),                  "change_flow");
  chai.add(chaiscript::fun(&CD2_solver_type::get_grid),                     "get_grid");
  chai.add(chaiscript::fun(&CD2_solver_type::get_lambda),                   "get_lambda");
  chai.add(chaiscript::fun(&CD2_solver_type::recenter_solution),            "recenter_solution");
  // chai.add(chaiscript::fun(&CD2_solver_type::add_space_in_back),            "add_space_in_back");
  // chai.add(chaiscript::fun(&CD2_solver_type::add_space_in_front),           "add_space_in_front");
  chai.add(chaiscript::fun(&CD2_solver_type::print_stats),                  "print_stats");
  chai.add(chaiscript::fun(&CD2_solver_type::change_filename),              "change_filename");
  chai.add(chaiscript::fun(&CD2_solver_type::reset_frame_number),           "reset_frame_number");
  chai.add(chaiscript::fun(&CD2_solver_type::change_lambda),                "change_lambda");
  // chai.add(chaiscript::fun(&CD2_solver_type::refine),                       "refine");
  chai.add(chaiscript::fun(&CD2_solver_type::change_frame_time),            "change_frame_time");
  chai.add(chaiscript::fun(&CD2_solver_type::change_CFL),                   "change_CFL");
  chai.add(chaiscript::fun(&CD2_solver_type::change_flame_location),        "change_flame_location");
  chai.add(chaiscript::fun(&CD2_solver_type::plot_limiter),                 "plot_limiter");
  chai.add(chaiscript::fun(&CD2_solver_type::plot_global_solution_vector),  "plot_global_solution_vector");
//----------------------------------------Member Variable-----------------------------------------//
  chai.add(chaiscript::fun(&CD2_solver_type::grid),                         "grid");

////////////////////////////////////////////////////////////////////////////////////////////////////
}

inline void execute_chaiscript_file(std::string filename) {

  chaiscript::ChaiScript chai;
  // chaiscript::ChaiScript chai(chaiscript::Std_Lib::library());

  register_RNS_Chaiscript(chai);

  try {
    chai.eval_file(filename);
  } catch (const chaiscript::exception::eval_error &ee) {
   std::cout << ee.pretty_print();
   if ( !ee.call_stack.empty() ) {
   std::cout << "during evaluation at (" << ee.call_stack[0].start().line << ", " << ee.call_stack[0].start().column << ")";
  }
  std::cout << '\n';
  }
}
