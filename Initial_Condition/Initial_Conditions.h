#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H
#include "../Low_Mach_Solver/RK4_Low_Mach_Solver.h"
#include "../High_Mach_Solver/RK4_High_Mach_Solver.h"
#include "../Usefull_Headers/Helper_Functions.h"

#include <fstream>

/////////////////////////////////////////////////////////////////////////
/// \brief
/// \param
template <typename grid_type>
  void set_initial_solution(grid_type& grid, const std::function<typename grid_type::global_solution_vector_type::value_type(typename grid_type::scalar_type)> ic) {
  for(typename grid_type::size_type i = 0; i < grid.number_of_cells(); ++i) {
    grid.global_solution_vector[i] = ic((i+0.5)*grid.dx());
  }
}

/////////////////////////////////////////////////////////////////////////
/// \brief
/// \param
template <typename grid_type, typename flow_type>
void set_initial_solution(grid_type& grid, const flow_type& flow, const std::string filename) {
  typename grid_type::global_solution_vector_type primitive_variable;
  typename grid_type::scalar_type old_x, dx;
  typename grid_type::global_solution_vector_type::value_type temp;
  std::ifstream fin;
  std::string line;
  fin.open(filename.c_str());
  if (!fin) {
    throw std::runtime_error("Could not load input file " + filename);
  }
    std::getline(fin, line);
    typename grid_type::scalar_type x, rho, u, T, Y, garb;
    fin >> x >> rho >> u >> garb >> T >> Y >> garb;
    grid.x_min = x;
    temp << rho, u, T, Y;
    primitive_variable.push_back(temp);
  while(fin >> x >> rho >> u >> garb >> T >> Y >> garb){
    dx = (x-old_x);
    old_x = x;
    temp << rho, u, T, Y;
    primitive_variable.push_back(temp);
  }
  for(int i = 0; i < primitive_variable.size(); ++i){
    grid.global_solution_vector[i] << primitive_variable[i][0],
                                      primitive_variable[i][0] * primitive_variable[i][1],
                                      primitive_variable[i][0] * primitive_variable[i][2] / (flow.gamma-1) +
                                      primitive_variable[i][0] * primitive_variable[i][1] * primitive_variable[i][1] * 0.5,
                                      primitive_variable[i][0] * primitive_variable[i][3];
  }
  grid.x_max = x+dx;
}

template <typename grid_type>
typename grid_type::global_solution_vector_type::value_type
make_RK4_solution_vector(RK4_Low_Mach_Solver low_mach_solution,
                         typename grid_type::scalar_type x,
                         typename grid_type::scalar_type gamma,
                         typename grid_type::scalar_type mf) {
  typename grid_type::global_solution_vector_type::value_type temp_vec;
  temp_vec <<   low_mach_solution.get_rho(x),
                low_mach_solution.get_rho(x) * low_mach_solution.get_U(x),
                low_mach_solution.get_rho(x) * low_mach_solution.get_T(x) /(gamma*mf*mf) /
                (gamma - 1.0) +
                low_mach_solution.get_rho(x) * low_mach_solution.get_U(x) *
                low_mach_solution.get_U(x) * 0.5,
                low_mach_solution.get_rho(x) * low_mach_solution.get_Y(x);
  return temp_vec;
}

template <typename grid_type>
typename grid_type::global_solution_vector_type::value_type
make_RK4_solution_vector(RK4_High_Mach_Solver<typename grid_type::scalar_type> low_mach_solution,
                         typename grid_type::scalar_type x,
                         typename grid_type::scalar_type gamma,
                         typename grid_type::scalar_type mf) {
  typename grid_type::global_solution_vector_type::value_type temp_vec;
  temp_vec <<   low_mach_solution.get_rho(x),
                low_mach_solution.get_rho(x) * low_mach_solution.get_U(x),
                low_mach_solution.get_rho(x) * low_mach_solution.get_T(x) /(gamma*mf*mf) /
                (gamma - 1.0) +
                low_mach_solution.get_rho(x) * low_mach_solution.get_U(x) *
                low_mach_solution.get_U(x) * 0.5,
                low_mach_solution.get_rho(x) * low_mach_solution.get_Y(x);
  return temp_vec;
}

// template <typename grid_type, typename flow_type>
// void deflagration_CJ_point(grid_type& grid, flow_type& flow) {
//   using scalar_type = typename grid_type::scalar_type;
//   using size_type = typename grid_type::size_type;
//   RK4_Low_Mach_Solver initial_low_mach = RK4_Low_Mach_Solver(flow.Le, flow.Q_low_mach, flow.theta_low_mach, flow.T_ignition_scalar);
//   flow.lambda = initial_low_mach.get_lambda();
//   scalar_type safety_factor = 0.999999999999;
//   flow.mf = safety_factor*sqrt(1 + flow.Q_low_mach + flow.Q_low_mach*flow.gamma - sqrt(flow.Q_low_mach*(1 + flow.gamma)*(2 + flow.Q_low_mach + flow.Q_low_mach*flow.gamma)));
//   scalar_type p_0 = 1/(flow.gamma*flow.mf*flow.mf);
//   scalar_type zeta = (1.0 + flow.mf*flow.mf*flow.mf*flow.mf - 2.0 * flow.mf*flow.mf*(1.0+flow.Q_low_mach+flow.gamma*flow.Q_low_mach)) /
//                 (pow(pow(flow.mf, 2)*(1.0+flow.gamma),2));
//   scalar_type rho_inf = 1.0/((1.0 + pow(flow.mf,2)*flow.gamma)/(pow(flow.mf,2)*(1.0 + flow.gamma)) -
//                    sqrt(zeta));
//   scalar_type p_inf = p_0*((1.0 + pow(flow.mf,2)*flow.gamma*(1.0-1.0/rho_inf)));
//   scalar_type u_inf = 1.0/rho_inf;
//
//   std::cout << "zeta:" << zeta << "p:" << p_inf << " rho: " << rho_inf << " u: " << u_inf << std::endl;
//   scalar_type flame_location = 0.5;
//   for (size_type i = 0; i < grid.number_of_cells(); ++i) {
//     if (i*grid.dx() < grid.domaine_length()*flame_location) {
//       grid.global_solution_vector[i] << 1.0,
//       1.0 * 1.0,
//       (1.0 / (flow.gamma * flow.mf * flow.mf)) /
//       (flow.gamma - 1.0) + 1.0 * 1.0 * 1.0 * 0.5,
//       1.0 * 1.0;
//     } else if ((i+1)*grid.dx() > grid.domaine_length()*flame_location + initial_low_mach.length()) {
//       grid.global_solution_vector[i] << rho_inf,
//       rho_inf * u_inf,
//       p_inf / (flow.gamma - 1.0) + rho_inf * u_inf * u_inf * 0.5,
//       rho_inf * 0;
//     } else {
//         grid.global_solution_vector[i] << make_RK4_solution_vector<grid_type>(initial_low_mach,
//                                                         (i+1)*grid.dx() - grid.domaine_length()*flame_location, flow.gamma, flow.mf);
//     }
//   }
// }
template <typename grid_type, typename flow_type>
void deflagration_CJ_point(grid_type& grid, flow_type& flow) {
  using scalar_type = typename grid_type::scalar_type;
  using size_type = typename grid_type::size_type;
  scalar_type safety_factor = 0.999999999999;
  flow.mf = safety_factor*sqrt(1 + flow.Q_low_mach + flow.Q_low_mach*flow.gamma - sqrt(flow.Q_low_mach*(1 + flow.gamma)*(2 + flow.Q_low_mach + flow.Q_low_mach*flow.gamma)));
  RK4_High_Mach_Solver<scalar_type> initial_low_mach = RK4_High_Mach_Solver<scalar_type>(flow.Le, flow.Q_low_mach,
                                                              flow.theta_low_mach, flow.gamma,
                                                              flow.mf, flow.Pr);
  flow.lambda = initial_low_mach.get_lambda();
  // scalar_type safety_factor = 0.999999;
  scalar_type p_0 = 1/(flow.gamma*flow.mf*flow.mf);
  scalar_type zeta = (1.0 + flow.mf*flow.mf*flow.mf*flow.mf - 2.0 * flow.mf*flow.mf*(1.0+flow.Q_low_mach+flow.gamma*flow.Q_low_mach)) /
                (pow(pow(flow.mf, 2)*(1.0+flow.gamma),2));
  scalar_type rho_inf = 1.0/((1.0 + pow(flow.mf,2)*flow.gamma)/(pow(flow.mf,2)*(1.0 + flow.gamma)) -
                   sqrt(zeta));
  scalar_type p_inf = p_0*((1.0 + pow(flow.mf,2)*flow.gamma*(1.0-1.0/rho_inf)));
  scalar_type u_inf = 1.0/rho_inf;

  std::cout << "zeta:" << zeta << "p:" << p_inf << " rho: " << rho_inf << " u: " << u_inf << std::endl;
  scalar_type flame_location = 0.5;
  for (size_type i = 0; i < grid.number_of_cells(); ++i) {
    if (i*grid.dx() < grid.domaine_length()*flame_location) {
      grid.global_solution_vector[i] << 1.0,
      1.0 * 1.0,
      (1.0 / (flow.gamma * flow.mf * flow.mf)) /
      (flow.gamma - 1.0) + 1.0 * 1.0 * 1.0 * 0.5,
      1.0 * 1.0;
    } else if ((i+1)*grid.dx() > grid.domaine_length()*flame_location + initial_low_mach.length()) {
      grid.global_solution_vector[i] << rho_inf,
      rho_inf * u_inf,
      p_inf / (flow.gamma - 1.0) + rho_inf * u_inf * u_inf * 0.5,
      rho_inf * 0;
    } else {
      // std::cout << i*grid.dx() << std::endl;
        grid.global_solution_vector[i] << make_RK4_solution_vector<grid_type>(initial_low_mach,
                                                        (i+1)*grid.dx() - grid.domaine_length()*flame_location, flow.gamma, flow.mf);
    }
  }
}
template <typename grid_type, typename flow_type>
void detonation_CJ_point(grid_type& grid, flow_type& flow) {
  using scalar_type = typename grid_type::scalar_type;
  using size_type = typename grid_type::size_type;

  RK4_Low_Mach_Solver initial_low_mach = RK4_Low_Mach_Solver(flow.Le, flow.Q_low_mach, flow.theta_low_mach, flow.T_ignition_scalar);
  flow.lambda = initial_low_mach.get_lambda();
  scalar_type safety_factor = 0.999999999999;

  flow.mf = safety_factor*sqrt(1 + flow.Q_low_mach + flow.Q_low_mach*flow.gamma + sqrt(flow.Q_low_mach*(1 + flow.gamma)*(2 + flow.Q_low_mach + flow.Q_low_mach*flow.gamma)));
  scalar_type p_0 = 1/(flow.gamma*flow.mf*flow.mf);
  scalar_type zeta = 0.0;//(1.0 + flow.mf*flow.mf*flow.mf*flow.mf - 2.0 * flow.mf*flow.mf*(1.0+flow.Q_low_mach+flow.gamma*flow.Q_low_mach)) /
                // (pow(pow(flow.mf, 2)*(1.0+flow.gamma),2));
  scalar_type rho_inf = 1.0/((1.0 + pow(flow.mf,2)*flow.gamma)/(pow(flow.mf,2)*(1.0 + flow.gamma)) +
                   sqrt(zeta));
  scalar_type p_inf = p_0*((1.0 + pow(flow.mf,2)*flow.gamma*(1.0-1.0/rho_inf)));
  scalar_type u_inf = 1.0/rho_inf;

  std::cout << "zeta:" << zeta << "p:" << p_inf << " rho: " << rho_inf << " u: " << u_inf << std::endl;
  scalar_type flame_location = 0.5;
  for (size_type i = 0; i < grid.number_of_cells(); ++i) {
    if (i*grid.dx() < grid.domaine_length()*flame_location) {
      grid.global_solution_vector[i] << 1.0,
      1.0 * 1.0,
      (1.0 / (flow.gamma * flow.mf * flow.mf)) /
      (flow.gamma - 1.0) + 1.0 * 1.0 * 1.0 * 0.5,
      1.0 * 1.0;
    } else {
      grid.global_solution_vector[i] << rho_inf,
      rho_inf * u_inf,
      p_inf / (flow.gamma - 1.0) + rho_inf * u_inf * u_inf * 0.5,
      rho_inf * 0;
    }
  }
}

// template <typename grid_type, typename flow_type>
// void RK4_mf_point(grid_type& grid, flow_type& flow) {
//   using scalar_type = typename grid_type::scalar_type;
//   using size_type = typename grid_type::size_type;
//   RK4_Low_Mach_Solver initial_low_mach = RK4_Low_Mach_Solver(flow.Le, flow.Q_low_mach, flow.theta_low_mach, flow.T_ignition_scalar);
//   flow.lambda = initial_low_mach.get_lambda();
//   scalar_type safety_factor = 0.999999;
//   scalar_type p_0 = 1/(flow.gamma*flow.mf*flow.mf);
//   scalar_type zeta = (1.0 + flow.mf*flow.mf*flow.mf*flow.mf - 2.0 * flow.mf*flow.mf*(1.0+flow.Q_low_mach+flow.gamma*flow.Q_low_mach)) /
//                 (pow(pow(flow.mf, 2)*(1.0+flow.gamma),2));
//   scalar_type rho_inf = 1.0/((1.0 + pow(flow.mf,2)*flow.gamma)/(pow(flow.mf,2)*(1.0 + flow.gamma)) -
//                    sqrt(zeta));
//   scalar_type p_inf = p_0*((1.0 + pow(flow.mf,2)*flow.gamma*(1.0-1.0/rho_inf)));
//   scalar_type u_inf = 1.0/rho_inf;
//
//   std::cout << "zeta:" << zeta << "p:" << p_inf << " rho: " << rho_inf << " u: " << u_inf << std::endl;
//   scalar_type flame_location = 0.5;
//   for (size_type i = 0; i < grid.number_of_cells(); ++i) {
//     if (i*grid.dx() < grid.domaine_length()*flame_location) {
//       grid.global_solution_vector[i] << 1.0,
//       1.0 * 1.0,
//       (1.0 / (flow.gamma * flow.mf * flow.mf)) /
//       (flow.gamma - 1.0) + 1.0 * 1.0 * 1.0 * 0.5,
//       1.0 * 1.0;
//     } else if ((i+1)*grid.dx() > grid.domaine_length()*flame_location + initial_low_mach.length()) {
//       grid.global_solution_vector[i] << rho_inf,
//       rho_inf * u_inf,
//       p_inf / (flow.gamma - 1.0) + rho_inf * u_inf * u_inf * 0.5,
//       rho_inf * 0;
//     } else {
//         grid.global_solution_vector[i] << make_RK4_solution_vector<grid_type>(initial_low_mach,
//                                                         (i+1)*grid.dx() - grid.domaine_length()*flame_location, flow.gamma, flow.mf);
//     }
//   }
// }

template <typename grid_type, typename flow_type>
void RK4_mf_point(grid_type& grid, flow_type& flow) {
  using scalar_type = typename grid_type::scalar_type;
  using size_type = typename grid_type::size_type;
  RK4_High_Mach_Solver<scalar_type> initial_low_mach = RK4_High_Mach_Solver<scalar_type>(flow.Le, flow.Q_low_mach,
                                                              flow.theta_low_mach, flow.gamma,
                                                              flow.mf, flow.Pr);
  flow.lambda = initial_low_mach.get_lambda();
  scalar_type safety_factor = 0.999999;
  scalar_type p_0 = 1/(flow.gamma*flow.mf*flow.mf);
  scalar_type zeta = (1.0 + flow.mf*flow.mf*flow.mf*flow.mf - 2.0 * flow.mf*flow.mf*(1.0+flow.Q_low_mach+flow.gamma*flow.Q_low_mach)) /
                (pow(pow(flow.mf, 2)*(1.0+flow.gamma),2));
  scalar_type rho_inf = 1.0/((1.0 + pow(flow.mf,2)*flow.gamma)/(pow(flow.mf,2)*(1.0 + flow.gamma)) -
                   sqrt(zeta));
  scalar_type p_inf = p_0*((1.0 + pow(flow.mf,2)*flow.gamma*(1.0-1.0/rho_inf)));
  scalar_type u_inf = 1.0/rho_inf;

  std::cout << "zeta:" << zeta << "p:" << p_inf << " rho: " << rho_inf << " u: " << u_inf << std::endl;
  scalar_type flame_location = 0.5;
  for (size_type i = 0; i < grid.number_of_cells(); ++i) {
    if (i*grid.dx() < grid.domaine_length()*flame_location) {
      grid.global_solution_vector[i] << 1.0,
      1.0 * 1.0,
      (1.0 / (flow.gamma * flow.mf * flow.mf)) /
      (flow.gamma - 1.0) + 1.0 * 1.0 * 1.0 * 0.5,
      1.0 * 1.0;
    } else if ((i+1)*grid.dx() > grid.domaine_length()*flame_location + initial_low_mach.length()) {
      grid.global_solution_vector[i] << rho_inf,
      rho_inf * u_inf,
      p_inf / (flow.gamma - 1.0) + rho_inf * u_inf * u_inf * 0.5,
      rho_inf * 0;
    } else {
      // std::cout << i*grid.dx() << std::endl;
        grid.global_solution_vector[i] << make_RK4_solution_vector<grid_type>(initial_low_mach,
                                                        (i+1)*grid.dx() - grid.domaine_length()*flame_location, flow.gamma, flow.mf);
    }
  }
}

template <typename grid_type, typename flow_type>
void manufactured_solution(grid_type& grid, flow_type& flow){
grid.x_min = 0.0;
grid.x_max = 2.0*atan(1.0)*4.0;
// grid.number_of_cells() = ;
// std::cout << dx << std::endl;
grid.global_solution_vector.resize(grid.number_of_cells());
// scalar_type dx = (x_max-x_min) / number_of_cells;
  for (size_t i = 0; i < grid.number_of_cells(); ++i) {
    typename grid_type::scalar_type x = (i+0.5)*grid.dx();
    grid.global_solution_vector[i] << (-0.45*tanh(4.0 * x - 10.0) + 0.55),
                          (-0.45*tanh(4.0 * x - 10.0) + 0.55)*(4.5*tanh(4.0 * x - 10.0) + 5.5),
                          2.0*tanh(4.0*x - 10.0) + 70000,
                          (-0.45*tanh(4.0 * x - 10.0) + 0.55)*(-0.5*tanh(x - 8.0/4.0) + 0.5);
  }
}
#endif //#ifndef INITIAL_CONDITIONS_H
