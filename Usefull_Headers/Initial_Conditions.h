#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H
#include <fstream>


template <typename T>
std::string tostring(T name) {
  return std::to_string(static_cast<int>(name));
}
/////////////////////////////////////////////////////////////////////////
/// \brief
/// \param
template <typename grid_type>
void set_initial_solution(grid_type& grid, const std::function<solution_vector_size(double)> ic) {
  for(typename grid::size_type i = 0; i < grid.number_of_cells(); ++i){
    grid.global_solution_vector[i] = ic((i+0.5)*grid.dx());
  }
}

/////////////////////////////////////////////////////////////////////////
/// \brief
/// \param
template <typename grid_type, typename flow_type>
void set_initial_solution(grid_type& grid, const flow_type& flow, const std::string filename) {
  typename grid_type::global_solution_vector_type primitive_variable;
  double old_x, dx;
  typename grid_type::global_solution_vector_type::value_type temp;
  std::ifstream fin;
  std::string line;
  fin.open(filename.c_str());
  if (!fin) {
    throw std::runtime_error("Could not load input file " + filename);
  }
    std::getline(fin, line);
    double x, rho, u, T, Y, garb;
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
  for(int i = 0; i < primitive_variable.size()-100; ++i){
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

template <typename flow_properties_type, typename grid_type>
void load_from_file(flow_properties_type& flow, grid_type& grid, std::string filename) {
  using global_solution_vector_type = typename grid_type::global_solution_vector_type;
  using solution_vector_type = typename grid_type::global_solution_vector_type::value_type;

  global_solution_vector_type primitive_variable;
  double old_x, dx;
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
    grid.x_min = x;
    temp << rho, u, T, Y;
    primitive_variable.push_back(temp);
  while(fin >> x >> rho >> u >> garb >> T >> Y >> garb){
    dx = (x-old_x);
    old_x = x;
    temp << rho, u, T, Y;
    primitive_variable.push_back(temp);
  }
  for(int i = 0; i < primitive_variable.size()-100; ++i){
    grid.global_solution_vector[i] << primitive_variable[i][0],
                                      primitive_variable[i][0] * primitive_variable[i][1],
                                      primitive_variable[i][0] * primitive_variable[i][2] / (flow.gamma-1) +
                                      primitive_variable[i][0] * primitive_variable[i][1] * primitive_variable[i][1] * 0.5,
                                      primitive_variable[i][0] * primitive_variable[i][3];
  }
  grid.x_max = x+dx;
}


template <typename flow_properties_type, typename grid_type>
void deflagration_CJ_point(flow_properties_type& flow, grid_type& grid) {
  using scalar_type = typename grid_type::scalar_type;
  using size_type = typename grid_type::size_type;
  RK4_Low_Mach_Solver initial_low_mach = RK4_Low_Mach_Solver(flow.Le, flow.Q_low_mach, flow.theta_low_mach, flow.T_ignition_scalar);
  flow.lambda = initial_low_mach.get_lambda();
  scalar_type safety_factor = 0.999999999999;
  flow.mf = safety_factor*sqrt(1 + flow.Q_low_mach + flow.Q_low_mach*flow.gamma - sqrt(flow.Q_low_mach*(1 + flow.gamma)*(2 + flow.Q_low_mach + flow.Q_low_mach*flow.gamma)));
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
    } else {
      grid.global_solution_vector[i] << rho_inf,
      rho_inf * u_inf,
      p_inf / (flow.gamma - 1.0) + rho_inf * u_inf * u_inf * 0.5,
      rho_inf * 0;
    }
  }
}

template <typename flow_properties_type, typename grid_type>
void detonation_CJ_point(flow_properties_type& flow, grid_type& grid) {
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

template <typename flow_properties_type, typename grid_type>
void RK4_mf_point(flow_properties_type& flow, grid_type& grid) {
  using scalar_type = typename grid_type::scalar_type;
  using size_type = typename grid_type::size_type;
  RK4_Low_Mach_Solver initial_low_mach = RK4_Low_Mach_Solver(flow.Le, flow.Q_low_mach, flow.theta_low_mach, flow.T_ignition_scalar);
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
        grid.global_solution_vector[i] << make_RK4_solution_vector<grid_type>(initial_low_mach,
                                                        (i+1)*grid.dx() - grid.domaine_length()*flame_location, flow.gamma, flow.mf);
    }
  }
}
//
// void straight_line(int number_of_cells, global_solution_vector_type &initial_solution, double &x_max, double &x_min,
//                   double mf, double gamma, double dx){
// x_min = 0.0;
// x_max = 1000;
// number_of_cells = x_max/dx;
//
// initial_solution.resize(number_of_cells);
//
// // double dx = (x_max-x_min) / number_of_cells;
//   for (int i = 0; i < number_of_cells; ++i) {
//     if (i < number_of_cells * 0.5) {
//       initial_solution[i] << 1.0,
//       1.0 * 1.0,
//       1.0 / (gamma * mf * mf) /
//       (gamma - 1.0) + 1.0 * 1.0 * 1.0 * 0.5,
//       1.0 * 1.0;
//     } else {
//       initial_solution[i] << 1.0,
//       1.0 * 1.0,
//       1.0 / (gamma * mf * mf) /
//       (gamma - 1.0) + 1.0 * 1.0 * 1.0 * 0.5,
//       1.0 * 1.0;
//     }
//     // initial_solution[i] << cos(x)+10, (cos(x)+10)*(cos(x)+10), cos(x)+10000, (cos(x)+10)*(cos(x)+10);
//   }
//   // for (int i = 1; i < number_of_cells-1; ++i) {
//   //   double x = x_min + (i+0.5)*dx;
//   //   initial_solution[i] << 0.5*cos(x)+10, 0.5*(cos(x)+10)*(cos(x)+10), 0.5*cos(x)+10000, 0.5*(cos(x)+10)*(cos(x)+10);
//   // }
// }
//
//
template <typename flow_properties_type, typename grid_type>
void manufactured_solution(flow_properties_type& flow, grid_type& grid){
grid.x_min = 0.0;
grid.x_max = 2.0*atan(1.0)*4.0;
// grid.number_of_cells() = ;
// std::cout << dx << std::endl;
grid.global_solution_vector.resize(grid.number_of_cells());
// double dx = (x_max-x_min) / number_of_cells;
  for (size_t i = 0; i < grid.number_of_cells(); ++i) {
    double x = (i+0.5)*grid.dx();
    grid.global_solution_vector[i] << (-0.45*tanh(4.0 * x - 10.0) + 0.55),
                          (-0.45*tanh(4.0 * x - 10.0) + 0.55)*(4.5*tanh(4.0 * x - 10.0) + 5.5),
                          2.0*tanh(4.0*x - 10.0) + 70000,
                          (-0.45*tanh(4.0 * x - 10.0) + 0.55)*(-0.5*tanh(x - 8.0/4.0) + 0.5);
  }
}
//
// void case_1(double &final_time, int number_of_cells, global_solution_vector_type &initial_solution, double gamma, double &x_max, double &x_min) {
//   // // ------------ Problem 1 ----------- //
//   double discontinuaty_location = 2.0; //  m
//   double rho_left   = 2.281;   //  kg/m^3
//   double rho_right  = 1.408;   //  kg/m^3
//   double u_left     = 164.83;  //  m/s
//   double u_right    = 0.0;     //  m/s
//   double p_left     = 201170;  //  Pa
//   double p_right    = 101100;   // Pa
//   final_time = 12.0e-3; //  s
//   x_max = 10.0;
//   for (int i = 0; i < number_of_cells; ++i) {
//     if ((i+0.5) * (x_max - x_min) / number_of_cells < discontinuaty_location) {
//       initial_solution[i] << rho_left, rho_left*u_left, p_left/(gamma-1) + rho_left * pow(u_left, 2) / 2, 0;
//     } else {
//       initial_solution[i] << rho_right, rho_right*u_right, p_right/(gamma-1) + rho_right * pow(u_right, 2) / 2, 0;
//     }
//   }
// }
// //
// void case_2(double &final_time, int number_of_cells, global_solution_vector_type &initial_solution, double gamma, double &x_max, double &x_min) {
//     // ------------ Problem 2 ----------- //
//   double discontinuaty_location = 2.0; //  m
//   double rho_left   = 1.045;   //  kg/m^3
//   double rho_right  = 3.483;   //  kg/m^3
//   double u_left     = 200.0;   //  m/s
//   double u_right    = 200.0;   //  m/s
//   double p_left     = 300000;  //  Pa
//   double p_right    = 300000;  // Pa
//   final_time = 25.0e-3; //  s
//   x_max = 10.0;
//   for (int i = 0; i < number_of_cells; ++i) {
//     if ((i+0.5) * (x_max - x_min) / number_of_cells < discontinuaty_location) {
//       initial_solution[i] << rho_left, rho_left*u_left, p_left/(gamma-1) + rho_left * pow(u_left, 2) / 2, 0;
//     } else {
//       initial_solution[i] << rho_right, rho_right*u_right, p_right/(gamma-1) + rho_right * pow(u_right, 2) / 2, 0;
//     }
//   }
// }
// //
// void case_3(double &final_time, int number_of_cells, global_solution_vector_type &initial_solution, double gamma, double &x_max, double &x_min) {
//   // ------------ Problem 3 ----------- //
//   double discontinuaty_location = 5.0; //  m
//   double rho_left   = 1.598;   //  kg/m^3
//   double rho_right  = 2.787;   //  kg/m^3
//   double u_left     = -383.64;  //  m/s
//   double u_right    = -216.97;     //  m/s
//   double p_left     = 91880;  //  Pa
//   double p_right    = 200000;   // Pa
//   final_time = 35.0e-3; //  s
//   x_max = 10.0;
//   for (int i = 0; i < number_of_cells; ++i) {
//     if ((i+0.5) * (x_max - x_min) / number_of_cells < discontinuaty_location) {
//       initial_solution[i] << rho_left, rho_left*u_left, p_left/(gamma-1) + rho_left * pow(u_left, 2) / 2, 0;
//     } else {
//       initial_solution[i] << rho_right, rho_right*u_right, p_right/(gamma-1) + rho_right * pow(u_right, 2) / 2, 0;
//     }
//   }
// }
//
// void case_4(double &final_time, int number_of_cells, global_solution_vector_type &initial_solution, double gamma, double &x_max, double &x_min) {
//   // ------------ Problem 4 ----------- //
//   double discontinuaty_location = 5.0; //  m
//   double rho_left   = 4.696;   //  kg/m^3
//   double rho_right  = 1.408;   //  kg/m^3
//   double u_left     = 0.0;  //  m/s
//   double u_right    = 0.0;     //  m/s
//   double p_left     = 404400;  //  Pa
//   double p_right    = 101100;   // Pa
//   final_time = 7.0e-3; //  s
//   x_max = 10.0;
//   for (int i = 0; i < number_of_cells; ++i) {
//     if ((i+0.5) * (x_max - x_min) / number_of_cells <= discontinuaty_location) {
//       initial_solution[i] << rho_left, rho_left*u_left, p_left/(gamma-1) + rho_left * pow(u_left, 2) / 2.0,0;
//     } else {
//       initial_solution[i] << rho_right, rho_right*u_right, p_right/(gamma-1) + rho_right * pow(u_right, 2) / 2.0,0;
//     }
//   }
// }
//
//
#endif //#ifndef INITIAL_CONDITIONS_H
