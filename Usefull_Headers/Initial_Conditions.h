#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H

template <typename T>
std::string tostring(T name) {
  return std::to_string(static_cast<int>(name));
}

solution_vector_type make_RK4_solution_vector(RK4_Low_Mach_Solver low_mach_solution, double x, double gamma, double mf) {
  solution_vector_type temp_vec;
  temp_vec <<   low_mach_solution.get_rho(x),
                low_mach_solution.get_rho(x) * low_mach_solution.get_U(x),
                low_mach_solution.get_rho(x) * low_mach_solution.get_T(x) /(gamma*mf*mf) /
                (gamma - 1.0) +
                low_mach_solution.get_rho(x) * low_mach_solution.get_U(x) *
                low_mach_solution.get_U(x) * 0.5,
                low_mach_solution.get_rho(x) * low_mach_solution.get_Y(x);
  return temp_vec;
}

// void RK4_low_mach_initial_conditions(double &lambda, int &number_of_cells,
//                                      global_solution_vector_type &initial_solution,
//                                      double Le, double Q, double theta, double T_ignition,
//                                      double gamma, double &x_max, double mf, double dx,
//                                      double domaine_length) {
//   RK4_Low_Mach_Solver initial_low_mach = RK4_Low_Mach_Solver(Le, Q, theta, T_ignition);
//   lambda = initial_low_mach.get_lambda();
//
//   // double domaine_length = 2000.0;
//   // double space_in_front = 1300.0;
//   // double space_in_back  = domaine_length - space_in_front - 1.0;
//   x_max = domaine_length;
//   number_of_cells = x_max/dx;
//   initial_solution.resize(number_of_cells);
// // #pragma omp parallel for
//   for (int i = 1; i <= number_of_cells; ++i) {
//     if (i*dx < domaine_length*0.5) {
//       initial_solution[number_of_cells - i] << 1.0,
//       1.0 * 1.0,
//       (1.0 / (gamma * mf * mf)) /
//       (gamma - 1.0) + 1.0 * 1.0 * 1.0 * 0.5,
//       1.0 * 1.0;
//     } else if ((i+1)*dx > domaine_length*0.5 + initial_low_mach.length()) {
//       initial_solution[number_of_cells - i] << initial_solution[number_of_cells-(i-1)];
//   // std::cout << "back" << std::endl;
//     } else {
//       initial_solution[number_of_cells - i] << make_RK4_solution_vector(initial_low_mach,
//                                                       (i+1)*dx - domaine_length*0.5, gamma, mf);
//     }
//     // std::cout << "sol" << std::endl;
//     initial_solution[number_of_cells-i][1] = -fabs(initial_solution[number_of_cells-i][1]);
//   }
//   // std::cout << "here" << std::endl;
// }
//
void RK4_low_mach_initial_conditions(double &lambda, int &number_of_cells,
                                     global_solution_vector_type &initial_solution,
                                     double Le, double Q, double theta, double T_ignition,
                                     double gamma, double &x_max, double mf, double dx,
                                     double domaine_length) {
  RK4_Low_Mach_Solver initial_low_mach = RK4_Low_Mach_Solver(Le, Q, theta, T_ignition);
  lambda = initial_low_mach.get_lambda();

  // double domaine_length = 2000.0;
  // double space_in_front = 1300.0;
  // double space_in_back  = domaine_length - space_in_front - 1.0;
  x_max = domaine_length;
  number_of_cells = x_max/dx;
  initial_solution.resize(number_of_cells);
// #pragma omp parallel for
  for (int i = 0; i < number_of_cells; ++i) {
    if (i*dx < domaine_length*0.7) {
      initial_solution[i] << 1.0,
      1.0 * 1.0,
      (1.0 / (gamma * mf * mf)) /
      (gamma - 1.0) + 1.0 * 1.0 * 1.0 * 0.5,
      1.0 * 1.0;
    } else if ((i+1)*dx > domaine_length*0.7 + initial_low_mach.length()) {
      initial_solution[i] << initial_solution[(i-1)];
  // std::cout << "back" << std::endl;
    } else {
      initial_solution[i] << make_RK4_solution_vector(initial_low_mach,
                                                      (i+1)*dx - domaine_length*0.7, gamma, mf);
    }
    // std::cout << "sol" << std::endl;
    // initial_solution[i][1] = -fabs(initial_solution[number_of_cells-i][1]);
  }
  // std::cout << "here" << std::endl;
}

void straight_line(int number_of_cells, global_solution_vector_type &initial_solution, double &x_max, double &x_min,
                  double mf, double gamma, double dx){
x_min = 0.0;
x_max = 1000;
number_of_cells = x_max/dx;

initial_solution.resize(number_of_cells);

// double dx = (x_max-x_min) / number_of_cells;
  for (int i = 0; i < number_of_cells; ++i) {
    if (i < number_of_cells * 0.5) {
      initial_solution[i] << 1.0,
      1.0 * 1.0,
      1.0 / (gamma * mf * mf) /
      (gamma - 1.0) + 1.0 * 1.0 * 1.0 * 0.5,
      1.0 * 1.0;
    } else {
      initial_solution[i] << 1.0,
      1.0 * 1.0,
      1.0 / (gamma * mf * mf) /
      (gamma - 1.0) + 1.0 * 1.0 * 1.0 * 0.5,
      1.0 * 1.0;
    }
    // initial_solution[i] << cos(x)+10, (cos(x)+10)*(cos(x)+10), cos(x)+10000, (cos(x)+10)*(cos(x)+10);
  }
  // for (int i = 1; i < number_of_cells-1; ++i) {
  //   double x = x_min + (i+0.5)*dx;
  //   initial_solution[i] << 0.5*cos(x)+10, 0.5*(cos(x)+10)*(cos(x)+10), 0.5*cos(x)+10000, 0.5*(cos(x)+10)*(cos(x)+10);
  // }
}


void manufactured_solution(int &number_of_cells, global_solution_vector_type &initial_solution, double &x_max, double &x_min, double dx){
x_min = 0.0;
x_max = 2.0*atan(1.0)*4.0;
number_of_cells = x_max/dx;
std::cout << dx << std::endl;
initial_solution.resize(number_of_cells);
// double dx = (x_max-x_min) / number_of_cells;
  for (size_t i = 0; i < number_of_cells; ++i) {
    double x = (i+0.5)*dx;
    initial_solution[i] << (-0.45*tanh(4.0 * x - 10.0) + 0.55),
                          (-0.45*tanh(4.0 * x - 10.0) + 0.55)*(4.5*tanh(4.0 * x - 10.0) + 5.5),
                          2.0*tanh(4.0*x - 10.0) + 70000,
                          (-0.45*tanh(4.0 * x - 10.0) + 0.55)*(-0.5*tanh(x - 8.0/4.0) + 0.5);
    // initial_solution[i] << cos(x)+10, (cos(x)+10)*(cos(x)+10), cos(x)+10000, (cos(x)+10)*(cos(x)+10);
  }
  // for (int i = 1; i < number_of_cells-1; ++i) {
  //   double x = x_min + (i+0.5)*dx;
  //   initial_solution[i] << 0.5*cos(x)+10, 0.5*(cos(x)+10)*(cos(x)+10), 0.5*cos(x)+10000, 0.5*(cos(x)+10)*(cos(x)+10);
  // }
}

void case_1(double &final_time, int number_of_cells, global_solution_vector_type &initial_solution, double gamma, double &x_max, double &x_min) {
  // // ------------ Problem 1 ----------- //
  double discontinuaty_location = 2.0; //  m
  double rho_left   = 2.281;   //  kg/m^3
  double rho_right  = 1.408;   //  kg/m^3
  double u_left     = 164.83;  //  m/s
  double u_right    = 0.0;     //  m/s
  double p_left     = 201170;  //  Pa
  double p_right    = 101100;   // Pa
  final_time = 12.0e-3; //  s
  x_max = 10.0;
  for (int i = 0; i < number_of_cells; ++i) {
    if ((i+0.5) * (x_max - x_min) / number_of_cells < discontinuaty_location) {
      initial_solution[i] << rho_left, rho_left*u_left, p_left/(gamma-1) + rho_left * pow(u_left, 2) / 2, 0;
    } else {
      initial_solution[i] << rho_right, rho_right*u_right, p_right/(gamma-1) + rho_right * pow(u_right, 2) / 2, 0;
    }
  }
}
//
void case_2(double &final_time, int number_of_cells, global_solution_vector_type &initial_solution, double gamma, double &x_max, double &x_min) {
    // ------------ Problem 2 ----------- //
  double discontinuaty_location = 2.0; //  m
  double rho_left   = 1.045;   //  kg/m^3
  double rho_right  = 3.483;   //  kg/m^3
  double u_left     = 200.0;   //  m/s
  double u_right    = 200.0;   //  m/s
  double p_left     = 300000;  //  Pa
  double p_right    = 300000;  // Pa
  final_time = 25.0e-3; //  s
  x_max = 10.0;
  for (int i = 0; i < number_of_cells; ++i) {
    if ((i+0.5) * (x_max - x_min) / number_of_cells < discontinuaty_location) {
      initial_solution[i] << rho_left, rho_left*u_left, p_left/(gamma-1) + rho_left * pow(u_left, 2) / 2, 0;
    } else {
      initial_solution[i] << rho_right, rho_right*u_right, p_right/(gamma-1) + rho_right * pow(u_right, 2) / 2, 0;
    }
  }
}
//
void case_3(double &final_time, int number_of_cells, global_solution_vector_type &initial_solution, double gamma, double &x_max, double &x_min) {
  // ------------ Problem 3 ----------- //
  double discontinuaty_location = 5.0; //  m
  double rho_left   = 1.598;   //  kg/m^3
  double rho_right  = 2.787;   //  kg/m^3
  double u_left     = -383.64;  //  m/s
  double u_right    = -216.97;     //  m/s
  double p_left     = 91880;  //  Pa
  double p_right    = 200000;   // Pa
  final_time = 35.0e-3; //  s
  x_max = 10.0;
  for (int i = 0; i < number_of_cells; ++i) {
    if ((i+0.5) * (x_max - x_min) / number_of_cells < discontinuaty_location) {
      initial_solution[i] << rho_left, rho_left*u_left, p_left/(gamma-1) + rho_left * pow(u_left, 2) / 2, 0;
    } else {
      initial_solution[i] << rho_right, rho_right*u_right, p_right/(gamma-1) + rho_right * pow(u_right, 2) / 2, 0;
    }
  }
}

void case_4(double &final_time, int number_of_cells, global_solution_vector_type &initial_solution, double gamma, double &x_max, double &x_min) {
  // ------------ Problem 4 ----------- //
  double discontinuaty_location = 5.0; //  m
  double rho_left   = 4.696;   //  kg/m^3
  double rho_right  = 1.408;   //  kg/m^3
  double u_left     = 0.0;  //  m/s
  double u_right    = 0.0;     //  m/s
  double p_left     = 404400;  //  Pa
  double p_right    = 101100;   // Pa
  final_time = 7.0e-3; //  s
  x_max = 10.0;
  for (int i = 0; i < number_of_cells; ++i) {
    if ((i+0.5) * (x_max - x_min) / number_of_cells <= discontinuaty_location) {
      initial_solution[i] << rho_left, rho_left*u_left, p_left/(gamma-1) + rho_left * pow(u_left, 2) / 2.0,0;
    } else {
      initial_solution[i] << rho_right, rho_right*u_right, p_right/(gamma-1) + rho_right * pow(u_right, 2) / 2.0,0;
    }
  }
}


#endif //#ifndef INITIAL_CONDITIONS_H
