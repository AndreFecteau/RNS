#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H

template <typename T>
std::string tostring(T name) {
  return std::to_string(static_cast<int>(name));
}

void bisection_lambda(double& lambda_min, double& lambda_max, double& lambda_run, double& T_sum_initial, double& T_sum) {
  if (T_sum < T_sum_initial){
    lambda_min = lambda_run;
    lambda_run = (lambda_min + lambda_max)*0.5;
  } else {
    lambda_max = lambda_run;
    lambda_run = (lambda_min + lambda_max)*0.5;
  }
}

double flame_position_algorithm(global_solution_vector_type global_solution_vector, double gamma) {
  double sum = 0.0;
  for (size_t i = 0; i < global_solution_vector.size(); ++i){
    Variable_Vector_Isolator<solution_vector_type> var_vec = Variable_Vector_Isolator<solution_vector_type>(global_solution_vector[i], gamma);
    sum += var_vec.T();
  }
  return sum;
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

void RK4_low_mach(double &lambda, int number_of_cells, global_solution_vector_type &initial_solution,
  double Le, double Q, double theta, double T_ignition, double gamma, double &x_max, double mf) {
  RK4_Low_Mach_Solver initial_low_mach = RK4_Low_Mach_Solver(Le, Q, theta, T_ignition);
  lambda = initial_low_mach.get_lambda();
  double domaine_length = 3.0;
  double space_in_front = 1.45;
  double space_in_back  = domaine_length - space_in_front - 1.0;
  x_max = initial_low_mach.length() * domaine_length;
  double dx = x_max / number_of_cells;
#pragma omp parallel for
  for (int i = 0; i < number_of_cells; ++i) {
    if (i < number_of_cells * space_in_front*0.86/domaine_length) {
      initial_solution[i] << 1.0,
      1.0 * 1.0,
      1.0 / (gamma * mf * mf) /
      (gamma - 1.0) + 1.0 * 1.0 * 1.0 * 0.5,
      1.0 * 1.0;
    } else if (i < number_of_cells * space_in_front/domaine_length && i > number_of_cells * (space_in_front*0.85)/domaine_length) {
      initial_solution[i] << make_RK4_solution_vector(initial_low_mach, (0.0) * dx + 0.5*dx, gamma, mf);
      // initial_solution[i] << 1.0,
      // 1.0 * 1.0,
      // 1.0 / (gamma * mf * mf) /
      // (gamma - 1.0) + 1.0 * 1.0 * 1.0 * 0.5,
      // 1.0 * 1.0 * 0.0;
    } else if (i > number_of_cells * (domaine_length - space_in_back)/domaine_length - 1) {
      initial_solution[i] << make_RK4_solution_vector(initial_low_mach, initial_low_mach.length(), gamma, mf);
    } else {
      initial_solution[i] << make_RK4_solution_vector(initial_low_mach, (i - space_in_front/domaine_length * number_of_cells) * dx + 0.5*dx, gamma, mf);
    }
  }
}

  // void RK4_low_mach2(double &lambda, int number_of_cells, global_solution_vector_type &initial_solution,
  //   double Le, double Q, double theta, double T_ignition, double gamma, double &x_max, double mf) {
  //     RK4_Low_Mach_Solver initial_low_mach = RK4_Low_Mach_Solver(Le, Q, theta, T_ignition);
  //     lambda = initial_low_mach.get_lambda();
  //     double domaine_length = 3.0;
  //     double space_in_front = 1.5;
  //     // double space_in_back  = domaine_length - space_in_front;
  //     x_max = initial_low_mach.length() * domaine_length;
  //     // double dx = x_max / number_of_cells;
  //     for (int i = 0; i < number_of_cells; ++i) {
  //       if (i < number_of_cells * space_in_front/domaine_length) {
  //         initial_solution[i] << 1.0,
  //         1.0 * 1.0,
  //         1.0 / (gamma * mf * mf) /
  //         (gamma - 1) + 1.0 * 1.0 * 1.0 * 0.5,
  //         1.0 * 1.0;
  //         // 1.0 - (1.0 - make_RK4_solution_vector(initial_low_mach, 0, gamma, mf)[3])/(space_in_front*x_max)*i*dx*domaine_length;
  //       } else {
  //         initial_solution[i] << make_RK4_solution_vector(initial_low_mach, initial_low_mach.length(), gamma, mf);
  //       // } else {
  //       //   initial_solution[i] << make_RK4_solution_vector(initial_low_mach, (i - space_in_front/domaine_length * number_of_cells) * dx + 0.5*dx, gamma, mf);
  //       // }
  //     }
  //   }
  // }


// void fixed_boundaries_evolution(double &final_time, int number_of_cells, global_solution_vector_type &initial_solution, double gamma, double &x_max, double &x_min) {
//   // // ------------ Problem 1 ----------- //
//   double discontinuaty_location = (x_max+x_min)/2; //  m
//   double rho_left   = 1.0;   //  kg/m^3
//   double rho_right  = 0.2;   //  kg/m^3
//   double u_left     = 1.0;  //  m/s
//   double u_right    = 10.0;     //  m/s
//   double p_left     = 1.0/(1.4*0.005*0.005);  //  Pa
//   double p_right    = 0.9997/(1.4*0.005*0.005);   // Pa
//   double Y_left     = 1.0;
//   double Y_right    = 0.0;
//   for (int i = 0; i < number_of_cells; ++i) {
//     if ((i+0.5) * (x_max - x_min) / number_of_cells < discontinuaty_location) {
//       initial_solution[i] << rho_left, rho_left*u_left, p_left/(gamma-1) + rho_left * pow(u_left, 2) / 2, rho_left*Y_left;
//     } else {
//       initial_solution[i] << rho_right, rho_right*u_right, p_right/(gamma-1) + rho_right * pow(u_right, 2) / 2, rho_right*Y_right;
//     }
//   }
// }
//
// void shock_evolution(double &final_time, int number_of_cells, global_solution_vector_type &initial_solution, double gamma, double &x_max, double &x_min) {
//   // // ------------ Problem 1 ----------- //
//   double discontinuaty_location = (x_max+x_min)/2; //  m
//   double Ms = 2;
//
//
//   double rho_left = 1.225;
//   double p_left   = 101325.0;
//   double c = sqrt(gamma * p_left / rho_left);
//   double u_left   = -Ms * c;
//
//   double rho_right = (gamma + 1) * Ms * Ms *rho_left / ((gamma - 1) * Ms * Ms + 2);
//   double u_right = 2 * (Ms * Ms - 1) * c / ((gamma + 1) * Ms) + u_left;
//   double p_right = 2 * gamma * (Ms * Ms - 1) * p_left / (gamma + 1) + p_left;
//
//   double Y_left     = 1.0;
//   double Y_right    = 0.0;
//   for (int i = 0; i < number_of_cells; ++i) {
//     if ((i+0.5) * (x_max - x_min) / number_of_cells < discontinuaty_location) {
//       initial_solution[i] << rho_left, rho_left*u_left, p_left/(gamma-1) + rho_left * pow(u_left, 2) / 2, rho_left*Y_left;
//     } else {
//       initial_solution[i] << rho_right, rho_right*u_right, p_right/(gamma-1) + rho_right * pow(u_right, 2) / 2, rho_right*Y_right;
//     }
//   }
// }
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
//
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
//
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
    if ((i+0.5) * (x_max - x_min) / number_of_cells < discontinuaty_location) {
      initial_solution[i] << rho_left, rho_left*u_left, p_left/(gamma-1) + rho_left * pow(u_left, 2) / 2, 0;
    } else {
      initial_solution[i] << rho_right, rho_right*u_right, p_right/(gamma-1) + rho_right * pow(u_right, 2) / 2, 0;
    }
  }
}


#endif //#ifndef INITIAL_CONDITIONS_H
