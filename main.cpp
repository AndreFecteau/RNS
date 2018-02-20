#include "/usr/include/eigen3/Eigen/Core"
#include "/usr/include/eigen3/Eigen/Dense"
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <iomanip>

typedef Eigen::Matrix<float, -1, 1> Vector_type;
// using Vector_type = Eigen::VectorXd;
using global_solution_vector_type = std::vector<Vector_type>;

template <typename T>
std::string tostring(T name) {
  return std::to_string(static_cast<int>(name));
}

int main(){

double per_FL = 16;
double domaine_length = 4000;
double CFL = 1e4;
int location = 0;

  std::string filename = "dat_saves/Implicit_CD_16_Domaine/Plot_" + tostring(domaine_length) + "_"
                                       + tostring(std::log10(CFL)) + "_2";

  global_solution_vector_type global_solution_vector;
  std::string line;
  std::ifstream gnu_input_file(filename + ".dat");
  // gnu_input_file.open(filename + ".dat", std::ios_base::app);
  if (gnu_input_file) {
    std::string first_line;
    std::getline(gnu_input_file,first_line);
    std::ofstream gnu_output_file(filename + "_N.dat");
    while(std::getline(gnu_input_file,line))
    {
      std::istringstream iss(line);
      std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                       std::istream_iterator<std::string>( ));
      Vector_type temp_vector;
      temp_vector.resize(results.size());
      for(int i = 0; i < results.size(); ++i) {
        // std::cout << results[i] << std::endl;
        temp_vector[i] = std::stold(results[i]);
      }
      // std::cout << std::setprecision(19) << temp_vector << std::endl;
      global_solution_vector.push_back(temp_vector);
    }
    int i = 0;
    while(location == 0 && i < global_solution_vector.size()) {
      // std::cout << global_solution_vector[i][1]<< std::endl;
      if(global_solution_vector[i][1] < 0.5) {
        location = i;
      } else {
      ++i;
    }
    }
    double temp = global_solution_vector[location][0] - (global_solution_vector[location][0]-global_solution_vector[location-1][0])/(global_solution_vector[location][1]-global_solution_vector[location-1][1])*((global_solution_vector[location][1]-0.5));
    for(int j = 0; j < global_solution_vector.size(); ++j){
      // std::cout << j << std::endl;
      global_solution_vector[j][0] -= temp;
    }
    gnu_output_file << first_line << std::endl;
    for(int j = 0; j < global_solution_vector.size(); ++j){
      // for (int n = 0; n < global_solution_vector[0].size(); ++n{
        gnu_output_file << global_solution_vector[j][0] << " " <<
                           global_solution_vector[j][1] << " " <<
                           global_solution_vector[j][2] << " " <<
                           global_solution_vector[j][3] << " " <<
                           global_solution_vector[j][4] << " " <<
                           global_solution_vector[j][5] << " " <<
                           global_solution_vector[j][6] << std::endl;
      // }
    }

    // std::cout << location << std::endl;
  } else {
    std::cout << "not a file: " << filename << ".dat" << std::endl;
  }
};
