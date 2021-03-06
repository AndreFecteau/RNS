#include <iomanip>
#include <iostream>
#include "Text_Output.h"
#include "ChaiScript/ChaiScript.h"

int main(int argc, char* argv[]) {

  /////////////////////////////////////////////////////////////////////
  // Im pretty now
  loading_text();
  /////////////////////////////////////////////////////////////////////
  // Enable floating-point exceptions (I cant do this due to the unkown
  // stable timestep i am trying and if fail re-starting)
  // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  std::cout << std::setprecision(20);
  /////////////////////////////////////////////////////////////////////
  // Parse command line arguments.
  if(argc > 2 && static_cast<std::string>(argv[2]) == "-omp"){
    std::string str_num(argv[3]);
    if(!std::all_of(str_num.begin(), str_num.end(), ::isdigit) || std::stoi(str_num) == 0) {
      throw std::invalid_argument("Error setting number of OpenMP threads, \""
      + str_num +
      "\" is not a valid positive integer.");
    }
    omp_set_num_threads(std::stoi(str_num));
  } else {
      omp_set_num_threads(1);
  }
  execute_chaiscript_file(argv[1]);
};
