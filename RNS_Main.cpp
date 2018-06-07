#define HYPERBOLIC
#define VISCOUS
#define SOURCE
#define RECENTER_FLAME
#define RIGHT_CST_EXTR
// #define LEFT_CST_EXTR
// #define MANUFACTURED

#include <iomanip>
#include <fenv.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
// #include "Solver/Solver.h"
// #include "Low_Mach_Solver/RK4_Low_Mach_Solver.h"
// #include "Usefull_Headers/Variable_Vector_Isolator.h"
// #include "Gnuplot_RNS/Gnuplot_Primitive_Variables.h"
// #include "Gnuplot_RNS/Gnuplot_Primitive_Variables_Reduced.h"
// #include "Solver/Implicit_Marching.h"
// #include "Usefull_Headers/Handle_Itterative_Lambda.h"
// #include "Physical_Property/Non_Dimensional_Navier_Stokes.h"
// #include "Grid/Grid1D.h"
// #include "Implicit_Flux_and_Sources/Implicit_Centered_Difference_2nd_Order.h"
// #include "Implicit_Flux_and_Sources/Implicit_HLLE.h"
// #include "Initial_Condition/Initial_Conditions.h"
// #include "Usefull_Headers/Read_from_file.h"
#include "ChaiScript/ChaiScript.h"

int main(int argc, char* argv[]) {

  /////////////////////////////////////////////////////////////////////
  // Enable floating-point exceptions
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  /////////////////////////////////////////////////////////////////////
  // Parse command line arguments.
  // std::string chaiscriptFilename(argv[1]);
  // if(( access( chaiscriptFilename.c_str(), F_OK ) != -1 )) {
    execute_chaiscript_file("test.chai");
    // execute_chaiscript_file(chaiscriptFilename);
  // }
};
