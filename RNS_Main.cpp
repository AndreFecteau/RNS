#define HYPERBOLIC
#define VISCOUS
#define SOURCE
#define RECENTER_FLAME
// #define RIGHT_CST_EXTR
// #define LEFT_CST_EXTR
// #define MANUFACTURED

#include "ChaiScript/ChaiScript.h"

int main(int argc, char* argv[]) {

  /////////////////////////////////////////////////////////////////////
  // Enable floating-point exceptions (I cant do this due to the unkown
  // stable timestep i am trying and if fail re-starting)
  // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  /////////////////////////////////////////////////////////////////////
  // Parse command line arguments.
  execute_chaiscript_file(argv[1]);
};
