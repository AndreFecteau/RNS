#include "RK4_High_Mach_Solver_Backwards.h"
// #include "RK4_High_Mach_Solver.h"

int main(){
  RK4_High_Mach_Solver<long double>(1, 9, 5*(1+9)*(1+9)/9, 1.4, 0.14, 0.75);
}
