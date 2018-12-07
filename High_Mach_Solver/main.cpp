#include "RK4_High_Mach_Solver_Backwards.h"
#include "RK4_High_Mach_Solver.h"

int main(){
  // auto a = RK4_High_Mach_Solver<long double>(1, 9, 5*(1+9)*(1+9)/9, 1.4, 0.14, 0.75);
  auto b = RK4_High_Mach_Solver_Backwards<double>(1, 9, 5*(1+9)*(1+9)/9, 1.4, 0.148777, 0.75, 17);
  // a.get_U(a.length())
}
