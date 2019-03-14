#ifndef HANDLE_ITTERATIVE_LAMBDA_H
#define HANDLE_ITTERATIVE_LAMBDA_H


void bisection_lambda(double& lambda_min, double& lambda_max, double& lambda_run, bool check) {
  if (check == 1){
    lambda_min = lambda_run;
    lambda_run = (lambda_min + lambda_max)*0.5;
  } else {
    lambda_max = lambda_run;
    lambda_run = (lambda_min + lambda_max)*0.5;
  }
}

void add_lambda_gap(bool& check, bool& old_check1, bool& old_check2, bool& old_check3,
                    double& lambda_min, double& lambda_max) {
  if(check == old_check1 && check == old_check2 && check == old_check3){
    if (check == 0) {
      lambda_min -= (lambda_max-lambda_min);
      std::cout << "added min" << std::endl;
    } else {
      lambda_max += (lambda_max-lambda_min);
      std::cout << "added max" << std::endl;
    }
  }
  old_check3 = old_check2;
  old_check2 = old_check1;
  old_check1 = check;
}

#endif //#ifndef HANDLE_ITTERATIVE_LAMBDA_H
