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

void add_lambda_gap(void& check, void& old_check1, void& old_check2, void& old_check3,
                    double& lambda_min, double& lambda_max) {
  if(old_check1 == old_check2 == old_check3 == check){
    if (check = 0) {
      lambda_min -= (lambda_max-lambda_min)*0.3;
    } else {
      lambda_max += (lambda_max-lambda_min)*0.3;
    }
  }
  old_check3 = old_check2;
  old_check2 = old_check1;
  old_check1 = check;
}

#endif //#ifndef HANDLE_ITTERATIVE_LAMBDA_H
