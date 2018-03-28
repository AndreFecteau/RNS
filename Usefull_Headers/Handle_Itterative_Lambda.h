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

#endif //#ifndef HANDLE_ITTERATIVE_LAMBDA_H
