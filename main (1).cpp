#include<iostream>
#include<math.h>
#include<vector>
#include<utility>
#include<memory>
#include<sstream>
#include<fstream>
#include</home/william/Include/eigen3/Eigen/Dense>


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                             ROE FLUX FUNCTION
////////////////////////////////////////////////////////////////////////////////
Eigen::Vector3d Flux_Roe(Eigen::Vector3d Ul, Eigen::Vector3d Ur){

   Eigen::Vector3d flux;
   double gamma  = 1.4;
  
   double rhol = Ul[0];
   double ul = Ul[1]/Ul[0];
   double pl = (Ul[2]-0.5*Ul[0]*ul*ul)*(gamma-1.0);
   double hl = gamma*pl/((gamma-1.0)*rhol) + ul*ul*0.5;
   double al = sqrt(gamma*pl/rhol);
   Eigen::Vector3d Fl;
   Fl[0] = rhol*ul;
   Fl[1] = rhol*ul*ul+pl;
   Fl[2] = ul*(gamma*pl/(gamma-1.0) + rhol*ul*ul*0.5);

   double rhor = Ur[0];
   double ur = Ur[1]/Ur[0];
   double pr = (Ur[2]-0.5*Ur[0]*ur*ur)*(gamma-1.0);
   double hr = gamma*pr/((gamma-1.0)*rhor) + ur*ur*0.5;
   double ar = sqrt(gamma*pr/rhor);
   Eigen::Vector3d Fr;
   Fr[0] = rhor*ur;
   Fr[1] = rhor*ur*ur+pr;
   Fr[2] = ur*(gamma*pr/(gamma-1.0) + rhor*ur*ur*0.5);
   
   double rhoh = sqrt(rhol*rhor);
   double uh = (sqrt(rhol)*ul + sqrt(rhor)*ur)/(sqrt(rhol)+sqrt(rhor));
   double hh = (sqrt(rhol)*hl + sqrt(rhor)*hr)/(sqrt(rhol)+sqrt(rhor));
   double ph = (hh - 0.5*uh*uh)*rhoh*(gamma-1.0)/gamma;
   double ah = sqrt(gamma*ph/rhoh);

   double wave_l;
   double delta_min = std::max(0.0, 4.0*((ur-ar)-(ul-al)));
   if(fabs(uh-ah) >= delta_min*0.5){
     wave_l = fabs(uh-ah);
   }else{
     wave_l = fabs(uh-ah)*fabs(uh-ah)/delta_min + 0.25*delta_min;
   }   
   
   double wave_r;
   double delta_max = std::max(0.0, 4.0*((ur+ar)-(ul+al)));
   if(fabs(uh+ah) >= delta_max*0.5){
     wave_r = fabs(uh+ah);
   }else{
     wave_r = fabs(uh+ah)*fabs(uh+ah)/delta_max + 0.25*delta_max;
   }   
   

   Eigen::Matrix3d R;
   R(0,0) = 1.0;
   R(0,1) = 1.0;
   R(0,2) = 1.0;
   R(1,0) = uh-ah;
   R(1,1) = uh;
   R(1,2) = uh+ah;
   R(2,0) = hh-uh*ah;
   R(2,1) = uh*uh*0.5;
   R(2,2) = hh+uh*ah;

   Eigen::Matrix3d R_inv = R.inverse();
   Eigen::Matrix3d Omega;
   Omega(0,0) = wave_l;
   Omega(0,1) = 0.0;
   Omega(0,2) = 0.0;
   Omega(1,0) = 0.0;
   Omega(1,1) = fabs(uh);
   Omega(1,2) = 0.0;
   Omega(2,0) = 0.0;
   Omega(2,1) = 0.0;
   Omega(2,2) = wave_r;

   flux = (Fl+Fr)*0.5 - 0.5*R*Omega*R_inv*(Ur-Ul);

  return flux;
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                             HLLE FLUX FUNCTION
////////////////////////////////////////////////////////////////////////////////
Eigen::Vector3d Flux_HLLE(Eigen::Vector3d Ul, Eigen::Vector3d Ur){

   Eigen::Vector3d flux;
   double gamma  = 1.4;
  
   double rhol = Ul[0];
   double ul = Ul[1]/Ul[0];
   double pl = (Ul[2]-0.5*Ul[0]*ul*ul)*(gamma-1.0);
   double hl = gamma*pl/((gamma-1.0)*rhol) + ul*ul*0.5;
   double al = sqrt(gamma*pl/rhol);
   Eigen::Vector3d Fl;
   Fl[0] = rhol*ul;
   Fl[1] = rhol*ul*ul+pl;
   Fl[2] = ul*(gamma*pl/(gamma-1.0) + rhol*ul*ul*0.5);

   double rhor = Ur[0];
   double ur = Ur[1]/Ur[0];
   double pr = (Ur[2]-0.5*Ur[0]*ur*ur)*(gamma-1.0);
   double hr = gamma*pr/((gamma-1.0)*rhor) + ur*ur*0.5;
   double ar = sqrt(gamma*pr/rhor);
   Eigen::Vector3d Fr;
   Fr[0] = rhor*ur;
   Fr[1] = rhor*ur*ur+pr;
   Fr[2] = ur*(gamma*pr/(gamma-1.0) + rhor*ur*ur*0.5);
   
   double rhoh = sqrt(rhol*rhor);
   double uh = (sqrt(rhol)*ul + sqrt(rhor)*ur)/(sqrt(rhol)+sqrt(rhor));
   double hh = (sqrt(rhol)*hl + sqrt(rhor)*hr)/(sqrt(rhol)+sqrt(rhor));
   double ph = (hh - 0.5*uh*uh)*rhoh*(gamma-1.0)/gamma;
   double ah = sqrt(gamma*ph/rhoh);

   double lambda_l = std::min(ul-al, uh-ah);
   double lambda_r = std::max(ur + ar, uh+ah);

   if(lambda_l > 0.0){
     flux = Fl;
   }else if(lambda_r < 0.0){
     flux =  Fr;
   }else{  
     flux = (lambda_r*Fl - lambda_l*Fr)/(lambda_r-lambda_l) + lambda_r*lambda_l*(Ur-Ul)/(lambda_r-lambda_l);
   }

  return flux;
}
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//              Limiter
/////////////////////////////////////////////////////////////////////////////////////
auto Limiter(std::vector<Eigen::Vector3d> U, double num_cells, double delta_x){

  std::vector<Eigen::Vector3d> Phi(num_cells);
  std::vector<Eigen::Vector3d> dUdx(num_cells);
  std::vector<Eigen::Vector3d> a(num_cells);
  std::vector<Eigen::Vector3d> b(num_cells);
  double e = 10e-6;

  
  for(int i = 0; i < num_cells; i++){

      if(i == 0){
          a[i] = (U[i]-U[i])/delta_x;
          b[i] = (U[i+1]-U[i])/delta_x;
      }else if(i == num_cells-1){
          a[i] = (U[i]-U[i-1])/delta_x;
          b[i] = (U[i]-U[i])/delta_x;
      }else{	
          a[i] = (U[i]-U[i-1])/delta_x;
          b[i] = (U[i+1]-U[i])/delta_x;
      }
      
   }// for
	

  for(int i = 0; i < num_cells; i++){
      
      for(int j = 0; j < 3; j++){	
          if(a[i][j]*b[i][j] < 0.0){
	    Phi[i][j] = 0.0;
	  }else{
            Phi[i][j] = 2.0*a[i][j]*b[i][j]/(a[i][j]*a[i][j]+b[i][j]*b[i][j]+e);
	  }  
      }

      dUdx[i] = 0.5*(a[i]+b[i]);
	  	  
   }// end for

   return std::make_tuple(Phi,dUdx);     
	
}
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//              MAIN
/////////////////////////////////////////////////////////////////////////////////////

 
int main() {
 
 double  CFL = 0.5;
 int num_nodes = 1000;
 int num_cells = num_nodes-1;
 double delta_x = 10.0/num_cells;

 // std::vector<double> C = {1.4, 2.0, 0.012, 2.281, 164.83, 201170, 1.408, 0.0, 101100};
 // std::vector<double> C = {1.4, 2.0, 0.025, 1.045, 200.0, 300000, 3.483, 200.0, 300000};
 // std::vector<double> C = {1.4, 5.0, 0.035, 1.598, -383.64, 91880, 2.787, -216.97, 200000};
  std::vector<double> C = {1.4, 5.0, 0.007, 4.696, 0.0, 404400, 1.408, 0.0, 101100};
 
 double gamma = C[0];
 double origin = C[1];
 double final_time = C[2];
 
 double rhol = C[3];
 double ul = C[4];
 double pl = C[5];
 
 double rhor = C[6];
 double ur = C[7];
 double pr = C[8];

 
 std::vector<Eigen::Vector3d> U(num_cells);
 std::vector<Eigen::Vector3d> U_hat(num_cells);
 std::vector<Eigen::Vector3d> dUdt(num_cells);
 std::vector<Eigen::Vector3d> dUdx(num_cells);
 std::vector<Eigen::Vector3d> Phi(num_cells);
 std::vector<double> X(num_cells);
 
 Eigen::Vector3d Ul;
 Eigen::Vector3d Ur;
 Eigen::Vector3d Fl;
 Eigen::Vector3d Fr;
 Eigen::Vector3d flux;
 
 
 
 for(int i = 0; i < num_cells; i++){
   X[i] = (i+0.5)*delta_x;
 }  
 
 
 for(int i = 0; i <num_cells; i++){
   if(X[i] >= 0.0 && X[i] <= origin){
      U[i][0] = rhol;
      U[i][1] = rhol*ul;
      U[i][2] = pl/(gamma-1.0) + rhol*ul*ul*0.5;
  }else{
      U[i][0] = rhor;
      U[i][1] = rhor*ur;
      U[i][2] = pr/(gamma-1.0) + rhor*ur*ur*0.5;
  }
 }
 
 
/////////////////////////////////////////////////////////////////////////////
//////       TIME MARCHING
/////////////////////////////////////////////////////////////////////////////
 
  double time = 0.0;
 
  while(time < final_time){
 
     double max_lambda = 0.0;
     double dummy_lambda = 0.0;
     
     for( int i = 0; i < num_cells; i++){
       dummy_lambda = fabs(U[i][1]/U[i][0]) + sqrt(gamma*(U[i][2]-0.5*U[i][0]*(U[i][1]/U[i][0])*(U[i][1]/U[i][0])*(gamma-1.0))/ U[i][0]);
         if(dummy_lambda > max_lambda){
             max_lambda = dummy_lambda;
         }
     }   
 		   
 
     double delta_t = CFL*delta_x/max_lambda; 
 
     if(time + delta_t > final_time){
        delta_t = final_time - time;
     }
	
 
 
     for(int i = 0; i < num_cells; i++){
           dUdt[i][0] = 0.0;
           dUdt[i][1] = 0.0;
           dUdt[i][2] = 0.0;
     }    

     std::tie(Phi,dUdx) = Limiter(U, num_cells, delta_x);
     
     for(int i = 0; i <= num_cells; i++){
 
         if(i == 0){
 	  Ul = U[i];
  	  Ur = U[i];
         }else if(i == num_cells){
 	  Ul = U[i-1];
 	  Ur = U[i-1];
         }else{
	   Ul = U[i-1] + dUdx[i-1].cwiseProduct(Phi[i-1])*0.5*delta_x;
	   Ur = U[i] - dUdx[i].cwiseProduct(Phi[i])*0.5*delta_x;
	 }  

	 flux = Flux_Roe(Ul, Ur);
	 
         if(i == 0){
 	  dUdt[i] += flux/delta_x;
         }else if(i == num_cells){
	  dUdt[i-1] -= flux/delta_x;
	 }else{  
 	  dUdt[i] += flux/delta_x;
 	  dUdt[i-1] -= flux/delta_x;
         }
 	
      }
 	
      for(int i = 0; i < num_cells; i++){ 
          U_hat[i] = U[i]+dUdt[i]*delta_t;
      }


      std::tie(Phi,dUdx) = Limiter(U_hat, num_cells, delta_x);
      

      for(int i = 0; i <= num_cells; i++){
 
         if(i == 0){
 	  Ul = U_hat[i];
  	  Ur = U_hat[i];
         }else if(i == num_cells){
 	  Ul = U_hat[i-1];
 	  Ur = U_hat[i-1];
         }else{
	   Ul = U_hat[i-1] + dUdx[i-1].cwiseProduct(Phi[i-1])*0.5*delta_x;
	   Ur = U_hat[i] - dUdx[i].cwiseProduct(Phi[i])*0.5*delta_x;
	 }  

	 flux = Flux_Roe(Ul, Ur);
	 
         if(i == 0){
 	  dUdt[i] += flux/delta_x;
         }else if(i == num_cells){
	  dUdt[i-1] -= flux/delta_x;
	 }else{  
 	  dUdt[i] += flux/delta_x;
 	  dUdt[i-1] -= flux/delta_x;
         }
 	
      }
 	
      for(int i = 0; i < num_cells; i++){ 
          U[i] = U[i]+0.5*dUdt[i]*delta_t;
      }
      time += delta_t;
 
 
   }// end while loop 
 
 
   
   
  std::ofstream fout;
  fout.open("A5_4.txt");


  for(int i = 0; i < num_cells; i++){
    fout <<X[i]<<" "<<U[i][0]<<" "<<U[i][1]/U[i][0]<<" "<<gamma*(U[i][2]-0.5*U[i][0]*(U[i][1]/U[i][0])*(U[i][1]/U[i][0]))*(gamma-1.0)<<std::endl;
  }

  fout.close();
 
 

 
  
return 0;

} // end main
