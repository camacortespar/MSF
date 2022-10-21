#include <iostream>
#include <fstream>
#include <cmath>


double fun1(double x1, double x2, double t, double L);
double fun2(double x1, double x2, double t, double L);
void pasoRungeKutta4(double & x1, double & x2, double & t, double dt, double L);

int main(void)
{
  std::cout.precision(8);
  std::cout.setf(std::ios::scientific);

//x2 es R y x1 es la derivada de R  
  double t=0;
  double x1=0, x2=0;
  double dt=0.001;
  double L=1;
  
  // std::ofstream myfile("datos.dat");
  
  for(L=0.1;L<=15; L+=0.1){
    for(t=0.001, x1=0, x2=1; t<=1;){
      pasoRungeKutta4(x1, x2, t, dt, L);    
    }
    std::cout << L  <<"\t"<< x2 <<std::endl;
    //myfile << L <<"\t" << x2 <<std::endl;
  }
  //myfile.close();
 
  //plot de gnuplot
  //std::cout <<"plot \"datos.dat\" u 1:2 w l" << std::endl;
  //std::cout <<"pause 7" << std::endl;
  return 0;
}

double fun1(double x1, double x2, double t, double L) 
{
  return (-x1/t-L*L*x2);
}
double fun2(double x1, double x2, double t, double L)
{
  return x1;
}
void pasoRungeKutta4(double & x1, double & x2, double & t, double dt, double L)
{
  double dx11, dx21, dx31, dx41;                      double dx12, dx22, dx32, dx42;
  dx11=dt*fun1(x1,x2,t,L);                             dx12=dt*fun2(x1,x2,t,L);
  dx21=dt*fun1(x1+0.5*dx11,x2+0.5*dx12,t+0.5*dt,L);    dx22=dt*fun2(x1+0.5*dx11,x2+0.5*dx12,t+0.5*dt,L);
  dx31=dt*fun1(x1+0.5*dx21,x2+0.5*dx22,t+0.5*dt,L);    dx32=dt*fun2(x1+0.5*dx21,x2+0.5*dx22,t+0.5*dt,L);
  dx41=dt*fun1(x1+dx31,x2+dx32,t+dt,L);                dx42=dt*fun2(x1+dx31,x2+dx32,t+dt,L);

  x1+=(dx11+2*dx21+dx31+dx41)/6;                       x2+=(dx12+2*dx22+dx32+dx42)/6;
  
  t+=dt;
}
