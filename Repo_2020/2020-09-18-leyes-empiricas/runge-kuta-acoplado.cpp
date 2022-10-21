#include <iostream>
#include <cmath>


double fun1(double x1, double x2, double t);
double fun2(double x1, double x2, double t);
void pasoRungeKutta4(double & x1, double & x2, double & t, double dt);

const double Omega=1;

int main(void)
{
  std::cout.precision(8);
  std::cout.setf(std::ios::scientific);
  
  double t=0;
  double x1=0, x2=0;
  double dt=0.1;
  
  for(t=0,x1=1, x2=0;t<10;){
    std::cout << t <<"\t" << x1 <<"\t"<< x2 <<std::endl;
    pasoRungeKutta4(x1, x2, t, dt);    
  }
  
    return 0;
}

double fun1(double x1, double x2, double t) 
{
  return -std::pow(Omega,2)*x2;
}
double fun2(double x1, double x2, double t)
{
  return x1;
}
void pasoRungeKutta4(double & x1, double & x2, double & t, double dt)
{
  double dx11, dx21, dx31, dx41;                      double dx12, dx22, dx32, dx42;
  dx11=dt*fun1(x1,x2,t);                             dx12=dt*fun2(x1,x2,t);
  dx21=dt*fun1(x1+0.5*dx11,x2+0.5*dx12,t+0.5*dt);    dx22=dt*fun2(x1+0.5*dx11,x2+0.5*dx12,t+0.5*dt);
  dx31=dt*fun1(x1+0.5*dx21,x2+0.5*dx22,t+0.5*dt);    dx32=dt*fun2(x1+0.5*dx21,x2+0.5*dx22,t+0.5*dt);
  dx41=dt*fun1(x1+dx31,x2+dx32,t+dt);                dx42=dt*fun2(x1+dx31,x2+dx32,t+dt);

  x1+=(dx11+2*dx21+dx31+dx41)/6;                     x2+=(dx12+2*dx22+dx32+dx42)/6;
  
  t+=dt; 
}
