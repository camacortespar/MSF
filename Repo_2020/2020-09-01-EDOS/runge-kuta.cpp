#include <iostream>
#include <cmath>


double fun(double x, double t);
void pasoRungeKutta4(double & x, double & t, double dt);


int main(void)
{
  std::cout.precision(8);
  std::cout.setf(std::ios::scientific);
  
  double t=0;
  double x=0;
  double dt=0.1;
  
  for(t=0,x=1;t<10;){
    std::cout << t <<"\t" << x <<std::endl;
    pasoRungeKutta4(x, t, dt);    
  }
  
    return 0;
}

double fun(double x, double t)
{
  return x/2.0;
}
void pasoRungeKutta4(double & x, double & t, double dt)
{
  double dx1, dx2, dx3, dx4;
  dx1=dt*fun(x,t);
  dx2=dt*fun(x+0.5*dx1, t+0.5*dt);
  dx3=dt*fun(x+0.5*dx2, t+0.5*dt);
  dx4=dt*fun(x+dx3, t+dt);

  x+=(dx1+2*dx2+dx3+dx4)/6; t+=dt;
}
