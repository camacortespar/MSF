#include <iostream>
#include <cmath>


double fun(double x, double t);
void pasoeuler(double & x, double & t, double dt);


int main(void)
{
  std::cout.precision(8);
  std::cout.setf(std::ios::scientific);
  
  double t=0;
  double x=0;
  double dt=0.1;
  
  for(t=0,x=1;t<10;){
    std::cout << t <<"\t" << x <<std::endl;
    pasoeuler(x, t, dt);    
  }


  
    return 0;
}


double fun(double x, double t)
{

  return -x/2.0;
}
void pasoeuler(double & x, double & t, double dt)
{
  double dx=0;
  dx=dt*fun(x,t);
  x+=dx;
  t+=dt;

}
