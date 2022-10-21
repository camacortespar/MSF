#include <iostream>
#include <cmath>


using fptr = double(double);  
double bisection(double xa, double xb, fptr f);
double fun(double x);


const double ERR= 1e-7;


int main(void)
{
  std::cout.precision(15);
  std::cout.setf(std::ios::scientific);
  
  double a=2, b=4; // b debe ser mayor que a


  std::cout << "la raiz está en x=" << bisection(a,b, fun) << std::endl;
  
  return 0;
}



double fun(double x)
{
  return std::sin(x)/x;
}


double bisection(double xa, double xb, fptr f)
{
  double fa =fun(xa);
  double  m=0;
  double fm =0;
while(xb-xa > ERR){
    //calcular m y f(m)
     m=(xa+xb)/2;
    fm = fun(m);
    //si fa y fm son del mismo signo, mover a hacia m
    if(fa*fm>0){
      xa=m;
      fa=fm;
    }
    //si no, mover b hacia m
    else{
      xb=m;
    }
 }
  return (xa+xb)/2;
}
