#include <iostream>
#include <cmath>


double fun(double t, double alpha, double x);
double simpson (double a, double b, int n, double alpha, double x);
double bessel (double a, double b, double n, double alpha, double x);

int main(void)
{
  std::cout.precision(8);
  std::cout.setf(std::ios::scientific);
  
  double a=0, b=M_PI, nmedios=2000;
  double alpha=0, x=0;
  
  for(x=0; x<=10; x+=0.01){
    std::cout << x<< "\t" << bessel(a,b,nmedios, alpha, x) << std::endl;
  }

    return 0;
}
//funcion para integrar
double fun(double t, double alpha, double x)
{
  return (std::cos(alpha*t-x*(std::sin(t))));
}
//integracion
double simpson (double a, double b, int n, double alpha, double x)
{
  n*=2;
  double h=(b-a)/n;
  double suma=0, t=0;
  
  for(int ii=0; ii<=n; ii++){
    t=a+h*ii;
    
    if(ii==0 || ii==n){
      suma+=fun(t,alpha,x);
    }
    else{ if(ii%2==0){suma+=2*fun(t,alpha,x);}
      else{suma+=4*fun(t,alpha,x);}
    }
  }
    return suma*h/3;
  }

double bessel (double a, double b, double n, double alpha, double x)
{
  return simpson(a, b, n, alpha, x)*1.0/M_PI;
}
