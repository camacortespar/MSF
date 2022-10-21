#include <iostream>
#include <cmath>


double fun(double x);
double simpson (double a, double b, int n);

int main(void)
{
  std::cout.precision(15);
  std::cout.setf(std::ios::scientific);
  
  double a=0, b=M_PI, nmedios=2400; 


  std::cout << "la integral es =" << simpson(a,b, nmedios) << std::endl;
  
  return 0;
}
//funcion para integrar
double fun(double x)
{
  return std::sin(x);
}
//integracion
double simpson (double a, double b, int n)
{
  n*=2;
  double h=(b-a)/n;
  double suma=0, x=0;
  
  for(int ii=0; ii<=n; ii++){
    x=a+h*ii;
    
    if(ii==0 || ii==n){
      suma+=fun(x);
    }
    else{ if(ii%2==0){suma+=2*fun(x);}
      else{suma+=4*fun(x);}
    }
  }
    return suma*h/3;
  }



