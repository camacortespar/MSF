#include <iostream>
#include <fstream>
#include <cmath>


double fun1(double x1, double x2, double t, double L);
double fun2(double x1, double x2, double t, double L);
void pasoRungeKutta4(double & x1, double & x2, double & t, double dt, double L);
template <typename F>
double bisection(double xl0, double xu0, double eps, F f);


struct Functor{
  double x0=0;
  double xf=0;
  double t=0;
  double dt=0.01;
 
  double operator() (double z){
   
    for(t=0.01, x0=0, xf=1; t<=1;){ //colocar las condiciones iniciales
      pasoRungeKutta4(x0, xf, t, dt, z);    
    }
    // std::cout<< xf << std::endl;
    return xf;
  };
};

int main(void)
{
  std::cout.precision(8);
  std::cout.setf(std::ios::scientific);

//x2 es R y x1 es la derivada de R  
  double t=0.01;
  double x1=0, x2=1; //condiciones iniciales
  double dt=0.01;
  double L=1;
  double eps=1e-10;
  
  double L1=0, L2=0, L3=0, L4=0;
  double x11=0,x21=0,x31=0,x41=0,x12=0,x22=0,x32=0,x42=0,t1=0,t2=0,t3=0,t4=0;

  Functor fun;
  fun.dt=dt;
  
  L1= bisection(2,4,eps,fun);
  L2= bisection(6,8,eps,fun);
  L3= bisection(9,11,eps,fun);
  L4= bisection(14,16,eps,fun);
  
  x11=x1, x21=x1, x31=x1, x41=x1;
  x12=x2, x22=x2, x32=x2, x42=x2;
  
  std::ofstream myfile("datos.dat");
  
  for(t1=t2=t3=t4=0.01, x1=0, x2=1; t4<=1;){

    pasoRungeKutta4(x11, x12, t1, dt, L1);
    pasoRungeKutta4(x21, x22, t2, dt, L2);
    pasoRungeKutta4(x31, x32, t3, dt, L3);
    pasoRungeKutta4(x41, x42, t4, dt, L4);
    
    //std::cout << t1 <<"\t" << x12 <<"\t" << x22 <<"\t" << x32 <<"\t" << x42 <<std::endl;
    myfile << t1 << "\t" << x12 <<"\t" << x22 <<"\t" << x32 <<"\t" << x42 <<std::endl; 
  }
  myfile.close();
  
  //plot de gnuplot
  //std::cout <<"plot \"datos.dat\" u 1:2 w l t \"L1\",\"\"u 1:3 w l t \"L2\",\"\"u 1:4 w l t \"L3\",\"\"u 1:5 w l t \"L4\"" << std::endl;
  //std::cout <<"pause 7" << std::endl;
    std::cout <<"#los lambda son: "<< L1 <<"\t" << L2 <<"\t" << L3 <<"\t" << L4 <<std::endl;
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

template <typename F>
double bisection(double xl0, double xu0, double eps, F f)
{
  int iter = 1;
  double mult=f(xl0)*f(xu0);
  double xl=xl0;
  double xu=xu0;
  
  while (mult < 0){
    double xr = (xl + xu)/2;
    //std::cout << iter << "\t"<<xl<<"\t" <<xu<<"\t" <<xr<<"\t"<< f(xr)<<"\n";
    if(std::fabs(f(xr))<=eps){
      return xr;
    }
    if(f(xr)*f(xu)<0){
      xl=xr;
    } else {
      xu=xr;
    }
    mult=f(xu)*f(xl);
    iter++;
  }
  return 0.5*(xl+xu);
}
