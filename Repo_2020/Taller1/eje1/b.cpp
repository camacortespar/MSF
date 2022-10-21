#include <iostream>
#include <fstream>
#include <cmath>

double fun1(double x1, double x2, double t);
double fun2(double x1, double x2, double t);
double fun3(double x1, double x2, double t);
void pasoRungeKutta4(double & x1, double & x2,double & x3, double & t, double dt);

//-------------------------------funciones de la solucion analitica----------------

double zita (double t,double S0, double R0);
double Infectados(double t,double S0, double R0);
double Suceptibles(double t,double S0, double R0);
double Recuperados(double t,double S0, double R0);

const double beta=0.35;
const double gama=0.08;
const double N=1;
const double ro=(N*gama)/beta;


int main(void)
{
  std::cout.precision(8);
  std::cout.setf(std::ios::scientific);
 
  double t=1;
  double s=0.999, i=0.001, r=0;
  double dt=1;
  double R0=r, I0=i, S0=s;
  double n=500;
  
  //std::ofstream myfile("datos.dat");
  
    for(t, s, i, r; t<= 150 ; ){

      // myfile << t <<"\t"<< x2 <<std::endl;
      std::cout << t <<"\t" << s <<"\t" << i <<"\t" << r <<"\t" << Suceptibles(t,S0,R0)
		<<"\t" << Infectados(t,S0,R0) <<"\t" << Recuperados(t,S0,R0) << std::endl;

      pasoRungeKutta4( s, i, r, t, dt);    
  }

    //myfile.close();
  
  return 0;
}

double fun1(double x1, double x2, double t) 
{
  return -beta*x2*x1;
}
double fun2(double x1, double x2, double t)
{
  return beta*x1*x2-gama*x2;
}
double fun3(double x1, double x2, double t)
{
  return gama*x2;
}
void pasoRungeKutta4(double & x1, double & x2,double & x3, double & t, double dt)
{
  double dx11, dx21, dx31, dx41;                      double dx12, dx22, dx32, dx42;                            double dx13, dx23, dx33, dx43;
  dx11=dt*fun1(x1,x2,t);                             dx12=dt*fun2(x1,x2,t);                                     dx13=dt*fun3(x1,x2,t);
  dx21=dt*fun1(x1+0.5*dx11,x2+0.5*dx12,t+0.5*dt);    dx22=dt*fun2(x1+0.5*dx11,x2+0.5*dx12,t+0.5*dt);            dx23=dt*fun3(x1+0.5*dx11,x2+0.5*dx12,t+0.5*dt);
  dx31=dt*fun1(x1+0.5*dx21,x2+0.5*dx22,t+0.5*dt);    dx32=dt*fun2(x1+0.5*dx21,x2+0.5*dx22,t+0.5*dt);            dx33=dt*fun3(x1+0.5*dx21,x2+0.5*dx22,t+0.5*dt);
  dx41=dt*fun1(x1+dx31,x2+dx32,t+dt);                dx42=dt*fun2(x1+dx31,x2+dx32,t+dt);                   dx43=dt*fun3(x1+dx31,x2+dx32,t+dt);
  
  x1+=(dx11+2*dx21+dx31+dx41)/6;                     x2+=(dx12+2*dx22+dx32+dx42)/6;                        x3+=(dx13+2*dx23+dx33+dx43)/6;
  t+=dt;
}



double zita (double t, double S0, double R0)
{
  double dt=0;
  double sum=0;
  double term=0;
  double N=t/0.01
    Infectado[N];
  Infectado[0]=N-S0-R0;


  for(int i=0, dt=0; dt<=t; dt+=0.01, i++){

    term=term + dt*(beta/N)*Infectado[i];
    sum+=term;
  }
  return sum;
}

double Suceptibles(double t, double S0, double R0)
{
  return S0*std::exp(-zita(t,S0,R0));
}

double Recuperados(double t, double S0, double R0)
{
  return R0+ro*zita(t,S0,R0);
}
double Infectados(double t, double S0, double R0)
{
  return (N - Suceptibles(t,S0,R0)-Recuperados(t,S0,R0));
}
