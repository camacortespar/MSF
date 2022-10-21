//Ejercicio 1.c

#include <iostream>
#include <cmath>

using namespace std;

//ctes de contagio y tiempo de recuperacion
const double Gamma = 0.08;  //fijo
const int N = 500; //numero de betas
const double B_min = 0.081;  //satisface la condicion R0*s(0)>1 
const double B_max = 0.5;
const double t_final = 7000; //tiempo de epidemia muy grande

//declaracion de funciones
double f1(double s, double i, double t, double Beta);  //ds/dt
double f2(double s, double i, double t, double Beta);  //di/dt
void RK4(double & s, double & i, double & t, double & dt, double Beta);  //runge-kutta de orden 4 para las EDO's acopladas

//implementacion de funciones
double f1(double s, double i, double t, double Beta)
{
  return -Beta*s*i;
}

double f2(double s, double i, double t, double Beta)
{
  return Beta*s*i - Gamma*i;
}

void RK4(double & s, double & i, double & t, double & dt, double Beta)
{
  //incrementos k para s y l para i
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;

  k1 = dt*f1(s,i,t,Beta);
  l1 = dt*f2(s,i,t,Beta);
  
  k2 = dt*f1(s + k1/2, i + l1/2, t + dt/2,Beta);
  l2 = dt*f2(s + k1/2, i + l1/2, t + dt/2,Beta);

  k3 = dt*f1(s + k2/2, i + l2/2, t + dt/2,Beta);
  l3 = dt*f2(s + k2/2, i + l2/2, t + dt/2,Beta);

  k4 = dt*f1(s + k3, i + l3, t + dt,Beta);
  l4 = dt*f2(s + k3, i + l3, t + dt,Beta);

  //actualizacion de variables
  s += (k1 + 2*k2 + 2*k3 + k4)/6;
  i += (l1 + 2*l2 + 2*l3 + l4)/6;
  t += dt;
}

int main(void){
  
  cout.precision(8);
  //condiciones iniciales
  double s = 0.999;
  double i = 0.001;
  double t = 0;
  double dt = 0.1;
  double Beta, dBeta;  //variacion de beta

  //comportamiento de s_inf vs R0
  double aux, R0;
  dBeta = (B_max-B_min)/N;
  for (int j = 0; j < N; ++j){
    Beta = B_min + j*dBeta;
    R0 = Beta/Gamma;
    while (t < t_final){
      aux = s;
      RK4(s, i, t, dt, Beta);
      if (aux == s){  //se detiene si ya se llego s_inf
	break;
      }
    }
    cout << R0 << " " << s << " " << Beta <<endl;
    //reiniciamos valores
    t = 0;
    s = 0.999;
  }
    return 0;
}
  
