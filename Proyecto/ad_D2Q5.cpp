/*
  Proyecto Final/Grupo 02:
  Lattice Boltzmann para la ecuación de Advección - Difusión en coordenadas cartesianas.
  Prueba con un pulso gaussiano como campo impuesto.
*/

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

//-----Constantes Globales-----//
const int Lx = 180;                 //tamaño de la simulacion
const int Ly = 400;

const int Q = 5;                    //numero de direcciones
const double W0 = 1.0/3.0;          //cte que define los pesos

const double C = 0.5;               //velocidad de onda
const double C2 = C*C;

const double tau = 0.7;             //valores en la funcion de colision - D = 0.05
const double Utau = 1.0/tau;
const double UmUtau = 1-Utau;

//-----Clase LatticeBoltzmann------//
class LatticeBoltzmann{
private:
  double w[Q];        //pesos por direccion
  int Vx[Q], Vy[Q];   //vectores de velocidad
  double *f, *fnew;   //funciones de distribucion
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int n(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;};   //indice de las celdas
  double rho(int ix, int iy, bool UseNew);
  double Jx(int ix, int iy, bool UseNew);
  double Jy(int ix, int iy, bool UseNew);
  double feq(double rho0, double Ux0, double Uy0, int i);
  void Start(double rho0, double Ux0, double Uy0);
  void Collision(void);
  void ImposeFields(double Ux0, double Uy0, int t);
  void Advection(void);
  void Print(const char * NameFile);
};
//Constructor
LatticeBoltzmann::LatticeBoltzmann(void){
  //Pesos
  w[0] = W0; w[1] = w[2] = w[3] = w[4] = (1.0-W0)/4.0;
  //Vectores de velocidad
  Vx[0] = 0;  Vx[1] = 1;  Vx[2] = 0;  Vx[3] = -1; Vx[4] = 0;
  Vy[0] = 0;  Vy[1] = 0;  Vy[2] = 1;  Vy[3] = 0;  Vy[4] = -1;
  //Crea arrays dinamicos
  int ArraySize = Lx*Ly*Q;
  f = new double [ArraySize];  fnew = new double [ArraySize];
}
//Destructor
LatticeBoltzmann::~LatticeBoltzmann(void){
  delete[] f;  delete[] fnew;
}
//Rho
double LatticeBoltzmann::rho(int ix, int iy, bool UseNew){
  double sum; int i, n0;
  for(sum=0, i=0; i<Q; i++){
    n0 = n(ix,iy,i);
    if(UseNew) sum += fnew[n0]; else sum += f[n0];
  }
  return sum;
} 
//Componentes del vector J 
double LatticeBoltzmann::Jx(int ix, int iy, bool UseNew){
  double sum; int i, n0;
  for(sum=0, i=0; i<Q; i++){
    n0 = n(ix,iy,i);
    if(UseNew) sum += Vx[i]*fnew[n0]; else sum += Vx[i]*f[n0];
  }
  return sum;
}
double LatticeBoltzmann::Jy(int ix, int iy, bool UseNew){
  double sum; int i, n0;
  for(sum=0, i=0; i<Q; i++){
    n0 = n(ix,iy,i);
    if(UseNew) sum += Vy[i]*fnew[n0]; else sum += Vy[i]*f[n0];
  }
  return sum;
}
//Funcion equilibrio
double  LatticeBoltzmann::feq(double rho0, double Ux0, double Uy0, int i){
  double UdotVi = Ux0*Vx[i]+Uy0*Vy[i], U2 = Ux0*Ux0+Uy0*Uy0;
  return rho0*w[i]*(1+(UdotVi/C2)+(pow(UdotVi,2)/(2*pow(C2,2)))-(U2/(2*C2)));  
} 
//Start
void LatticeBoltzmann::Start(double rho0, double Ux0, double Uy0){
  int ix, iy, i, n0;
  for(ix=0; ix<Lx; ix++)    //para cada celda
    for(iy=0; iy<Ly; iy++)
      for(i=0; i<Q; i++){   //en cada direccion
        n0 = n(ix,iy,i);
        f[n0] = feq(rho0,Ux0,Uy0,i);
      }
}  
//Colision
void LatticeBoltzmann::Collision(void){
  int ix, iy, i, n0; double rho0, Ux0, Uy0;
  for(ix=0; ix<Lx; ix++)      //para cada celda
    for(iy=0; iy<Ly; iy++){
      //Calcule los campos macroscopicos en la celda
      rho0 = rho(ix,iy,false); Ux0 = Jx(ix,iy,false)/rho0; Uy0 = Jy(ix,iy,false)/rho0;
      for(i=0; i<Q; i++){     //para cada vector de velocidad
        n0 = n(ix,iy,i);
        fnew[n0] = UmUtau*f[n0]+Utau*feq(rho0,Ux0,Uy0,i);
      }
    }  
}
//Imponer campos
void LatticeBoltzmann::ImposeFields(double Ux0, double Uy0, int t){
  int i, ix, iy, ix0, iy0, n0;
  double rho_gauss, x0_bar, y0_bar;
  //Constantes Gaussian Hill
  double A = 1.0, D = C2*(tau-0.5), Sigma = 8.0, Sigma2 = Sigma*Sigma;
  double Coef = A/(1+(2*D*t/Sigma2));
  //Fuente gaussiana en el medio
  ix0 = Lx/2; iy0 = Ly/2;
  for(ix=0; ix<Lx; ix++){      //para cada celda
    for(iy=0; iy<Ly; iy++){
      x0_bar = ix0+Ux0*t; y0_bar = iy0+Uy0*t;
      rho_gauss = Coef*exp(-(pow(ix-x0_bar,2)+pow(iy-y0_bar,2))/(2*(Sigma2+2*D*t)));   
      //K=0.05 es usada para evitar NaN.
      for(i=0; i<Q; i++){
        n0 = n(ix,iy,i); 
        fnew[n0] = feq(rho_gauss,Ux0,Uy0,i);
      }
    }
  }
}
//Adveccion
void LatticeBoltzmann::Advection(void){
  int ix, iy, i, ixnext, iynext, n0, n0next;
  for(ix=0; ix<Lx; ix++)      //para cada celda
    for(iy=0; iy<Ly; iy++)
      for(i=0; i<Q; i++){     //en cada direccion
        ixnext = (ix+Vx[i]+Lx)%Lx; iynext = (iy+Vy[i]+Ly)%Ly;
        n0 = n(ix,iy,i); n0next = n(ixnext,iynext,i);
        f[n0next] = fnew[n0];     //fronteras periodicas
      }
}
//Print
void LatticeBoltzmann::Print(const char * NameFile){
  ofstream MyFile(NameFile); double rho0; int ix, iy;
  for(ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      rho0=rho(ix,iy,true);
      MyFile<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}

//-----Programa Principal-----//

int main(void){
  LatticeBoltzmann Air;
  int t, tmax = 10;
  double rho0 = 1.0, Ux0 = 1.0, Uy0 = -1.5; 

  Air.Start(rho0, Ux0, Uy0);
  //Evolucione
  for(t=0; t<tmax; t++){
    Air.Collision();
    Air.ImposeFields(Ux0, Uy0, t);
    Air.Advection();
  }
  //Imprima
  Air.Print("ad.dat");

  return 0;
}
