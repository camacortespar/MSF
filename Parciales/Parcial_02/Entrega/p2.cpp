/*
  Parcial 02/Punto 02: Difusion 2D - Flujo de Poiseuille
  Camilo Andrés Cortés Parra
  camacortespar@unal.edu.co
*/

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

//-----Constantes Globales-----//
const int Lx = 256;                  //tamaño de la simulacion
const int Ly = 64;

const int Q = 9;                     //numero de direcciones

const double tau = 0.55;             //valores en la funcion de colision
const double Utau = 1.0/tau;
const double UmUtau = 1-Utau;

//-----Clase LatticeBoltzmann------//
class LatticeBoltzmann{
private:
  double w[Q];        //pesos por direccion
  int Vx[Q], Vy[Q];   //vectores de velocidad
  double *f, *fnew;   //funciones de distribucion - * es un apuntador que lo lleva a la direccion de memoria
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
  //void ImposeFields(void);
  void Advection(void);
  double Varianza(void);
  void Print(const char * NameFile, double Ufan);
};
//Constructor
LatticeBoltzmann::LatticeBoltzmann(void){   
   //Pesos
  w[0]=4.0/9;  w[1]=w[2]=w[3]=w[4]=1.0/9;  w[5]=w[6]=w[7]=w[8]=1.0/36;
  //Vectores de velocidad
  Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;
  Vx[5]=1;  Vx[6]=-1; Vx[7]=-1; Vx[8]=1;
  Vy[5]=1;  Vy[6]=1;  Vy[7]=-1; Vy[8]=-1;
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
  return 0.0;
}
double LatticeBoltzmann::Jy(int ix, int iy, bool UseNew){
  return 0.0;
} 
//Funcion equilibrio
double  LatticeBoltzmann::feq(double rho0, double Ux0, double Uy0, int i){
  double UdotVi = Ux0*Vx[i]+Uy0*Vy[i], U2 = Ux0*Ux0+Uy0*Uy0;
  return w[i]*rho0*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
  } 
//Start
void LatticeBoltzmann::Start(double rho0, double Ux0, double Uy0){
  int ix, iy, i, n0;
  double sigma = Ly/9, mu = Lx/2, coef = 1/(sigma*sqrt(2*M_PI));  //informacion gaussiana - modificacion
  double rho_gauss;
  for(ix=0; ix<Lx; ix++)    //para cada celda
    for(iy=0; iy<Ly; iy++)
      for(i=0; i<Q; i++){   //en cada direccion
        n0 = n(ix,iy,i);
        rho_gauss = coef*exp(-0.5*pow((ix-mu)/sigma,2));
        f[n0] = feq(rho_gauss,0,0,i);
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
//Adveccion
void LatticeBoltzmann::Advection(void){
  int ix, iy, i, ixnext, iynext, n0, n0next;
  for(ix=0; ix<Lx; ix++)      //para cada celda
    for(iy=0; iy<Ly; iy++)
      for(i=0; i<Q; i++){     //en cada direccion
        ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
        n0=n(ix,iy,i); n0next=n(ixnext,iynext,i);
        f[n0next]=fnew[n0];   //fronteras periodicas
    }
}
//Varianza
double LatticeBoltzmann::Varianza(void)
{
  int ix, iy;
  double N, xprom, var, sum_n, sum_x, sum_v;
  //Calcule N
  for(ix=0; ix<Lx; ix++)
    for(iy=0; iy<Ly; iy++)
      sum_n += rho(ix,iy,true);
  N = sum_n;
  //Calcule xprom
  for(ix=0; ix<Lx; ix++)
    for(iy=0; iy<Ly; iy++)
      sum_x += rho(ix,iy,true)*ix;
  xprom = sum_x/N;
  //Calcule varianza
  for(ix=0; ix<Lx; ix++)
    for(iy=0; iy<Ly; iy++)
      sum_v += rho(ix,iy,true)*pow(ix-xprom,2);
  var = sum_v/N;

  return var;
}
//Print
void LatticeBoltzmann::Print(const char * NameFile, double Ufan){
  ofstream MyFile(NameFile); double rho0, Ux0, Uy0; int ix, iy;
  for(ix=0; ix<Lx; ix+=4){
    for(iy=0; iy<Ly; iy+=4){
      rho0 = rho(ix,iy,true);
      MyFile<<ix<<"\t"<<iy<<"\t"<<rho0<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}

//-----Programa Principal-----//

int main(void){
  LatticeBoltzmann Air;
  int t, tmax = 1000;
  double rho0 = 1.0, Ufan0 = 0.1; 

  Air.Start(rho0, Ufan0, 0);
  //Evolucione
  for(t=0; t<tmax; t++){
    Air.Collision();
    Air.Advection();
    cout<<t<<"\t"<<Air.Varianza()<<endl;
  }
  //Imprima
  //Air.Print("p2.dat", Ufan0);

  return 0;
}