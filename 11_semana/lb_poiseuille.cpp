/*
  Semana 11:
  Simulacion de fluidos
*/

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

//-----Constantes Globales-----//
const int Lx = 256;                 //tamaño de la simulacion
const int Ly = 64;

const int Q = 9;                    //numero de direcciones

const double tau = 0.55;             //valores en la funcion de colision
const double Utau = 1.0/tau;
const double UmUtau = 1-Utau;

//-----Clase LatticeBoltzmann------
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
  double feq(double rho0,double Ux0, double Uy0, int i);
  void Start(double rho0, double Ux0, double Uy0);
  void Collision(void);
  void ImposeFields(double Ufan);
  void Advection(void);
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
  return rho0*w[i]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
  } 
//Start
void LatticeBoltzmann::Start(double rho0, double Ux0, double Uy0){
  int ix,iy,i,n0;
  for(ix=0; ix<Lx; ix++)    //para cada celda
    for(iy=0; iy<Ly; iy++)
      for(i=0; i<Q; i++){   //en cada direccion
    n0 = n(ix,iy,i);
    f[n0] = feq(rho0,Ux0,Uy0,i);
      }
}  
//Colision
void LatticeBoltzmann::Collision(void){
  int ix,iy,i,n0; double rho0,Ux0,Uy0;
  for(ix=0;ix<Lx;ix++)      //para cada celda
    for(iy=0;iy<Ly;iy++){
      //Calcule los campos macroscopicos en la celda
      rho0=rho(ix,iy,false); Ux0=Jx(ix,iy,false)/rho0; Uy0=Jy(ix,iy,false)/rho0;
      for(i=0;i<Q;i++){     //para cada vector de velocidad
        n0 = n(ix,iy,i);
        fnew[n0] = UmUtau*f[n0]+Utau*feq(rho0,Ux0,Uy0,i);
      }
    }  
}
//Imponer campos, en este caso se impone geometría de ventilador + obstaculo
void LatticeBoltzmann::ImposeFields(double Ufan){
  int i,ix,iy,n0;
  double rho0; int ixc = Lx/8, iyc = Ly/2, R = Ly/5; double R2 = R*R;   //ctes que definen obstaculo circular
  for(ix=0;ix<Lx;ix++)      //para cada celda
    for(iy=0;iy<Ly;iy++){
      rho0 = rho(ix, iy, false);
      //Ventilador
      if(ix==0)
        for(i=0;i<Q;i++){n0=n(ix,iy,i); fnew[n0]=feq(rho0,Ufan,0,i);}
      //Obstaculo
      else if(pow((ix-ixc),2)+pow((iy-iyc),2)<=R2)
        for(i=0;i<Q;i++){n0=n(ix,iy,i); fnew[n0]=feq(rho0,0,0,i);}
      //Hay que agregar un punto de perturbacion para romper la isotropia
      else if(ix==ixc && iy==iyc+R+1)
        for(i=0;i<Q;i++){n0=n(ix,iy,i); fnew[n0]=feq(rho0,0,0,i);}
    }
}

//Adveccion
void LatticeBoltzmann::Advection(void){
  int ix,iy,i,ixnext,iynext,n0,n0next;
  for(ix=0;ix<Lx;ix++)      //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){     //en cada direccion
    ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
    n0=n(ix,iy,i); n0next=n(ixnext,iynext,i);
    f[n0next]=fnew[n0];     //fronteras periodicas
      }
}
//Print
void LatticeBoltzmann::Print(const char * NameFile, double Ufan){
  ofstream MyFile(NameFile); double rho0,Ux0,Uy0; int ix,iy;
  for(ix=0;ix<Lx;ix+=4){
    for(iy=0;iy<Ly;iy+=4){
      rho0=rho(ix,iy,true); Ux0=Jx(ix,iy,true)/rho0; Uy0=Jy(ix,iy,true)/rho0;
      MyFile<<ix<<" "<<iy<<" "<<Ux0/Ufan*4<<" "<<Uy0/Ufan*4<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}

//-----Programa Principal-----

int main(void){
  LatticeBoltzmann Air;
  int t, tmax = 10000;
  double rho0 = 1.0, Ufan0 = 0.1; 

  Air.Start(rho0, Ufan0, 0);
  //Evolucione
  for(t=0; t<tmax; t++){
    Air.Collision();
    Air.ImposeFields(Ufan0);
    Air.Advection();
  }
  //Imprima
  Air.Print("Air.dat", Ufan0);

  return 0;
}
