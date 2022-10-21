#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "Random64.h"
using namespace std;

const int Lx=200;
const int Ly=200;

const int Q=5;
const double W0=1.0/3;

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

//--------------------- Clase LatticeBoltzmann ------------
class LatticeBoltzmann{
private:
  vector <double> w;
  vector<vector<int>> V;
  vector<vector<vector<double>>> f;
  vector<vector<vector<double>>> fnew;
  
public:
  LatticeBoltzmann(void);
  double Ccelda(int ix,int iy);
  double rho(int ix,int iy,bool UseNew);
  double Jx(int ix,int iy,bool UseNew);
  double Jy(int ix,int iy,bool UseNew);
  double feq(double rho0,double Jx0,double Jy0,int i, int ix, int iy);
  void Colisione(void);
  void ImponerCampos(int t);
  void Adveccione(void);
  void Inicie(double rho0,double Jx0,double Jy0);
  void Imprimase(const char * NombreArchivo);
};
double LatticeBoltzmann::Ccelda(int ix,int iy)
{
  double C=0.5;
    
  return C;  
}
LatticeBoltzmann::LatticeBoltzmann(void){
  int q=Q, llx=Lx, lly=Ly;
  w.resize(Q);
  V.resize(2,vector<int>(Q)); //V[0][i]=Vi_x ,  V[1][i]=Vi_y
  f.resize(Lx, vector<vector<double>>(Ly, vector<double>(Q)));;
  fnew.resize(Lx, vector<vector<double>>(Ly, vector<double>(Q))); // f[ix][iy][i]


  
  //cargar los pesos
  w[0]=W0; w[1]=w[2]=w[3]=w[4]=(1.0-W0)/4;
  //cargar los vectores velocidad
  V[0][0]=0;  V[0][1]=1;  V[0][2]=0;  V[0][3]=-1; V[0][4]=0;
  V[1][0]=0;  V[1][1]=0;  V[1][2]=1;  V[1][3]=0;  V[1][4]=-1;
}
double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  double suma; int i;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]; else suma+=f[ix][iy][i];
  return suma;
}  
double LatticeBoltzmann::Jx(int ix,int iy,bool UseNew){
  double suma; int i;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=V[0][i]*fnew[ix][iy][i]; else suma+=V[0][i]*f[ix][iy][i];
  return suma;
}  
double LatticeBoltzmann::Jy(int ix,int iy,bool UseNew){
  double suma; int i;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=V[1][i]*fnew[ix][iy][i]; else suma+=V[1][i]*f[ix][iy][i];
  return suma;
}  
double  LatticeBoltzmann::feq(double rho0,double Jx0,double Jy0,int i, int ix, int iy){
  if(i>0)
    return 3*w[i]*(std::pow(Ccelda(ix,iy),2)*rho0+V[0][i]*Jx0+V[1][i]*Jy0);
  else
    return rho0*(1-3*std::pow(Ccelda(ix,iy),2)*(1-W0));
}  
void LatticeBoltzmann::Colisione(void){
  int ix,iy,i,ixc=50,iyc=100; double rho0,Jx0,Jy0,R2=std::pow(100,2);
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++){
      //calcular los campos macroscópicos en la celda    
        //espejo
      if((ix-ixc)*(ix-ixc)+(iy-iyc)*(iy-iyc)>R2 && ix>ixc){
	Jx0=-1*Jx(ix,iy,false); Jy0=-1*Jy(ix,iy,false);}
      else{
      rho0=rho(ix,iy,false); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
      for(i=0;i<Q;i++) //en cada dirección
	fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(rho0,Jx0,Jy0,i,ix,iy);}
    }  
}
void LatticeBoltzmann::ImponerCampos(int t){
  int i,ix,iy;
  double lambda,omega,rho0,Jx0,Jy0; lambda=10; omega=2*M_PI/lambda*Ccelda(ix,iy);
  ix=0;
  rho0=10*sin(omega*t); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
  for(iy=0;iy<Ly;iy++)
    for(i=0;i<Q;i++)
      fnew[ix][iy][i]=feq(rho0,Jx0,Jy0,i,ix,iy);

}
void LatticeBoltzmann::Adveccione(void){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++) //en cada dirección
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i]; //fronteras periódicas
}
void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++) //en cada dirección
	f[ix][iy][i]=feq(rho0,Jx0,Jy0,i,ix,iy);
}  
void LatticeBoltzmann::Imprimase(const char * NombreArchivo){
  ofstream MiArchivo(NombreArchivo); double rho0;
  for(int ix=0;ix<200;ix++){
    for(int iy=0;iy<200;iy++){
      rho0=rho(ix,iy,true);
      MiArchivo<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    MiArchivo<<endl;
  }
  MiArchivo.close();
}

//------------------- Funciones Globales ------------

int main(void){
  LatticeBoltzmann Ondas;
  int t,tmax=100000;
  double rho0=0,Jx0=0,Jy0=0;
  
  //Inicie
  Ondas.Inicie(rho0,Jx0,Jy0);
  //Corra
  for(t=0;t<tmax;t++){
    Ondas.Colisione();
    Ondas.ImponerCampos(t);
    Ondas.Adveccione();
  }
  Ondas.Imprimase("Ondas.dat");
  
  return 0;
}  
