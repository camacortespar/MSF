//Un pulso gaussiano plano en Lattice-Boltzmann para Electrodinamica (Cambio de medio)
#include<iostream>
#include<cmath>
#include "Vector.h"
using namespace std;

//------------------------CONSTANTES-------------------------------
const int Lx = 1;   //
const int Ly = 1;   //
const int Lz = 200; //
//-------------------
const double Tau = 0.5;
const double UTau = 1/Tau;
const double UmUTau=1-1/Tau;
//-------------------
const double Epsilon0=1, Mu0=2;
const double Sigma1=0.0;
const double C=1.0/sqrt(2.0);

//------------------FUNIONES PARA EL MEDIO------------------------------
double mur(int ix,int iy,int iz){
  return 1.0;
}
double epsilonr(int ix,int iy,int iz){
  return 1.75+0.75*tanh((double)(iz-(Lz/2)));
}
double sigma(int ix,int iy,int iz){
  return 0.0;
}
double prefactor(double Epsilonr,double Sigma){
  return Sigma/(1+(Sigma*Mu0)/(4*Epsilonr));
}
//----------------------------- CLASES -----------------------------
class LatticeBoltzmann{
private:
  int V[3][3][4], V0[3]; /*V[xyz][p][i]*/  vector3D v[3][4],v0; //v[p][i]
  vector3D e[3][4][2], e0; //e[p][i][j]
  vector3D b[3][4][2], b0; //b[p][i][j]
  double f[Lx][Ly][Lz][2][3][4][2],fnew[Lx][Ly][Lz][2][3][4][2];//f[ix][iy][iz][r][p][i][j]
  double f0[Lx][Ly][Lz],f0new[Lx][Ly][Lz];//f0[ix][iy][iz] (r=0)
public:
  LatticeBoltzmann(void);
  //Campos de sumas directas
  double rhoc(int ix,int iy,int iz,bool UseNew);
  vector3D D(int ix,int iy,int iz,bool UseNew);
  vector3D B(int ix,int iy,int iz,bool UseNew);
  //Campos deducidos por constantes
  vector3D E(vector3D & D0,double Epsilonr);
  vector3D H(vector3D & B0,double Mur);
  //Campos primados
  vector3D Jprima(vector3D & E0,double Prefactor);
  vector3D Eprima(vector3D & E0,vector3D & Jprima0,double Epsilonr);
  //Funciones de equilibrio
  double feq(vector3D & Jprima0,vector3D & Eprima0,vector3D & B0,
	     double Epsilonr,double Mur,
	     int r,int p,int i,int j);
  double feq0(double rhoc0);
  //Funciones de simulación
  void Inicie(void);
  void Colisione(void);
  void Adveccione(void);
  void Muestre(void);
};

LatticeBoltzmann::LatticeBoltzmann(void){
  int ix,iy,iz,alpha,r,p,i,j;
  //Los vectores Velocidad V[alpha][p][i]=V^p_i_(alpha) en componentes
  V0[0]=V0[1]=V0[2]=0;

  V[0][0][0]=V[0][1][0]=V[1][2][0]=1;
  V[1][0][0]=V[2][1][0]=V[2][2][0]=1;
  V[2][0][0]=V[1][1][0]=V[0][2][0]=0;

  V[0][0][1]=V[0][1][1]=V[1][2][1]=-1;
  V[1][0][1]=V[2][1][1]=V[2][2][1]=1;
  V[2][0][1]=V[1][1][1]=V[0][2][1]=0;

  V[0][0][2]=V[0][1][2]=V[1][2][2]=-1;
  V[1][0][2]=V[2][1][2]=V[2][2][2]=-1;
  V[2][0][2]=V[1][1][2]=V[0][2][2]=0;

  V[0][0][3]=V[0][1][3]=V[1][2][3]=1;
  V[1][0][3]=V[2][1][3]=V[2][2][3]=-1;
  V[2][0][3]=V[1][1][3]=V[0][2][3]=0;
  //Los vectores Velocidad V[p][i]=V^p_i como vectores
  v0.cargue(V0[0],V0[1],V0[2]);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++){
      v[p][i].cargue(V[0][p][i],V[1][p][i],V[2][p][i]);
  }
  //Los vectores Electricos e[p][i][j]=e^p_{ij}
  e0.cargue(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++){
      e[p][i][0]=v[p][(i+1)%4]*0.5;
      e[p][i][1]=v[p][(i+3)%4]*0.5;
  }

  //Los vectores Magneticos b[p][i][j]=b^p_{ij}=v^p_i x e^p_{ij}
  b0.cargue(0,0,0);  
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
	b[p][i][j]=(v[p][i]^e[p][i][j]);
}
//-------------------CAMPOS-------------
//----------Campos de sumas directas
double LatticeBoltzmann::rhoc(int ix,int iy,int iz,bool UseNew){
  int p,i,j; double suma=f0[ix][iy][iz];
  for(p=0;p<2;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
	if(UseNew) 
	  suma+=fnew[ix][iy][iz][0][p][i][j]; 
	else 
	  suma+=f[ix][iy][iz][0][p][i][j];
  return suma;
}
vector3D LatticeBoltzmann::D(int ix,int iy,int iz,bool UseNew){
  int p,i,j; vector3D suma; suma.cargue(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
	if(UseNew) 
	    suma+=e[p][i][j]*fnew[ix][iy][iz][0][p][i][j];
	else
	    suma+=e[p][i][j]*f[ix][iy][iz][0][p][i][j];
  return suma;
}
vector3D LatticeBoltzmann::B(int ix,int iy,int iz,bool UseNew){
  int p,i,j; vector3D suma; suma.cargue(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
	if(UseNew) 
	  suma+=b[p][i][j]*fnew[ix][iy][iz][1][p][i][j];
	else
	  suma+=b[p][i][j]*f[ix][iy][iz][1][p][i][j];
  return suma;
}
//Campos deducidos por constantes
vector3D LatticeBoltzmann::E(vector3D & D0,double Epsilonr){
  return D0*(1.0/Epsilonr);
}
vector3D LatticeBoltzmann::H(vector3D & B0,double Mur){
  return B0*(1.0/Mur);
}
//Campos primados
vector3D Jprima(vector3D & E0,double Prefactor){
  return E0*Prefactor;
}
vector3D Eprima(vector3D & E0,vector3D & Jprima0,double Epsilonr){
  return E0-Jprima0*(Mu0/(4*Epsilonr));
}
//---------------FUNCIONES DE EQUILIBRIO-------------
double LatticeBoltzmann::feq(vector3D & Jprima0,vector3D & Eprima0,vector3D & B0,
			     double Epsilonr,double Mur,
			     int r,int p,int i,int j){
  double VdotJp=(v[p][i]*Jprima0),Epdote=(e[p][i][j]*Eprima0),Bdotb=(b[p][i][j]*B0),aux;
  if(r==0)
    aux=0.25*(0.25*VdotJp+Epsilonr*Epdote+0.5/Mur*Bdotb);
  if(r==1)
    aux=0.25*(0.25*VdotJp+Epdote+0.5*Bdotb);
  return aux;
}
double LatticeBoltzmann::feq0(double rhoc0){
  return rhoc0;
}

//-------------------FUNCIONES DE SIMULACION----------------------------
void LatticeBoltzmann::Inicie(void){
  int ix,iy,iz,r,p,i,j; double Sigma,Mur,Epsilonr,Prefactor;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0;
  double E00=1.0,alp0=0.01,iz0=40; double B00=E00/C;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
	//Calculo las constantes
	Sigma=sigma(ix,iy,iz); Mur=mur(ix,iy,iz); Epsilonr=epsilonr(ix,iy,iz);
	Prefactor=prefactor(Epsilonr,Sigma);
	//Impongo los campos
	rhoc0=0; Jprima0.cargue(0,0,0);
	B0.cargue(0,B00*exp(-alp0*(iz-iz0)*(iz-iz0)),0);
	Eprima0.cargue(E00*exp(-alp0*(iz-iz0)*(iz-iz0)),0,0);
	//Hago que los valores iniciales sean los de equilibrio
	for(r=0;r<2;r++)
	  for(p=0;p<3;p++)
	    for(i=0;i<4;i++)
	      for(j=0;j<2;j++)
		fnew[ix][iy][iz][r][p][i][j]=f[ix][iy][iz][r][p][i][j]=
		  feq(Jprima0,Eprima0,B0,Epsilonr,Mur,r,p,i,j);
	f0new[ix][iy][iz]=f0[ix][iy][iz]=feq0(rhoc0);
      }
}
void LatticeBoltzmann::Colisione(void){
  int ix,iy,iz,r,p,i,j; double Sigma,Mur,Epsilonr,Prefactor;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
	//Calculo las constantes
	Sigma=sigma(ix,iy,iz); Mur=mur(ix,iy,iz); Epsilonr=epsilonr(ix,iy,iz);
	Prefactor=prefactor(Epsilonr,Sigma);
	//Calculo los campos
	rhoc0=rhoc(ix,iy,iz,false); D0=D(ix,iy,iz,false); B0=B(ix,iy,iz,false);
	E0=E(D0,Epsilonr); H0=H(B0,Mur);
	Jprima0=E0*Prefactor; Eprima0=E0-Jprima0*(Mu0/(4*Epsilonr)); 
	//Hago la evolucion de la ec. de Boltzmann BGK
	for(r=0;r<2;r++)
	  for(p=0;p<3;p++)
	    for(i=0;i<4;i++)
	      for(j=0;j<2;j++)
		fnew[ix][iy][iz][r][p][i][j]=UmUTau*f[ix][iy][iz][r][p][i][j]
		  +UTau*feq(Jprima0,Eprima0,B0,Epsilonr,Mur,r,p,i,j);
	f0new[ix][iy][iz]=UmUTau*f0[ix][iy][iz]+UTau*feq0(rhoc0);
      }
}
void LatticeBoltzmann::Adveccione(void){
  int ix,iy,iz,r,p,i,j;
  for(ix=0;ix<Lx;ix++) 
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){//para cada celda
	for(r=0;r<2;r++)
	  for(p=0;p<3;p++)
	    for(i=0;i<4;i++)
	      for(j=0;j<2;j++)
		f[(ix+V[0][p][i]+Lx)%Lx][(iy+V[1][p][i]+Ly)%Ly][(iz+V[2][p][i]+Lz)%Lz][r][p][i][j]=
		  fnew[ix][iy][iz][r][p][i][j];
	f0new[(ix+V0[0]+Lx)%Lx][(iy+V0[1]+Ly)%Ly][(iz+V0[2]+Lz)%Lz]=
	  f0[ix][iy][iz];
      }
}
void LatticeBoltzmann::Muestre(void){
  int ix=0,iy=0,iz,r,p,i,j; double Sigma,Mur,Epsilonr,Prefactor;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0; double E2,B2;
  for(iz=0;iz<Lz;iz++){
    //Calculo las constantes
    Sigma=sigma(ix,iy,iz); Mur=mur(ix,iy,iz); Epsilonr=epsilonr(ix,iy,iz);
    Prefactor=prefactor(Epsilonr,Sigma);
    //Calculo los campos
    rhoc0=rhoc(ix,iy,iz,true); D0=D(ix,iy,iz,true); B0=B(ix,iy,iz,true);
    E0=E(D0,Epsilonr); H0=H(B0,Mur);
    Jprima0=E0*Prefactor; Eprima0=E0-Jprima0*(Mu0/(4*Epsilonr)); 
    //Imprimo
    E2=norma2(Eprima0); B2=norma2(B0);
    cout<<iz<<" "<<0.5*(Epsilonr*E2+B2/Mur)<<endl;
  }
}


int main(){
  LatticeBoltzmann Pulso;
  int t, tmax=140;
  
  Pulso.Inicie();
  
  for(t=0;t<tmax;t++){
    Pulso.Colisione();
    Pulso.Adveccione();
  }
  
  Pulso.Muestre();
  
  return 0;
}
