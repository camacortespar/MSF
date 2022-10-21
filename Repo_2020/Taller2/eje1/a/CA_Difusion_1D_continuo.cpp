
 #include <cmath>
#include "Random64.h"
using namespace std;

const int Lx=256;
const int Ly=256;
const double p=0.25;
const double p0=0.25;
const int Q=4;

//--------------------- Clase LatticeGas ------------
class LatticeGas{
private:
  
  double w[Q];
  int V[2][Q]; //V[0][i]=Vi_x ,  V[1][i]=Vi_y                                   
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; // f[ix][iy][i]
  
public:
  LatticeGas(void);
  void Inicie(double mu0,double sigma0,int N);
  double rho(int ix, int iy, bool ShowNew);
  void Colisione(void);
  void Adveccione(void);
  void Show(bool ShowNew);
  double Sigma2(void);
};  
LatticeGas::LatticeGas(void){
  //cargar los pesos                                                            
  w[0]=p0; w[1]=w[3]=p; w[2]=1-p0-2*p;
  
  //cargar los vectores velocidad                                               
  V[0][0]=1;  V[0][1]=0;  V[0][2]=-1;  V[0][3]=0; 
  V[1][0]=0;  V[1][1]=1;  V[1][2]=0;   V[1][3]=-1;
//derecha     arriba      izquierda    abajo
  
  for(int ix=0;ix<Lx;ix++)                                   
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
	f[ix][iy][i]=fnew[ix][iy][i]=0;
}
void LatticeGas::Inicie(double mu0,double sigma0,int N){
  for(int ix=0;ix<Lx;ix++)
     for(int iy=0;iy<Ly;iy++)
       for(int i=0;i<Q;i++)
	 f[ix][iy][i]=N/(sigma0*sqrt(2*M_PI))*exp(-0.5*pow((ix-mu0)/sigma0,2)-0.5*pow((iy-mu0)/sigma0,2));
}  
double LatticeGas::rho(int ix,int iy, bool ShowNew){
  double suma; int i;
  for(suma=0,i=0;i<Q;i++)
    if(ShowNew)
      suma+=fnew[ix][iy][i];
    else
      suma+=f[ix][iy][i];
  return suma;
}  
void LatticeGas::Colisione(void){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++) //en cada dirección
	fnew[ix][iy][i]=w[0]*f[ix][iy][i]+w[1]*f[ix][iy][(i+1)%4]+w[2]*f[ix][iy][(i+2)%4]+w[3]*f[ix][iy][(i+3)%4];
}
void LatticeGas::Adveccione(void){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++)
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i]; //fronteras periódicas
}
void LatticeGas::Show(bool ShowNew){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      cout<<ix<<"\t"<<iy<<"\t"<<rho(ix,iy,ShowNew)<<endl;
}
  
double LatticeGas::Sigma2(void){
  double N,rprom,suma; int ix;
  //Rectificar N
  for(suma=0,ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      suma+=rho(ix,iy,false);
  N=suma;
  //Calcular ixprom
  for(suma=0,ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      suma+=std::hypot(ix*rho(ix,iy,false),iy*rho(ix,iy,false));
  rprom=suma/N;
  //Calcular sigma2
  for(suma=0,ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++) 
      suma+=pow(hypot(ix,iy)-rprom,2)*rho(ix,iy,false);
  return suma/(N-1);
}
//------------------- Funciones Globales ------------

int main(void){
  LatticeGas Difusion;
  int t,tmax=350;
  int N=2400;

  //Inicie
  Difusion.Inicie(0.5*Lx,16,N);
  //Corra
  
  for(t=0;t<tmax;t++){
    Difusion.Colisione();
    Difusion.Adveccione();
    cout<<t<<" "<<Difusion.Sigma2()<<endl;
  }
  //  Difusion.Show(false);
  
  return 0;
}  


