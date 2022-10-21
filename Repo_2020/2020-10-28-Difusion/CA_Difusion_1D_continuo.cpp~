#include <cmath>
#include "Random64.h"
using namespace std;

const int Lx=1024;
const double p=0.5;

//--------------------- Clase LatticeGas ------------
class LatticeGas{
private:
  int V[2]; //V[i], i=0 (derecha), i=1 (izquierda)
  double f[Lx][2], fnew[Lx][2];
public:
  LatticeGas(void);
  void Inicie(double mu0,double sigma0,int N);
  double rho(int ix,bool ShowNew);
  void Colisione(void);
  void Adveccione(void);
  void Show(bool ShowNew);
  double Sigma2(void);
};  
LatticeGas::LatticeGas(void){
  V[0]=1;  V[1]=-1;
  for(int ix=0;ix<Lx;ix++)
    f[ix][0]=f[ix][1]=fnew[ix][0]=fnew[ix][1]=0;
}
void LatticeGas::Inicie(double mu0,double sigma0,int N){
  for(int ix=0;ix<Lx;ix++)
    f[ix][0]=f[ix][1]=N/(sigma0*sqrt(2*M_PI))*exp(-0.5*pow((ix-mu0)/sigma0,2));
}  
double LatticeGas::rho(int ix, bool ShowNew){
  if(ShowNew)
    return fnew[ix][0]+fnew[ix][1];
  else
    return f[ix][0]+f[ix][1];
}  
void LatticeGas::Colisione(void){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int i=0;i<2;i++) //en cada dirección
      fnew[ix][i]=p*f[ix][i]+(1-p)*f[ix][(i+1)%2];
}
void LatticeGas::Adveccione(void){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int i=0;i<2;i++)
      f[(ix+V[i]+Lx)%Lx][i]=fnew[ix][i]; //fronteras periódicas
}
void LatticeGas::Show(bool ShowNew){
  for(int ix=0;ix<Lx;ix++)
    cout<<ix<<" "<<rho(ix,ShowNew)<<endl;
}
  
double LatticeGas::Sigma2(void){
  double N,ixprom,suma; int ix;
  //Rectificar N
  for(suma=0,ix=0;ix<Lx;ix++)
    suma+=rho(ix,false);
  N=suma;
  //Calcular ixprom
  for(suma=0,ix=0;ix<Lx;ix++)
    suma+=ix*rho(ix,false);
  ixprom=suma/N;
  //Calcular sigma2
  for(suma=0,ix=0;ix<Lx;ix++)
    suma+=pow(ix-ixprom,2)*rho(ix,false);
  return suma;
}
//------------------- Funciones Globales ------------

int main(void){
  LatticeGas Difusion;
  int t,tmax=40000;
  int N=1;

  //Inicie
  Difusion.Inicie(0.5*Lx,0.02*Lx,N);
  //Corra
  
  for(t=0;t<tmax;t++){
    Difusion.Colisione();
    Difusion.Adveccione();
    cout<<t<<" "<<Difusion.Sigma2()<<endl;
  }
  
  //  Difusion.Show(false);
  
  return 0;
}  
