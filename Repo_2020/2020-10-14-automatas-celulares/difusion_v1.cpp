#include <cmath>
#include "Random64.h"
using namespace std;

const int Lx=256;
const double p=0.5;

//--------------------- Clase LatticeGas ------------
class LatticeGas{
private:
  int V[2]; //V[i], i=0 (derecha), i=1 (izquierda)
  int n[Lx][2];  int nnew[Lx][2];
public:
  LatticeGas(void);
  void Inicie(Crandom & ran64);
  void Show(bool ShowNew);
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  int DondeEstaLaBolita(void);
};
LatticeGas::LatticeGas(void){
  V[0]=1;  V[1]=-1;
  for(int ix=0;ix<Lx;ix++)
    n[ix][0]=n[ix][1]=nnew[ix][0]=nnew[ix][1]=0;
}
void LatticeGas::Inicie(Crandom & ran64){
  int ix=Lx/2;
  int i=(int) 2*ran64.r();
  n[ix][i]=1;
}  
void LatticeGas::Show(bool ShowNew){
  for(int i=0;i<2;i++){
    for(int ix=0;ix<Lx;ix++)
      if(ShowNew) cout<<nnew[ix][i]; else cout<<n[ix][i];
    cout<<endl;
  }
  cout<<endl;
}  
void LatticeGas::Colisione(Crandom & ran64){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    if (ran64.r()<p) //con probabilidad p
      {nnew[ix][0]=n[ix][0]; nnew[ix][1]=n[ix][1];}//dejelo igual;
    else  //con probabilidad 1-p
      {nnew[ix][0]=n[ix][1]; nnew[ix][1]=n[ix][0];}//intercámbielos;      
}
void LatticeGas::Adveccione(void){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int i=0;i<2;i++)
      n[(ix+V[i]+Lx)%Lx][i]=nnew[ix][i]; //fronteras periódicas
}
int LatticeGas::DondeEstaLaBolita(void){
  int ix=0;
  while(n[ix][0]==0 && n[ix][1]==0)
    ix++;
  return ix;
}
//------------------- Funciones Globales ------------
const int N=100;

double Sigma2(LatticeGas * Difusion){
  double ixprom,suma; int c;
  //Calcular ixprom
  for(suma=0,c=0;c<N;c++)
    suma+=Difusion[c].DondeEstaLaBolita();
  ixprom=suma/N;
  //Calcular ixprom
  for(suma=0,c=0;c<N;c++)
    suma+=pow(Difusion[c].DondeEstaLaBolita()-ixprom,2);
  return suma/(N-1);
}
  
int main(void){
  LatticeGas Difusion[N];
  Crandom ran64(2);
  int c;
  int t,tmax=10000;
  
  //Inicie
  for(c=0;c<N;c++)
    Difusion[c].Inicie(ran64);
  //Corra
  for(t=0;t<tmax;t++){
    for(c=0;c<N;c++){
      Difusion[c].Colisione(ran64);
      Difusion[c].Adveccione();
    }
    cout<<t<<" "<<Sigma2(Difusion)<<endl;
  }
  
  return 0;
}  
