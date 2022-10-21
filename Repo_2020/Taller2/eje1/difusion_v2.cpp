

#include <cmath>
#include "Random64.h"
using namespace std;

const int Q=4;
const int Lx=256;
const int Ly=256;
const double p=0.25;
const double p0=0.25;

//--------------------- Clase LatticeGas ------------
class LatticeGas{
private:
  double w[Q];
  int V[2][Q]; //V[i], i=0 (derecha), i=1 (izquierda)
  int n[Lx][Ly][Q];  int nnew[Lx][Ly][Q];
public:
  LatticeGas(void);
  void Inicie(int N,double mu,double sigma,Crandom & ran64);
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  double GetSigma2(void);
};
LatticeGas::LatticeGas(void){
  //cargar los pesos                                          
  w[0]=p0; w[1]=w[3]=p; w[2]=1-p0-2*p;
  //cargar los vectores v                                    
  V[0][0]=1;  V[0][1]=0;  V[0][2]=-1;  V[0][3]=0;
  V[1][0]=0;  V[1][1]=1;  V[1][2]=0;   V[1][3]=-1;
//derecha     arriba      izquierda    abajo 
  for(int ix=0;ix<Lx;ix++)
    for(int iy;iy<Ly;iy++)
      n[ix][iy][0]=n[ix][iy][1]=nnew[ix][iy][0]=nnew[ix][iy][1]=0;
}
void LatticeGas::Inicie(int N,double mu,double sigma,Crandom & ran64){
  int ix,iy,i;
  do{
    ix=(int) ran64.gauss(mu,sigma); if(ix<0) ix=0; if(ix>(Lx-1)) ix=Lx-1;
    iy=(int) ran64.gauss(mu,sigma); if(iy<0) iy=0; if(iy>(Ly-1)) iy=Ly-1;
    i=(int) Q*ran64.r();
    if(n[ix][iy][i]==0)
      {n[ix][iy][i]=1; N--;}
    }while(N>0);
}   
void LatticeGas::Colisione(Crandom & ran64){
  double proba=ran64.r();
  for(int ix = 0; ix < Lx; ix++){
    for(int iy = 0; iy < Ly; iy++){                                                 //para cada celda                                                        
      if(proba >= 0 && proba < p0){                           //con probabilidad p0                                                            
                                nnew[ix][iy][0] = n[ix][iy][0];
				nnew[ix][iy][1] = n[ix][iy][1];
                                nnew[ix][iy][2] = n[ix][iy][2];
				nnew[ix][iy][3] = n[ix][iy][3];
                        }                                                                                                       //dejelo igual.                          
                        else if(proba >= p0 && proba < p){                      //con probabilidad p                                                             
                                nnew[ix][iy][0] = n[ix][iy][3];
                                nnew[ix][iy][1] = n[ix][iy][0];
                                nnew[ix][iy][2] = n[ix][iy][1];
                                nnew[ix][iy][3] = n[ix][iy][2];
                        }                                                                                                               //girar 180 grados.              
			else if(proba >= p && proba < (1 - 2*p - p0)){  //con 1 -2*p - p0                                                                        
                                nnew[ix][iy][0] = n[ix][iy][2];
                                nnew[ix][iy][1] = n[ix][iy][3];
                                nnew[ix][iy][2] = n[ix][iy][0];
                                nnew[ix][iy][3] = n[ix][iy][1];
                        }                                                                                                       //girar 270 grados.                      
                        else{  //con probabilidad 1-p                      
                                nnew[ix][iy][0] = n[ix][iy][1];
                                nnew[ix][iy][1] = n[ix][iy][2];
                                nnew[ix][iy][2] = n[ix][iy][3];
                                nnew[ix][iy][3] = n[ix][iy][0];
                        }
                }
        }
}
   
void LatticeGas::Adveccione(void){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int iy = 0; iy < Ly; iy++)
      for(int i=0;i<Q;i++)
	n[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=nnew[ix][iy][i]; //fronteras periódicas
}
double LatticeGas::GetSigma2(void){
  double ixprom,suma,N; int ix,iy=0;
  //Calcular N
  for(N=0,ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      N+=(n[ix][iy][0]+n[ix][iy][1]+n[ix][iy][2]+n[ix][iy][3]);
  //Calcular ixprom
  for(suma=0,ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      suma+=(n[ix][iy][0]+n[ix][iy][1])*ix(+n[ix][iy][2]+n[ix][iy][3])*iy;
  ixprom=suma/N;
  //Calcular sigma2
  for(suma=0,ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      suma+=(n[ix][iy][0]+n[ix][iy][1])*pow(ix-ixprom,2)+(n[ix][iy][2]+n[ix][iy][3])*pow(iy-ixprom,2);
  return suma/(N-1);
}  
//------------------- Funciones Globales ------------

int main(void){
  LatticeGas Difusion;
  Crandom ran64(0);
  int t,tmax=350;
  int N=2400; double mu=Lx/2,sigma=16;

  //Inicie
  Difusion.Inicie(N,mu,sigma,ran64);
  //Corra
  for(t=0;t<tmax;t++){
    Difusion.Colisione(ran64);
    Difusion.Adveccione();
    if(t%10==0)
      cout<<t<<" "<<Difusion.GetSigma2()<<endl;
  }
  
  return 0;
}  
