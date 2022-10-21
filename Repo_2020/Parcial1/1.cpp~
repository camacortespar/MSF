// Simular el movimiento de N granos 2D por PEFRL con grvedad y choques
#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h" 
using namespace std;


//------------------------------Declarar constantes------------------
const double g=9.8, K=1.0e4, Gamma=50, Kcundall=10, MU=0.4;
const double Lx=60, Ly=60;
const int Nx=1, Ny=1;
const int N=Nx*Ny;

const double E=0.1786178958448091e00;
const double L=-0.2123418310626054e0;
const double X=-0.6626458266981849e-1;
const double coeficiente1=(1-2*L)/2;
const double coeficiente2=(1-2*(X+E))*2/2;

//-----------------------------Declarar clases----------------------
class Cuerpo;
class Colisionador;

//-----------------------------Implementar clases--------------------

//-----------------------------Clase Cuerpo-------------------------
class Cuerpo{
private:
  vector3D  r, V, F;   double m,R, theta, omega, tau, I;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0, double theta0, double omega0);
  void BorreFuerza(void){F.cargue(0,0,0); tau=0;};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void AdicioneTorque(double tau0){tau+=tau0;};
  void Mueva_r(double dt, double coeficiente);
  void Mueva_V(double dt, double coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();};  //inline
  double Gety(void){return r.y();}; //inline
   double Gettheta(void){return theta;};

  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0, double theta0, double omega0){
  r.cargue(x0,y0,0); V.cargue(Vx0,Vy0,0); m=m0;  R=R0; theta=theta0; omega=omega0; I=2.0*m*R*R/5;
}

void Cuerpo::Mueva_r(double dt, double coeficiente){
  r+=V*dt*coeficiente;  theta+=omega*dt*coeficiente;
}

void Cuerpo::Mueva_V(double dt, double coeficiente){
  V+=(F*dt*coeficiente)/m; omega+=(tau*dt*coeficiente)/I;
}

void Cuerpo::Dibujese(void){ 
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";    //dibujar la circunferencia
  cout<<" , "<<r.x()<<"+"<<R*cos(theta)/7<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7<<"*t";  //dibujar una linea 
}

//---------------------------Clase Colisionador-----------------------

class Colisionador{

private:

public:
  void CalculeFuerzas(Cuerpo * Grano);
  void CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2);
};

void Colisionador::CalculeFuerzas(Cuerpo * Grano){

  int i, j; vector3D Fg;
  //borrar todas las fuerzas
  
  for(i=0; i<N+4; i++)
  Grano[i].BorreFuerza();
  
//sumar el peso 
  
  for(i=0; i<N; i++){
    Fg.cargue(0,-Grano[i].m*g,0);
    Grano[i].AdicioneFuerza(Fg);
  }
  
  //Calcular todas las fuerzas entre granos
  
  for (i=0; i<N+4; i++){
    for (j=i+1; j<N+4; j++){
      CalculeFuerzaEntre(Grano[i], Grano[j]);
    }
  }
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2){
  vector3D r21=Grano2.r-Grano1.r;
  double d=norma(r21);
  double s=(Grano1.R+Grano2.R)-d;
  vector3D n= r21*(1.0/d);

  if (s>0){
    vector3D F2= K*std::pow(s,1.5)*n;
  Grano2.AdicioneFuerza(F2); Grano1.AdicioneFuerza(-1*F2);
  }
}


//-------------------------- Funciones de Animacion -------------------

void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'Granos2D.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;
  cout<<"set yrange[-10:"<<Ly+10<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
  cout<<"plot 0,0 ";   
  cout<<" , "<<Lx/7<<"*t,0";           //pared de abajo
  cout<<" , "<<Lx/7<<"*t,"<<Ly;       //pared de arriba
  cout<<" ,0,"<<Ly/7<<"*t";          //pared izquierda
  cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared derecha
}
void TermineCuadro(void){
  cout<<endl;
}

//-----------------------  Programa Principal ------------------------

int main(void){
  Cuerpo Grano[N+4];
  Colisionador Hertz;
  Crandom ran64(1);
  
  double m0=1, R0=2, kT=10, V0=std::sqrt(2*kT/m0);
  int ix,iy, i;
  double dx=Lx/(Nx+1), dy=Ly/(Ny+1);
  double theta;
  
    double t, tdibujo, tmax=10*(Lx/V0), tcuadro=tmax/500, dt=0.001;

  
 //----------------------------Inicializar las parredes-----------------
  double Rpared=10000*Lx;
  double Mpared=100*m0;
  
    //----------------(x0, y0, Vx0, Vy0, m0, R0, theta0, omega0)
  //arriba
  Grano[N].Inicie  (Lx/2, Ly+Rpared, 0, 0, Mpared, Rpared, 0, 0);
  //abajo
  Grano[N+1].Inicie(Lx/2,   -Rpared, 0, 0, Mpared, Rpared, 0, 0);
  //derecha
  Grano[N+2].Inicie(Lx+Rpared, Ly/2, 0, 0, Mpared, Rpared, 0, 0);
  //izquierda
  Grano[N+3].Inicie(-Rpared,   Ly/2, 0, 0, Mpared, Rpared, 0, 0); 
 
  //---------------------------Inicializar los granos -------------------
  
  //---------------------------(x0, y0, Vx0, Vy0, m0, R0, theta0, omega0)
  for(ix=0; ix<Nx; ix++){
    for(iy=0; iy<Ny; iy++){
      theta =2*M_PI*ran64.r();
      
      Grano[Nx*iy+ix].Inicie(dx*(ix+1), dy*(iy+1),V0*std::cos(theta),
			     V0*std::sin(theta),m0,R0, 0, 1);      
    }
  }
  
  InicieAnimacion(); //Dibujar
  
  for(t=0, tdibujo=0; t<tmax; t+=dt, tdibujo+=dt){
    
    
    //Dibujar animacion
    if(tdibujo>tcuadro){
      
      InicieCuadro();
      for(int i=0; i<N; i++)
	Grano[i].Dibujese();
      TermineCuadro();
      
      
      // hacer un plot
      //std::cout<< Grano[0].Getx()<<"\t"<<Grano[0].Gety()<<"\t"<< Grano[1].Getx()<<"\t"<<Grano[1].Gety()<<endl;
       tdibujo=0;
    }
    //muevase por OMELYAN PEFRL
    
    for(i=0; i<N; i++)Grano[i].Mueva_r(dt,E);
    Hertz.CalculeFuerzas(Grano);  for(i=0; i<N; i++)Grano[i].Mueva_V(dt,coeficiente1);
    for(i=0; i<N; i++)Grano[i].Mueva_r(dt,X);
    Hertz.CalculeFuerzas(Grano);  for(i=0; i<N; i++)Grano[i].Mueva_V(dt,L);
    for(i=0; i<N; i++)Grano[i].Mueva_r(dt,coeficiente2);
    Hertz.CalculeFuerzas(Grano);  for(i=0; i<N; i++)Grano[i].Mueva_V(dt,L);
    for(i=0; i<N; i++)Grano[i].Mueva_r(dt,X);
    Hertz.CalculeFuerzas(Grano);  for(i=0; i<N; i++)Grano[i].Mueva_V(dt,coeficiente1);
    for(i=0; i<N; i++)Grano[i].Mueva_r(dt,E);
    
  }   
  return 0;
}


