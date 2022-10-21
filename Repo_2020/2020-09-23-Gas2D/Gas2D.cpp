// Simular el movimiento de 2 planetasN moléculas en un gas 2D por PEFRL
#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h" 
using namespace std;


//------------------------------Declarar constantes------------------
const double K=6.0e4;
const double Lx=60, Ly=60;
const int Nx=10, Ny=10;
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
  vector3D  r, V, F;   double m,R;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void BorreFuerza(void){F.cargue(0,0,0);};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void Mueva_r(double dt, double coeficiente);
  void Mueva_V(double dt, double coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();};  //inline
  double Gety(void){return r.y();}; //inline

  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
  r.cargue(x0,y0,0); V.cargue(Vx0,Vy0,0); m=m0;  R=R0; 
}

void Cuerpo::Mueva_r(double dt, double coeficiente){
  r+=V*dt*coeficiente;
}

void Cuerpo::Mueva_V(double dt, double coeficiente){
  V+=(F*dt*coeficiente)/m;
}

void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}

//---------------------------Clase Colisionador-----------------------

class Colisionador{

private:

public:
  void CalculeFuerzas(Cuerpo * Molecula);
  void CalculeFuerzaEntre(Cuerpo & Molecula1, Cuerpo & Molecula2);
};

void Colisionador::CalculeFuerzas(Cuerpo * Molecula){

  int i, j;
  //borrar todas las fuerzas
  
  for(i=0; i<N+4; i++)
  Molecula[i].BorreFuerza();
  
  //Calcular todas las fuerzas
  
  for (i=0; i<N+4; i++){
    for (j=i+1; j<N+4; j++){
      CalculeFuerzaEntre(Molecula[i], Molecula[j]);
    }
  }
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Molecula1, Cuerpo & Molecula2){
  vector3D r21=Molecula2.r-Molecula1.r;
  double d=norma(r21);
  double s=(Molecula1.R+Molecula2.R)-d;
  vector3D n= r21*(1.0/d);
  if (s>0){
    vector3D F2= K*std::pow(s,1.5)*n;
  Molecula2.AdicioneFuerza(F2); Molecula1.AdicioneFuerza(-1*F2);
  }
}


//-------------------------- Funciones de Animacion -------------------

void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'Gas2D.gif'"<<endl;
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
  Cuerpo Molecula[N+4];
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
  
    //----------------(x0, y0, Vx0, Vy0, m0, R0)
  //arriba
  Molecula[N].Inicie  (Lx/2, Ly+Rpared, 0, 0, Mpared, Rpared);
  //abajo
  Molecula[N+1].Inicie(Lx/2,   -Rpared, 0, 0, Mpared, Rpared);
  //derecha
  Molecula[N+2].Inicie(Lx+Rpared, Ly/2, 0, 0, Mpared, Rpared);
  //izquierda
  Molecula[N+3].Inicie(-Rpared,   Ly/2, 0, 0, Mpared, Rpared); 
 
  //---------------------------Inicializar las moléculas----------------
  
  //---------------------------(x0, y0, Vx0, Vy0, m0, R0)
  for(ix=0; ix<Nx; ix++){
    for(iy=0; iy<Ny; iy++){
      theta =2*M_PI*ran64.r();
      
      Molecula[Nx*iy+ix].Inicie(dx*(ix+1), dy*(iy+1),V0*std::cos(theta),
				      V0*std::sin(theta),m0,R0);      
    }
  }
  
  InicieAnimacion(); //Dibujar
  
  for(t=0, tdibujo=0; t<tmax; t+=dt, tdibujo+=dt){
    
    
    //Dibujar animacion
    if(tdibujo>tcuadro){
      
      InicieCuadro();
      for(int i=0; i<N; i++)
	Molecula[i].Dibujese();
      TermineCuadro();
      
      
      // hacer un plot
      //std::cout<< Molecula[0].Getx()<<"\t"<<Molecula[0].Gety()<<"\t"<< Molecula[1].Getx()<<"\t"<<Molecula[1].Gety()<<endl;
       tdibujo=0;
    }
    //muevase por OMELYAN PEFRL
    
    for(i=0; i<N; i++)Molecula[i].Mueva_r(dt,E);
    Hertz.CalculeFuerzas(Molecula);  for(i=0; i<N; i++)Molecula[i].Mueva_V(dt,coeficiente1);
    for(i=0; i<N; i++)Molecula[i].Mueva_r(dt,X);
    Hertz.CalculeFuerzas(Molecula);  for(i=0; i<N; i++)Molecula[i].Mueva_V(dt,L);
    for(i=0; i<N; i++)Molecula[i].Mueva_r(dt,coeficiente2);
    Hertz.CalculeFuerzas(Molecula);  for(i=0; i<N; i++)Molecula[i].Mueva_V(dt,L);
    for(i=0; i<N; i++)Molecula[i].Mueva_r(dt,X);
    Hertz.CalculeFuerzas(Molecula);  for(i=0; i<N; i++)Molecula[i].Mueva_V(dt,coeficiente1);
    for(i=0; i<N; i++)Molecula[i].Mueva_r(dt,E);
    
  }   
  return 0;
}


