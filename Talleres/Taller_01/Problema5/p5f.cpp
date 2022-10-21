// SIMULACIÓN DEL MOVIMIENTO DE N PARTÍCULAS BAJO EL INFLUJO DE UNA FUERZA DE LENNARD JONES Y CONSIDERANDO FUERZAS REPULSIVAS GENERADAS POR LAS PAREDES  

#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
using namespace std;

//---- declarar constantes ---
const double K=1.0e4;
const double Lx=60, Ly=120;
const int Nx=5, Ny=5, N=Nx*Ny;

//Constantes de la Fuerza de Lennard Jones
const double E=1.0, r0=10;

const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//--- declarar clases -----
class Cuerpo;
class Colisionador;

//---- interface e implementacion de clases ----
//---- clase cuerpo ---
class Cuerpo{
private:
  vector3D r,V,F; double m,R; 
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void BorreFuerza(){F.load(0,0,0);};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
  double GetV(void){return V.x();};
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0;
} 
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt);  
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m); 
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
  // cout<<" , "<<r.x()<<"+"<<R*cos(theta)/7<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7<<"*t";
}


//--- clase Colisionador ----
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Particula);
  void CalculeFuerzaEntre(Cuerpo & Particula1, Cuerpo & Particula2);
  void CalculeFuerzaPared(Cuerpo & Particula1);
};

void Colisionador::CalculeFuerzas(Cuerpo * Particula){
  int i,j;
  
  //--- Borrar todas las fuerzas ---
  for(i=0;i<N;i++)
    Particula[i].BorreFuerza();
  
  //--- Sumar la fuerza de la pared ---
  for(i=0;i<N;i++){
    CalculeFuerzaPared(Particula[i]);
  }
  //--- Calcular Fuerzas entre pares de granos ---
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++){
      CalculeFuerzaEntre(Particula[i], Particula[j]);
    }
}


//Implementación de la fuerzas repulsivas generadas por cada pared

void Colisionador::CalculeFuerzaPared(Cuerpo &Particula1){
  double x=Particula1.Getx(), y=Particula1.Gety();
  vector3D r=Particula1.r;
  double h, d=r.norm(), R=Particula1.R;
  vector3D n;
  
  //Pared de la izquierda
  
  if(x<R){
    h=R-x;
    n.load(1,0,0);
    vector3D F=n*(K*pow(h,1.5));  
    Particula1.AdicioneFuerza(F);
  }
  
  //Pared de la derecha
  
  if((x+R)>Lx){
    h=R-(Lx-x);
    n.load(1,0,0);
    vector3D F=n*(K*pow(h,1.5));  
    Particula1.AdicioneFuerza(F*(-1));
  }
  
  //Pared de abajo
  
  if(y<R){
	 h=R-y;
	 n.load(0,1,0);
	 vector3D F=n*(K*pow(h,1.5));  
    Particula1.AdicioneFuerza(F);
  }
  
  
  //Pared de arriba
  
  if((y+R)>Ly){
    h=R-(Ly-y);
    n.load(0,1,0);
    vector3D F=n*(K*pow(h,1.5));  
    Particula1.AdicioneFuerza(F*(-1));
  }
  
}



//Fuerza de Lennard Jones entre dos moleculas

void Colisionador::CalculeFuerzaEntre(Cuerpo & Particula1,Cuerpo & Particula2){
  vector3D r21=Particula2.r-Particula1.r;
  double d=r21.norm();
    vector3D n=r21*(1.0/d);
    vector3D F2=n*12*E/d*(pow(r0/d,12)-pow(r0/d,6));
    Particula2.AdicioneFuerza(F2);   Particula1.AdicioneFuerza(F2*(-1));
  
 
}



//----------------- Funciones de Animacion ----------

void InicieAnimacion(void){
  // cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'Gas2D.gif'"<<endl;
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
    cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}

void TermineCuadro(void){
    cout<<endl;
}

//-----------  Programa Principal --------------  
int main(void){
  Cuerpo Particula[N];
  Colisionador Hertz;
  Crandom ran64(1);
  double m0=1.0, R0=2.5, kT=20.0, V0=sqrt(2*kT/m0);
  int i,ix,iy;
  double t,tdibujo,tmax=280, tcuadro=tmax/1000,dt=0.001;
  double Theta; double suma=0;
 
  //Inicializar las moléculas
  
  for(ix=0;ix<Nx;ix++)
    for(iy=0;iy<Ny;iy++){
      Theta=2*M_PI*ran64.r();//el ángulo de cada molécula respecto a x es aleatorio 
				 //--------------------(       x0,        y0,           Vx0,           Vy0, m0, R0)
				 Particula[Nx*iy+ix].Inicie((ix+1)*10, (iy+1)*10, V0*cos(Theta), V0*sin(Theta), m0, R0);
    }
  
  for(t=0; t<tmax; t+=dt){
    
    if (80<=t){
    for(i=0;i<N;i++) {
      cout<<t<<"\t"<<Particula[i].GetV()<<endl;
    }
    
    }
;
    
    //--- Muevase por PEFRL ---
    for(i=0;i<N;i++)Particula[i].Mueva_r(dt,epsilon);
    Hertz.CalculeFuerzas(Particula);
    for(i=0;i<N;i++)Particula[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Particula[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Particula);
    for(i=0;i<N;i++)Particula[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Particula[i].Mueva_r(dt,chiepsilon);
    Hertz.CalculeFuerzas(Particula);
    for(i=0;i<N;i++)Particula[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Particula[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Particula);
    for(i=0;i<N;i++)Particula[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Particula[i].Mueva_r(dt,epsilon);  

  }   

 
  
  return 0;
}

  
