// SIMULACIÓN DEL MOVIMIENTO DE N PARTÍCULAS BAJO EL INFLUJO DE UNA FUERZA DE LENNARD JONES Y CONSIDERANDO FUERZAS REPULSIVAS GENERADAS POR LAS PAREDES  

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include "vector.h"
#include "Random64.h"
#include <sstream>
#include <fstream>
#include <vector>
using namespace std;
ofstream outfile;

//---- declarar constantes ---
const double K=1.0e4;
const double Lx=60, Ly=60;
const int Nx=5, Ny=5, N=Nx*Ny;

//Constantes de la Fuerza de Lennard Jones
const double E=1.0, r0=10;

const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);
int choque=0;

//--- declarar clases -----
class Cuerpo;
class Colisionador;

//---- interface e implementacion de clases ----
//---- clase cuerpo ---
class Cuerpo{
private:
  vector3D r,V,F; double m,R; 
public:
  int choquex=0, choquey=0; 
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void BorreFuerza(){F.load(0,0,0);};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
  double GetVx(void){return V.x();};
  double GetVy(void){return V.y();};
  void Choquex(void){choquex=1;};
  void Choquey(void){choquey=1;};
  void ReiniciarChoque(void){
    choquex=0;
    choquey=0;
  };
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
  outfile <<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
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
    Particula1.Choquex();
  }
  
  //Pared de la derecha
  
  if((x+R)>Lx){
    h=R-(Lx-x);
    n.load(1,0,0);
    vector3D F=n*(K*pow(h,1.5));  
    Particula1.AdicioneFuerza(F*(-1));
    Particula1.Choquex();
  }
  
  //Pared de abajo
  
  if(y<R){
    h=R-y;
    n.load(0,1,0);
    vector3D F=n*(K*pow(h,1.5));  
    Particula1.AdicioneFuerza(F);
    Particula1.Choquey();
  }
  
  //Pared de arriba
  
  if((y+R)>Ly){
    h=R-(Ly-y);
    n.load(0,1,0);
    vector3D F=n*(K*pow(h,1.5));  
    Particula1.AdicioneFuerza(F*(-1));    
    Particula1.Choquey();
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
  outfile <<"set terminal gif animate"<<endl;
  outfile <<"set output 'gas2d.gif'"<<endl;
  outfile <<"unset key"<<endl;
  outfile <<"set xrange[-10:"<<Lx+10<<"]"<<endl; 
  outfile <<"set yrange[-10:"<<Ly+10<<"]"<<endl;
  outfile <<"set size ratio -1"<<endl;
  outfile <<"set parametric"<<endl;
  outfile <<"set trange [0:7]"<<endl;
  outfile <<"set isosamples 12"<<endl;
}

void InicieCuadro(void){
    outfile <<"plot 0,0 ";
    outfile <<" , "<<Lx/7<<"*t,0";        //pared de abajo
    outfile <<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    outfile <<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    outfile <<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}

void TermineCuadro(void){
    outfile <<endl;
}

//-----------  Programa Principal --------------  
int main(void){
  Cuerpo Particula[N];
  Colisionador Hertz;
  Crandom ran64(1);
  vector<double> kbtvect{2,3,4,5,6,7,8,10,15,20};
  double m0=1.0, R0=2.5, kT, V0;
  int i,ix,iy;
  double t,tdibujo,tmax=120, tcuadro=tmax/1000,dt=0.001;
  double Theta;
  double total_change=0;
  int count[2*N]={0};
  int count3[N]={0}; //save wall hits for each

  outfile.open("T1_P5_i.dat");
  //InicieAnimacion(); 
 
  //Inicializar las moléculas
  for(double u : kbtvect){
    kT=u;
    V0=sqrt(2*kT/m0);
    total_change=0;
  for(ix=0;ix<Nx;ix++){
    for(iy=0;iy<Ny;iy++){
      Theta=2*M_PI*ran64.r();//el ángulo de cada molécula respecto a x es aleatorio 
      //--------------------(       x0,        y0,           Vx0,           Vy0, m0, R0)
      Particula[Nx*iy+ix].Inicie((ix+1)*10, (iy+1)*10, V0*cos(Theta), V0*sin(Theta), m0, R0);
    }
  }
  
  for(t=0 ; t<tmax ; t+=dt){

    double Vy[N], Vx[N];
    for(i=0;i<N;i++){
      Vy[i]=Particula[i].GetVy();
      Vx[i]=Particula[i].GetVx();
    }

    //cout<<Vx[0]<<endl;
    for(i=0;i<N;i++)Particula[i].ReiniciarChoque();
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

    double change=0;
    for(i=0;i<N;i++){
      if(Particula[i].choquex==1){
           
	if(t>80){change+=fabs(Vx[i]-Particula[i].GetVx());}
      }
      if(Particula[i].choquey==1){
	
	if(t>80){change+=fabs(Vy[i]-Particula[i].GetVy());}
      }
      else{count[N+i]=0;}
    }
    
    total_change+=change;
  }
  
  double Press=(total_change/(40*60));
  outfile<<u<<" "<<Press<<endl;
}
  outfile.close();
  return 0;
}
