#include <iostream>
#include <cmath>
#include "../../vector.h"
#include <fstream>
/////////////////////////////////////////////////////////////////
// Dos Cuerpos Celestes - metodo de integracion PEFRL          //
//                                                             //
// En este programa se simula y anima el movimiento ente dos   //
// planetas con un observador ajeno al sistema, a este último  //
// se le llama colisionador o interaccionador.                 //
// Los objetos celestes son el Sol y Júpiter a quienes vamos a //
// asignar las masas y radios como razones entre ellos:        //
// m_sol = 1047 m_jupiter    ,  r_sol = 9.73 r_jupiter         //
// m_jupiter = 1             ,  r_jupiter = 1                  //
// El código esta "configurado" para realizar la animación del //
// movimiento de ambos cuerpos. Si se desea graficar sin       //
// animación comentar lineas 129 a 135 y descomentar lines     //
// 137 y 138. Para el gif activar las lineas 96 y 97.          //
/////////////////////////////////////////////////////////////////
const double G=1.0; // Constante gravitacional
const int N=2;      // Nos permitirá crear el arreglo Sol-Júpiter

// Coeficientes de PEFRL
const double Zeta   = 0.1786178958448091e00;
const double Lambda = -0.2123418310626054e00;
const double Chi    = -0.6626458266981849e-1;
const double Coeficiente1 =(1-2*Lambda)/2;
const double Coeficiente2 =1-2*(Chi+Zeta);

// Clases
class Cuerpo;
class Colisionador;
//------------------ Clase Cuerpo ------------------------------//
class Cuerpo{
private:
  vector3D r, V, F;
  double m,R;
public:
  void Inicie(double x0,double y0,double z0, double Vx0, double Vy0, double Vz0,double m0,double R0);
  void BorreFuerza(void);
  void SumeFuerza(vector3D F0);
  void Mueva_r(double dt, double coeficiente);
  void Mueva_v(double dt, double coeficiente);
  void Dibujese(void);             
  double Getx(void){return r.x();}; //Inline
  double Gety(void){return r.y();}; //Inline
  double Getz(void){return r.z();}; //Inline
  friend class Colisionador; //Cualquier especimen de la clase colisionador
                             // le puede meter la mano a la clase cuerpo
};
void Cuerpo::Inicie(double x0,double y0,double z0, double Vx0, double Vy0, double Vz0,double m0,double R0){
  r.load(x0,y0,z0);
  V.load(Vx0,Vy0,Vz0);
  m=m0; R=R0;
}
void Cuerpo::BorreFuerza(void){
  F.load(0,0,0);
}
void Cuerpo::SumeFuerza(vector3D F0){
  F+=F0;
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  vector3D r_new;
  r+=V*dt*coeficiente;
} 
void Cuerpo::Mueva_v(double dt,double coeficiente){
  vector3D r_new;
  V+=F*dt*coeficiente/m;
}
void Cuerpo::Dibujese(void){
  std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)"; 
}
//------------------ Class Colisionador ------------------------//
class Colisionador {
private:
public:
  void CalculeFuerza(Cuerpo * Planeta);
  void CalculeFuerzaEntre(Cuerpo &  Planeta1, Cuerpo & Planeta2);
};
void Colisionador::CalculeFuerza(Cuerpo * Planeta){
  int i,j;
  for(i=0;i<N;i++)
    Planeta[i].BorreFuerza();
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++)
      CalculeFuerzaEntre(Planeta[i],Planeta[j]);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Planeta1, Cuerpo &Planeta2){
  vector3D r21,n,F1; double d21,F;
  r21 = Planeta2.r - Planeta1.r; d21=r21.norm(); n=r21/d21;
  F = G*Planeta1.m*Planeta2.m*pow(d21,-2.0);
  F1=F*n;
  Planeta1.SumeFuerza(F1); Planeta2.SumeFuerza(F1*(-1));
}
//------------------ Funciones globales ------------------------//
void InicieAnimacion(void){
  std::cout<<"set terminal gif animate"<<std::endl;  // quitar si quiero imprimir en el termminal
  std::cout<<"set output 'DosPlanetas.gif'"<<std::endl; // quitar si quiero imprimir en el termminal
  std::cout<<"set grid"<<std::endl;
  std::cout<<"show grid"<<std::endl;
  std::cout<<"unset key"<<std::endl;   // quitar si quiero imprimir un gif
  std::cout<<"set xrange[-1100:1100]"<<std::endl;
  std::cout<<"set yrange[-1100:1100]"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set parametric"<<std::endl; // la variable parametrica es t 
  std::cout<<"set trange[0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;
}
void InicieCuadro(void){
  std::cout<<"plot 0,0 ";
}
void TermineCuadro(void){
  std::cout<<std::endl;
}
int main(){
  Cuerpo Planeta[N];
  Colisionador Newton;
  double m1=1047, m2=1, r=1000;
  double M=m1+m2, x1=-m2*r/M, x2=m1*r/M;
  double omega = sqrt(G*M/pow(r,3.0)), T=2*M_PI/omega, V1=omega*x1, V2=omega*x2;
  double t,tmax=1.0*T, dt=0.01;
  double tdibujo, tcuadro=T/100;
  int i;
  std::ofstream outfile;

  outfile.open("S-J_original_axis.dat");
  Planeta[0].Inicie(x1,0,0,0,V1,0,m1,97.3);
  Planeta[1].Inicie(x2,0,0,0,V2,0,m2,10.0);
  
  InicieAnimacion();
  for(t=0,tdibujo=0; t<=tmax; t+=dt, tdibujo+=dt){
    //Dibujar
    if(tdibujo>tcuadro){
      InicieCuadro();
      for(i=0;i<N;i++) 	Planeta[i].Dibujese();
      TermineCuadro();
      tdibujo=0; 
      }
  
    //outfile<<Planeta[0].Getx()<<" "<<Planeta[0].Gety()<<" "
    //   <<Planeta[1].Getx()<<" "<<Planeta[1].Gety()<<std::endl;
    // Haga el movimiento y los calculos por PEFRL
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt, Zeta);            // 1 
    Newton.CalculeFuerza(Planeta);
    for(i=0;i<N;i++) Planeta[i].Mueva_v(dt,Coeficiente1);     // 2
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt, Chi);             // 3
    Newton.CalculeFuerza(Planeta);     
    for(i=0;i<N;i++) Planeta[i].Mueva_v(dt,Lambda);           // 4
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Coeficiente2);     // 5
    Newton.CalculeFuerza(Planeta);     
    for(i=0;i<N;i++) Planeta[i].Mueva_v(dt,Lambda);           // 4
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt, Chi);             // 3
    Newton.CalculeFuerza(Planeta);
    for(i=0;i<N;i++) Planeta[i].Mueva_v(dt,Coeficiente1);     // 2
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt, Zeta);            // 1
  }
  //std::cout<<T<<std::endl;
  return 0;
}
