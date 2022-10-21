// Simular el movimiento de un balón por PEFRL 
#include <iostream>
#include <cmath>
#include "Vector.h" 
using namespace std;


const double GM=1.0;
const double E=0.1786178958448091e00;
const double L=-0.2123418310626054e0;
const double X=-0.6626458266981849e-1;

const double coeficiente1=(1-2*L)/2;
const double coeficiente2=(1-2*(X+E))/2;



class Cuerpo{
private:
  vector3D  r, V, F;   double m,R;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void CalculeFuerza(void);
  void Mueva_r(double dt, double coeficiente);
  void Mueva_V(double dt, double coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();};  //inline
  double Gety(void){return r.y();}; //inline
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
  r.cargue(x0,y0,0); V.cargue(Vx0,Vy0,0); m=m0;  R=R0; 
} 
void Cuerpo::CalculeFuerza(void){
  double aux = -GM*m*std::pow(norma2(r),-1.5);
  F=aux*r;
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
//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'Planeta.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-7:15]"<<endl;
  cout<<"set yrange[-8:8]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
  cout<<"plot 0,0 ";
}
void TermineCuadro(void){
  cout<<endl;
}

//-----------  Programa Principal --------------  
int main(void){
  Cuerpo Planeta;
  double r0=10, m0=1;
  double omega, V0, T;
  
  omega=std::sqrt(GM/(r0*r0*r0)); T=2*M_PI/omega; V0=r0*omega;
  
  double t, tdibujo, tmax=3.3*T, tcuadro=T/40, dt=0.1;
  
  //------------(x0,y0,Vx0,Vy0, m0,R0)
  Planeta.Inicie(r0, 0, 0 ,0.5*V0 , m0, 0.5);
  
  //InicieAnimacion(); //Dibujar
  
  for(t=0, tdibujo=0; t<tmax; t+=dt, tdibujo+=dt){
    
    
    //Dibujar animacion
    if(tdibujo>tcuadro){
      /*  
	  InicieCuadro();
	  Planeta.Dibujese();
	  TermineCuadro();
      */
      
      // hacer un plot
      std::cout<< Planeta.Getx()<<"\t"<<Planeta.Gety()<<endl;
      //tdibujo=0;
    }
    //muevase por OMELYAN PEFRL
    Planeta.Mueva_r(dt,E);
    Planeta.CalculeFuerza();   Planeta.Mueva_V(dt,coeficiente1);
    Planeta.Mueva_r(dt,X);
    Planeta.CalculeFuerza();   Planeta.Mueva_V(dt,L);
    Planeta.Mueva_r(dt,coeficiente2);
    Planeta.CalculeFuerza();   Planeta.Mueva_V(dt,L);
    Planeta.Mueva_r(dt,X);
    Planeta.CalculeFuerza();   Planeta.Mueva_V(dt,coeficiente1);
    Planeta.Mueva_r(dt,E);
  }   
  return 0;
}


