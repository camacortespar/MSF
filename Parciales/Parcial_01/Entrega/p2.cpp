/*  Rebote con oscilacion - solucion:
    Se simula el rebote de un bola de ping pong cuando la raqueta (pared inferior) oscula.
    
    Autor: Camilo Cortés Parra
*/
#include <iostream>
#include <cmath>
#include "Vector.h"
using namespace std;

//-----Decalaracion de constantes globales-----
const double g = 9.8, K = 1.0e4, Gamma = 10;
const double Lx = 10, Ly = 60;
const int Nx = 1, Ny = 1;
const int N = Nx*Ny;
const double A = 1, T = 5, W = 2*M_PI/T;      //raqueta oscilante

//Constantes de PEFRL
const double zeta=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double coef1=(1.0-2.0*lambda)/2.0;
const double coef2 =1.0-2.0*(chi+zeta);

//-----Decalaracion e implementacion de clases-----
class Cuerpo;
class Colisionador;

//-----Cuerpo-----
class Cuerpo{
private:
  vector3D  r, V, F;   double m, R;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void BorreFuerza(void){F.load(0,0,0);};
  void AdicioneFuerza(vector3D F0){F += F0;};
  void Mueva_r(double dt, double coeficiente);
  void Mueva_V(double dt, double coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();};     //inline
  double Gety(void){return r.y();};     //inline
  void Oscile(double t);                //funciones para la oscilacion de la raqueta (pared inferior)
  void Stop(double t);

  friend class Colisionador;
};

void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0); m = m0;  R = R0; 
}

void Cuerpo::Mueva_r(double dt, double coeficiente){
  r += V*dt*coeficiente; 
}

void Cuerpo::Mueva_V(double dt, double coeficiente){
  V += (F*dt*coeficiente)/m;
}

void Cuerpo::Oscile(double t){
  vector3D oscilar;
  oscilar.load(0,A*sin(W*t),0);
  r += oscilar;
}
void Cuerpo::Stop(double t){
  vector3D oscilar;
  oscilar.load(0,A*sin(W*t),0);
  r -= oscilar;
 }

void Cuerpo::Dibujese(void){ 
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";     //dibujo de la bola de ping pong
}

//-----Colisionador-----
class Colisionador{

private:
public:
  void CalculeFuerzas(Cuerpo * Particula);
  void CalculeFuerzaPared(Cuerpo & Particula1, Cuerpo & Particula2);
};

void Colisionador::CalculeFuerzas(Cuerpo * Particula){
  int i, j; vector3D Fg;

  //Borre todas las fuerzas
  for(i=0; i<N+1; i++){
    Particula[i].BorreFuerza();
  }
 
  //Agregue fuerza de gravedad a la bola de pinp pong
  for(i=0; i<N; i++){
    Fg.load(0,-Particula[i].m*g,0);
    Particula[i].AdicioneFuerza(Fg);
  }
  
  //Calcule la fuerza debida a la raqueta (pared inferior)
  for (i=0; i<N+1; i++){
    for (j=i+1; j<N+1; j++){
      CalculeFuerzaPared(Particula[i], Particula[j]);
    }
  }
}

//Fuerzas debido al contacto con la raqueta (pared inferior)
void Colisionador::CalculeFuerzaPared(Cuerpo & Particula1, Cuerpo & Particula2){
  vector3D r21 = Particula2.r-Particula1.r;
  double d = r21.norm();
  double s = (Particula1.R+Particula2.R)-d;     //distancia interpenetracion
  vector3D n = r21*(1.0/d);
  double F = 0;

  if (s > 0){
    vector3D u; u.load(0,1,0);
    double Vy = Particula1.V*u;                                //componente vertical de la velocidad
    F = K*pow(s,1.5) - Gamma*Particula1.m*Vy*pow(s,0.5);       //fuerza de oscilacion
    if(F > 0){
    vector3D Fr = F*n;
    Particula1.AdicioneFuerza((-1)*Fr); Particula2.AdicioneFuerza(Fr);
    }
  }
}

//-----Animacion-----
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;           //realiza gif 
  cout<<"set output 'rebote_osc.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;        //define rangos
  cout<<"set yrange[-10:"<<Ly<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;                  //proporcion 1:1 entre las unidades de x e y
  cout<<"set parametric"<<endl;                     //se realiza una grafica parametrica
  cout<<"set trange [0:7]"<<endl;                   //rango de la animacion
  cout<<"set isosamples 12"<<endl;                  //realiza 12 divisiones dentro del rango
}

void InicieCuadro(void){
  cout<<"plot 0,0 ";   
  //cout<<" , "<<Lx/7<<"*t,0";             //raqueta(pared inferior) || NO se como hacer oscilar la raqueta
}

void TermineCuadro(void){
  cout<<endl;
}

//-----Main-----
int main(void){

  Cuerpo Particula[N+1];        //bola de ping pong + raqueta (pensada como un grano muuuuuy grande)
  Colisionador PingPong;
  
  double m0 = 1, R0 = 2;
  double t, tdibujo, tmax = 200, tcuadro = tmax/400, dt = 0.0001;
  int ix, iy, i;
  
  //Inicializacion de la raqueta (pared inferior)    
  double Rpared = 10000*Lx;
  double Mpared = 100*m0;
  //-----------------(x0,y0,Vx0,Vy0,m0,R0)
  Particula[N].Inicie(Lx/2,-Rpared, 0, 0, Mpared, Rpared);
  
  //Inicializacion de la bola de ping pong
  //-----------------(x0,y0,Vx0,Vy0,m0,R0)
  Particula[0].Inicie(5,30,0,0,m0,R0);      
  
  InicieAnimacion();
  
  for(t=0, tdibujo=0; t<tmax; t+=dt, tdibujo+=dt){
    //Realizar animacion
    if(tdibujo > tcuadro){
      InicieCuadro();
      for(int i=0; i<N; i++){
	    Particula[i].Dibujese();
      }
      TermineCuadro();      
      tdibujo=0;
    }

    //Datos de altura vs tiempo
    //cout<< t<<"\t"<<Particula[0].Gety()<<endl;

    Particula[N].Oscile(t);         //raqueta empiece a oscilar

    //Movimiento por PEFRL    
    for(i=0; i<N; i++)Particula[i].Mueva_r(dt,zeta);
    PingPong.CalculeFuerzas(Particula);  
    for(i=0; i<N; i++)Particula[i].Mueva_V(dt,coef1);
    for(i=0; i<N; i++)Particula[i].Mueva_r(dt,chi);
    PingPong.CalculeFuerzas(Particula);
    for(i=0; i<N; i++)Particula[i].Mueva_V(dt,lambda);
    for(i=0; i<N; i++)Particula[i].Mueva_r(dt,coef2);
    PingPong.CalculeFuerzas(Particula);
    for(i=0; i<N; i++)Particula[i].Mueva_V(dt,lambda);
    for(i=0; i<N; i++)Particula[i].Mueva_r(dt,chi);
    PingPong.CalculeFuerzas(Particula);
    for(i=0; i<N; i++)Particula[i].Mueva_V(dt,coef1);
    for(i=0; i<N; i++)Particula[i].Mueva_r(dt,zeta);

    Particula[N].Stop(t);       //raqueta deje de oscilar
  }   
  return 0;
}