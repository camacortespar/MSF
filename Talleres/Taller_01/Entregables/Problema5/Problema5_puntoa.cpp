
//SIMULACIÓN DE UNA PARTÍCULA BIDIMENSIONAL BAJO EL INFLUJO DE UNA FUERZA DE LENNARD-JONES, IMPLEMENTADA COMO UNA FUERZA CENTRAL

#include <iostream>
#include <cmath>
using namespace std;

//Constantes globales
const double E=1.0, KT=0.5;

//Declaración de las clases
class Cuerpo;


//---------- Clase Cuerpo --------------

class Cuerpo{
private:
  double x,y,Vx,Vy,Fx,Fy,m,R;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void CalculeFuerza(double r0);
  void Muevase(double dt);
  double Getx(void){return x;}; //Inline
  double Gety(void){return y;}; //Inline
};

void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
  x=x0; y=y0; Vx=Vx0; Vy=Vy0; m=m0; R=R0;
}


//Implementación de la fuerza de Lennard Jones

void Cuerpo::CalculeFuerza(double r0){
  double aux=(-12*E/(x*x+y*y))*(pow(r0/sqrt(x*x+y*y),12)-pow(r0/sqrt(x*x+y*y),6));
  Fx=-aux*x;  Fy=-aux*y;
}


void Cuerpo::Muevase(double dt){
  x+=Vx*dt;     y+=Vy*dt;
  Vx+=Fx/m*dt;  Vy+=Fy/m*dt;
}

//----------- Funciones Globales -----------
int main(){
  Cuerpo Particula;
  double r0=10,m0=1,V0;
  double t,dt=0.0001;
  V0=sqrt(2*KT/m0);
  
  //------------(x0, y0, Vx0, Vy0, m0, R0)
  Particula.Inicie(r0,  0,   V0,   0, m0,2.5);
  
  for(t=0;t<100;t+=dt){
    cout<<t<<" "<<Particula.Getx()<<endl;
    Particula.CalculeFuerza(r0);
    Particula.Muevase(dt);
  }
  
  return 0;
}
