// Simular el movimiento de un balón por Dinámica Molecular
#include <iostream>
#include <cmath>
using namespace std;


const double GM=1.0;

class Cuerpo{
private:
  double x,y,Vx,Vy,Fx,Fy,m,R;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void CalculeFuerza(void);
  void Muevase(double dt);
  void Dibujese(void);
  double Getx(void){return x;}; //inline
  double Gety(void){return y;}; //inline
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
  x=x0; y=y0; Vx=Vx0; Vy=Vy0;  m=m0;  R=R0; 
} 
void Cuerpo::CalculeFuerza(void){
  double aux = -GM*m*std::pow(x*x+y*y,-1.5);
  Fx=aux*x; Fy=aux*y;
} 
void Cuerpo::Muevase(double dt){
  x+=Vx*dt;      y+=Vy*dt;
  Vx+=Fx/m*dt;  Vy+=Fy/m*dt;
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<x<<"+"<<R<<"*cos(t),"<<y<<"+"<<R<<"*sin(t)";
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
  
  double t, tdibujo, tmax=3.3*T, tcuadro=T/40, dt=0.00001;
  
  //------------(x0,y0,Vx0,Vy0, m0,R0)
  Planeta.Inicie(r0, 0, 0 ,0.5*V0 , m0, 0.5);

  InicieAnimacion(); //Dibujar
  
  for(t=0, tdibujo=0; t<tmax; t+=dt, tdibujo+=dt){

    
    //Dibujar animacion
    if(tdibujo>tcuadro){
      
    InicieCuadro();
    Planeta.Dibujese();
    TermineCuadro();
    
    
    // hacer un plot
    //std::cout<< Planeta.Getx()<<"\t"<<Planeta.Gety()<<endl;
    tdibujo=0;
  }
    Planeta.CalculeFuerza();
    Planeta.Muevase(dt);
  }   
  return 0;
}

  
