//simular el movimiento de un balon
#include <iostream>
#include <cmath>

const double g=9.8;

class Cuerpo{
private:
  double x, y, Vx, Vy, Fx, Fy, m, R;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
  void CalculeFuerza(void);
  void Muevase(double dt);
  void Dibujese(void);
  double Getx(void){return x;}; //in line
  double Gety(void){return y;}; //in line 
};


void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
  x=x0; y=y0; Vx=Vx0; Vy=Vy0; m=m0; R=R0;
};

void Cuerpo::CalculeFuerza(void){
  Fx=0;
  Fy=-m*g;
}
void Cuerpo::Muevase(double dt){
  x+=Vx*dt;
  y+=Vy*dt;
  Vx+=Fx*dt/m;
  Vy+=Fy*dt/m;
}
void Cuerpo::Dibujese(void){
  std::cout<<" , "<<x<<"+"<<R<<"*cos(t),"<<y<<"+"<<R<<"*sin(t)";
}
//----------------------funciones de animacion--------------------
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output 'Balon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-1:21]"<<endl;
  cout<<"set yrange[-1:5]"<<endl;
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


int main(void)
{
  std::cout.precision(8);
  std::cout.setf(std::ios::scientific);
  
  Cuerpo Balon;
  double t, tmax=2,tdibujo=0.05, dt=0.01;
  
  //--------(x0,y0,Vx0,Vy0,m0,R0)
  Balon.Inicie(0,0,10,8,0.475,0.15);

  InicieAnimacion();

  for(t=0,tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    
    //dibujar el balon
    InicieCuadro();
    Balon.Dibujese();
    TermineCuadro();
    
    // std::cout << t <<"\t"<< Balon.Getx() <<"\t" << Balon.Gety() << std::endl;
    Balon.CalculeFuerza();
    Balon.Muevase(dt);
  }
  
  return 0;
}
