// Simular el movimiento de 2 planetas por PEFRL
#include <iostream>
#include <cmath>
#include "Vector.h" 
using namespace std;


//------------------------------Declarar constantes------------------

const int N=3;  //numero de cuerpos
const double G=1.0;

const double E=0.1786178958448091e00;
const double L=-0.2123418310626054e0;
const double X=-0.6626458266981849e-1;

const double coeficiente1=(1-2*L)/2;
const double coeficiente2=(1-2*(X+E));

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

  void CalculeFuerzas(Cuerpo * Planeta);
  void CalculeFuerzaEntre(Cuerpo & Planeta1, Cuerpo & Planeta2);
};

void Colisionador::CalculeFuerzas(Cuerpo * Planeta){

  int i, j;
  //borrar todas las fuerzas 
  for(i=0; i<N; i++){
  Planeta[i].BorreFuerza();
  }
  //Calcular todas las fuerzas 
  for (i=0; i<N; i++){
    for (j=i+1; j<N; j++){
      CalculeFuerzaEntre(Planeta[i], Planeta[j]);
    }
  }
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Planeta1, Cuerpo & Planeta2){
  vector3D r21=Planeta2.r-Planeta1.r;
  double aux=G*Planeta2.m*Planeta1.m*std::pow(norma2(r21),-1.5);
  vector3D F1= r21*aux;
  Planeta1.AdicioneFuerza(F1); Planeta2.AdicioneFuerza(-1*F1); 
}

//-------------------------- Funciones de Animacion -------------------

void InicieAnimacion(void){
  //cout<<"set terminal gif animate"<<endl; 
  // cout<<"set output 'Planeta.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-120:120]"<<endl;
  cout<<"set yrange[-120:120]"<<endl;
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

//-----------------------  Programa Principal ------------------------

int main(void){
  Cuerpo Planeta[N];
  Colisionador Newton;
  int i;
  double m0=1047, m1=1, m2=0.005, r=1000;
  double M=m0+m1, Mu=m0*m1/M;
  double x1=r-r*m1/M, x0=-r*m1/M;
  double omega=std::sqrt(G*M*std::pow(r,-3.0));
  double V1=omega*x1,V0=omega*x0, T=2*M_PI/omega;
  double R0=10, R1=5;
  double theta0=M_PI/3;
  double x2=x1*cos(theta0),y2=x1*sin(theta0),vx2=-V1*sin(theta0),vy2=V1*cos(theta0);
  double xrotado0, yrotado0, xrotado1, yrotado1,xrotado2, yrotado2;
  
  double t, tdibujo, tmax=20*T, tcuadro=T/1000, dt=0.0001;

  
  //---------------(x0, y0,Vx0,   Vy0, m0,R0)
  Planeta[0].Inicie(x0, 0, 0 ,   V0 , m0, R0);
  Planeta[1].Inicie(x1, 0, 0 ,   V1 , m1, R1);
  Planeta[2].Inicie(x2, y2, vx2*0.9999 ,vy2*0.9999 , m2,R0);
  
  //InicieAnimacion(); //Dibujar
  
  for(t=0, tdibujo=0; t<tmax; t+=dt, tdibujo+=dt){
    
    xrotado2=Planeta[2].Getx()*std::cos(omega*t)+Planeta[2].Gety()*std::sin(omega*t);
    yrotado2=-Planeta[2].Getx()*std::sin(omega*t)+Planeta[2].Gety()*std::cos(omega*t);
    
    // hacer un plot
    //std::cout<< Planeta[0].Getx()<<"\t"<<Planeta[0].Gety()<<"\t"<< Planeta[1].Getx()<<"\t"<<Planeta[1].Gety()<<"\t"<<Planeta[2].Getx()<<"\t"<<Planeta[2].Gety()<<endl;
    std::cout<< xrotado2 <<"\t"<< yrotado2 <<"\t" << t << endl;

    //Dibujar animacion        
    /*
      if(tdibujo>tcuadro){
      InicieCuadro();
      for(int i=0; i<N; i++)
      Planeta[i].Dibujese();
      TermineCuadro();
      //tdibujo=0;}
      */
    
    
    //muevase por OMELYAN PEFRL
    
    for(i=0; i<N; i++)Planeta[i].Mueva_r(dt,E);
    Newton.CalculeFuerzas(Planeta);  for(i=0; i<N; i++)Planeta[i].Mueva_V(dt,coeficiente1);
    for(i=0; i<N; i++)Planeta[i].Mueva_r(dt,X);
    Newton.CalculeFuerzas(Planeta);  for(i=0; i<N; i++)Planeta[i].Mueva_V(dt,L);
    for(i=0; i<N; i++)Planeta[i].Mueva_r(dt,coeficiente2);
    Newton.CalculeFuerzas(Planeta);  for(i=0; i<N; i++)Planeta[i].Mueva_V(dt,L);
    for(i=0; i<N; i++)Planeta[i].Mueva_r(dt,X);
    Newton.CalculeFuerzas(Planeta);  for(i=0; i<N; i++)Planeta[i].Mueva_V(dt,coeficiente1);
    for(i=0; i<N; i++)Planeta[i].Mueva_r(dt,E);
    
  }   
  return 0;
}


