// Simular el movimiento del atomo de hidrogen por PEFRL 
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "Vector.h" 
using namespace std;


const double E=0.1786178958448091e00;
const double L=-0.2123418310626054e0;
const double X=-0.6626458266981849e-1;

const double coeficiente1=(1-2*L)/2;
const double coeficiente2=(1-2*(X+E))/2;

const double K=100;

class Cuerpo{
private:
  vector3D  r, V, F, El, B, P;   double m=1, R=0.1, e=-1, dt;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,std::vector<double> &campoe, std::vector<double> &campob, double Dt);
  void CalculeFuerza(void);
  void Mueva_r(double dt, double coeficiente);
  void Mueva_V(double dt, double coeficiente);
  void Muevase(void);
  double Getx(void){return r.x();};  //inline
  double Gety(void){return r.y();}; //inline
  double Getz(void){return r.z();}; //inline
  double Getvx(void){return V.x();};  //inline
  double Getvy(void){return V.y();}; //inline
  double Getvz(void){return V.z();}; //inline
  double GetPx(void){return P.x();}; //inline
  double GetPy(void){return P.y();}; //inline
  double GetPz(void){return P.z();}; //inline
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0, std::vector<double> &campoe, std::vector<double> &campob, double Dt){
  r.cargue(x0,y0,0); V.cargue(Vx0,Vy0,Vy0/110); dt=Dt;
  El.cargue(campoe[0],campoe[1],campoe[2]); B.cargue(campob[0],campob[1],campob[2]);
} 
void Cuerpo::CalculeFuerza(void){

  double aux = -K*std::pow(norma2(r),-1.5);
  F=e*((V^B)+El)+aux*r;
  P=V-0.5*(B^r);
}
void Cuerpo::Mueva_r(double dt, double coeficiente){
  r+=V*dt*coeficiente;
}
void Cuerpo::Mueva_V(double dt, double coeficiente){
  V+=(F*dt*coeficiente)/m;
}

void Cuerpo::Muevase(void){
  Mueva_r(dt,E);
  CalculeFuerza();    Mueva_V(dt,coeficiente1);
  Mueva_r(dt,X);
  CalculeFuerza();    Mueva_V(dt,L);
  Mueva_r(dt,coeficiente2);
  CalculeFuerza();    Mueva_V(dt,L);
  Mueva_r(dt,X);
  CalculeFuerza();    Mueva_V(dt,coeficiente1);
  Mueva_r(dt,E);
}


std::vector <double> radio0(int n);

//-----------  Programa Principal ---------------------  
int main(void){

  Cuerpo Atomo;  int n=10;

  double t, tmax=10, tdibujo=0, tcuadro=0.0001, dt=0.00001;
   
//------------colocar el campo electromagnetico---------
  
  std::vector<double>Electrico{10,0,0};
  std::vector<double>Magnetico{0,0,16};
  std::vector<double>a{0,0}; a=radio0(n);

//------------parametros del problema-------------------

   
  Atomo.Inicie(a[0], 0, 0 ,a[1], Electrico, Magnetico,dt);
  
  std::ofstream fout;
  std::ofstream out;
  fout.open("espacio.dat");
  out.open("plano.dat");
   
  for(t=0; t<tmax; t+=dt, tdibujo+=dt)
    {
      if(tdibujo>tcuadro){
	
	// hacer un plot
	out<< Atomo.Getx()<<"\t"<< Atomo.Gety()<<"\t"<< Atomo.Getz()<<"\t"
	   <<Atomo.GetPx()<<"\t"<< Atomo.GetPy()<<"\t"<< Atomo.GetPz()<<endl;
	
	fout<<Atomo.Getx()<<"\t"<<Atomo.Gety()<<"\t"<<Atomo.Getz()<<endl;
	
	tdibujo=0;	
      }
      Atomo.Muevase();
    }
  fout.close();
  out.close();

  return 0;
}
//------------calcular el radio para un n dado-------

std::vector <double> radio0(int n)
{
  std::vector<double> a(2,0);
  double omega=0;
  
  a[0]=n*n/(2*K);
  
  omega=std::sqrt(K/(a[0]*a[0]*a[0]));
  
  a[1]=a[0]*omega;

  return a;
}
