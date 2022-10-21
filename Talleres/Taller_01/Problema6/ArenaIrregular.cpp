#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
using namespace std;

//---- Constantes ---
const double K=1.0e4;
const double Lx=160, Ly=60;
const int N=20, Ns=80, Ntot=N+Ns+3; //Granos que caen, del fondo y totales

const double g=9.8, Gamma=150, Kcundall=500, mu=0.4;

const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//--- Clases -----
class Cuerpo;
class Colisionador;

//---- Clase cuerpo ---
class Cuerpo{
private:
  vector3D r,V,F; double m,R; double theta,omega,tau,I;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,
	      double theta0,double omega0,double m0,double R0);
  void BorreFuerza(){F.load(0,0,0); tau=0;};
  void AdicioneFuerza(vector3D F0,double tau0){F+=F0; tau+=tau0;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,
		    double theta0,double omega0,double m0,double R0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0;
  theta=theta0; omega=omega0; I=2.0/5*m*R*R;
}
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt); theta+=omega*(Coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m); omega+=tau*(Coeficiente*dt/I);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
      <<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";
}

//--- Clase Colisionador ----
class Colisionador{
private:
  //Matriz para guardar el valor del resorte de Kundall entre cada par de granos
  double xCundall[Ntot][Ntot],sold[Ntot][Ntot];
public:
  void Inicie(void);
  void CalculeFuerzas(Cuerpo * Grano,int Nlive,double dt);
  void CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2
			  ,double & x_Cundall,double & s_old,double dt);
};
void Colisionador::Inicie(void){
  int i,j; //j>i
  for(i=0;i<Ntot;i++)
    for(j=0;j<Ntot;j++)
      xCundall[i][j]=sold[i][j]=0;
}
void Colisionador::CalculeFuerzas(Cuerpo * Grano,int Nlive,double dt){
  int i,j;
  //--- Borrar todas las fuerzas ---
  for(i=0;i<Ntot;i++)
    Grano[i].BorreFuerza();
  //--- Añadir fuerza de gravedad ---
  vector3D F0;
  for(i=0;i<Nlive;i++){ //Se añade únicamente a los granos, se omiten las paredes
    F0.load(0,-Grano[i].m*g,0);
    Grano[i].AdicioneFuerza(F0,0);
  }
  //--- Calcular Fuerzas entre pares de granos ---
    for(i=0;i<Nlive;i++)
    for(j=i+1;j<Nlive;j++)
      CalculeFuerzaEntre(Grano[i],Grano[j],xCundall[i][j],sold[i][j],dt);
  //--- Calcular Fuerzas entre los granos y las paredes ---
    for(i=0;i<Nlive;i++)
    for(j=N;j<Ntot;j++)
      CalculeFuerzaEntre(Grano[i],Grano[j],xCundall[i][j],sold[i][j],dt);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2
				      ,double & x_Cundall,double & s_old,double dt){

  //Cantidades generales para saber si hay contacto
  vector3D r21=Grano2.r-Grano1.r; double R1=Grano1.R, R2=Grano2.R;
  double d=r21.norm(), s=R1+R2-d;

  if(s>0){//si hay contacto

    //Variables a calcular
    vector3D F2,F1,tau2,tau1;

    //Vectores unitarios
    vector3D n=r21*(1.0/d),t,k;
    t.load(n.y(),-n.x(),0); //Se gira el verctor n 90° de forma horaria
    k.load(0,0,1); //Vector perpendicular al plano 2D

    //Velocidades relativas
    vector3D V21=Grano2.V-Grano1.V;
    vector3D Rw; Rw.load(0,0,R1*Grano1.omega+R2*Grano2.omega);
    vector3D Vc=V21-(Rw^n); //Velocidad de contacto (diapositivas)
    double Vn=Vc*n, Vt=Vc*t; //Componentes normal y tangencial

    //Fuerza normal (Fuerza de Hertz-Kuwabara-Kono)
    double m12=(Grano1.m*Grano2.m)/(Grano1.m+Grano2.m);
    double Fn=K*pow(s,1.5)-Gamma*m12*sqrt(s)*Vn;
    if(Fn<0) Fn=0; //Esto asegura que la fuerza sea positiva

    //Fuerza tangencial (Fuerza de Cundall)
    x_Cundall+=Vt*dt; double Ft=-Kcundall*x_Cundall;
    double Ftmax=mu*fabs(Fn); //Límite de la fuerza estática
    if(fabs(Ft)>Ftmax) Ft=Ft/fabs(Ft)*Ftmax; //Condición para empezar fuerza dinámica (Ftmax) en dirección Ft

    //Calcula y Cargue las fuerzas
    F2=n*Fn+t*Ft; tau2=((n*(-R2))^F2); F1=F2*(-1); tau1=((n*R1)^F1);
    Grano2.AdicioneFuerza(F2,tau2*k);   Grano1.AdicioneFuerza(F1,tau1*k);
  }

  //Cuando se cumple esta condición es porque no hay más contacto y el resorte vuelve al eq.
  if(s_old>=0 && s<0) x_Cundall=0;
  s_old=s;
}

//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl;
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
  Cuerpo Grano[Ntot];
  Colisionador Hertz;
  Crandom ran64(1);
  double m0=1, R0=0.0;//, kT=10, V0=sqrt(2*kT/m0);
  int i,Nlive;
  double cuadros=5,t,tdibujo,dt=1e-3,
        tmax=cuadros*sqrt(Ly/g),tcuadro=tmax/(10*cuadros);
  double Omega0, OmegaMax=8.0;
  double Rpared=100*Lx, Mpared=100*m0, Rs=Lx/(Ns*2);

  InicieAnimacion(); //Dibujar

  //Inicializar las paredes
  for(i=0;i<Ns;i++)
    Grano[N+i].Inicie(Rs*(i*2+1),  0,  0,  0,     0,     0,Mpared,Rs); //Granos de suelo
  //------------------(  x0,       y0,Vx0,Vy0,theta0,omega0,m0,R0)
  Grano[N+Ns].Inicie(Lx/2,Ly+Rpared,  0,  0,     0,     0,Mpared,Rpared); //Pared de arriba
  Grano[N+Ns+1].Inicie(Lx+Rpared,Ly/2,  0,  0,     0,     0,Mpared,Rpared); //Pared derecha
  Grano[N+Ns+2].Inicie(  -Rpared,Ly/2,  0,  0,     0,     0,Mpared,Rpared); //Pared izquierda
  //Inicializar los granos
  for(i=0;i<N;i++){
      Omega0=OmegaMax*(2*ran64.r()-1); //Modificamos resultado del generador para dar un # entre -1 y 1
      //-------------(  x0,     y0,Vx0,Vy0, theta0,omega0,m0,R0)
      R0=1.6+ran64.r();
      Grano[i].Inicie(Lx/2,Ly-2*R0,  0,  0,      0,Omega0,m0,R0);
    }
  //Soltamos cada grano uno por uno
  for(Nlive=1;Nlive<=N;Nlive++)
    for(t=0,tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    //Dibujar
    if(tdibujo>tcuadro){
      InicieCuadro();
      for(i=N;i<N+Ns;i++)  Grano[i].Dibujese(); //Dibujar suelo
      for(i=0;i<=Nlive;i++) Grano[i].Dibujese(); //Dibujar granos
      TermineCuadro();

      tdibujo=0;
    }

    //--- Muevase por PEFRL ---
    for(i=0;i<Nlive;i++)Grano[i].Mueva_r(dt,epsilon);
    Hertz.CalculeFuerzas(Grano,Nlive,dt);
    for(i=0;i<Nlive;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<Nlive;i++)Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Grano,Nlive,dt);
    for(i=0;i<Nlive;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<Nlive;i++)Grano[i].Mueva_r(dt,chiepsilon);
    Hertz.CalculeFuerzas(Grano,Nlive,dt);
    for(i=0;i<Nlive;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<Nlive;i++)Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Grano,Nlive,dt);
    for(i=0;i<Nlive;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<Nlive;i++)Grano[i].Mueva_r(dt,epsilon);

  }
  return 0;
}
