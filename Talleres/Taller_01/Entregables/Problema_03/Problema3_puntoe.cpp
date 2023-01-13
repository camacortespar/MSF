#include <iostream>
#include <cmath>
#include "planetas_lib.h"
#include <fstream>
/////////////////////////////////////////////////////////////////
// Cuerpos Celestes - metodo de integracion PEFRL              //
//                                                             //
// En este programa se simula y anima el movimiento ente dos   //
// planetas con un observador ajeno al sistema, a este último  //
// se le llama colisionador o interaccionador.                 //
// Los objetos celestes son el Sol y Júpiter a quienes vamos a //
// asignar las masas y radios como razones entre ellos:        //
//                                                             //
// m_sol = 1047 m_jupiter     ,  r_sol = 9.73 r_jupiter        //
// m_jupiter = 1              ,  r_jupiter = 1                 //
//                                                             //
// Tambien se agrega un tercer objeto al cual se le denomina   //
// troyano de Júpiter y va a estar ubicado en el punto de      //
// Lagrange L4 el cuel esta a pi/3 en adelanto a Júpiter. Su   //
// masa y radio se escogen:                                    //
// m_troyano = 0.005 m_jupiter ,  r_troyano = 1_jupiter/10     //
/////////////////////////////////////////////////////////////////
const double G=1.0;
// Para el sistema Sol-Jupiter-Troyano hacemos N=3
const int N=3;

// Coeficientes de PEFRL
const double Zeta   = 0.1786178958448091e00;
const double Lambda = -0.2123418310626054e00;
const double Chi    = -0.6626458266981849e-1;
const double Coeficiente1 =(1-2*Lambda)/2;
const double Coeficiente2 =1-2*(Chi+Zeta);

int main(){
  Cuerpo Planeta[N];
  Colisionador Newton;
  double m1=1047, m2=1, m3=0.005, r=1000;
  double theta=M_PI/3;
  double M=m1+m2, x1=-m2*r/M, x2=m1*r/M, x3=x2*cos(theta),y3=x2*sin(theta);
  double omega = sqrt(G*M/pow(r,3.0)), T=2*M_PI/omega, V1=omega*x1, V2=omega*x2;
  //double Vx3=-V2*cos((M_PI/2)-theta), Vy3=V2*sin((M_PI/2)-theta); //otra forma de verlo
  double Vx3=-V2*sin(theta), Vy3=V2*cos(theta);
  double t, dt=0.01;
  double tmax=1.1*T;
  double tdibujo, tcuadro=T/100;
  int i;
  std::ofstream outfile;

  outfile.open("S-J-T_perturbated_periodos.dat");

  Planeta[0].Inicie(x1, 0,0,  0, V1,0,m1,97.3); // Sol
  Planeta[1].Inicie(x2, 0,0,  0, V2,0,m2,10.0); // Júpiter
  Planeta[2].Inicie(x3,y3,0,Vx3*0.995,Vy3*0.995,0,m3, 1.0); // Planeta Troyano Perturbado


  for(t=0,tdibujo=0; t<tmax; t+=dt, tdibujo+=dt){

    double xrotado_0= Planeta[0].Getx()*cos(omega*t)+Planeta[0].Gety()*sin(omega*t);
    double yrotado_0=-Planeta[0].Getx()*sin(omega*t)+Planeta[0].Gety()*cos(omega*t);

    double xrotado_1= Planeta[1].Getx()*cos(omega*t)+Planeta[1].Gety()*sin(omega*t);
    double yrotado_1=-Planeta[1].Getx()*sin(omega*t)+Planeta[1].Gety()*cos(omega*t);

    double xrotado_2= Planeta[2].Getx()*cos(omega*t)+Planeta[2].Gety()*sin(omega*t);
    double yrotado_2=-Planeta[2].Getx()*sin(omega*t)+Planeta[2].Gety()*cos(omega*t);

    outfile<<t<<" "<<xrotado_2<<std::endl;

    // Haga el movimiento y los calculos por PEFRL
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt, Zeta);            // 1
    Newton.CalculeFuerza(Planeta,N,G);
    for(i=0;i<N;i++) Planeta[i].Mueva_v(dt,Coeficiente1);     // 2
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt, Chi);             // 3
    Newton.CalculeFuerza(Planeta,N,G);
    for(i=0;i<N;i++) Planeta[i].Mueva_v(dt,Lambda);           // 4
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Coeficiente2);     // 5
    Newton.CalculeFuerza(Planeta,N,G);
    for(i=0;i<N;i++) Planeta[i].Mueva_v(dt,Lambda);           // 4
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt, Chi);             // 3
    Newton.CalculeFuerza(Planeta,N,G);
    for(i=0;i<N;i++) Planeta[i].Mueva_v(dt,Coeficiente1);     // 2
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt, Zeta);            // 1
  }
  return 0;
}
