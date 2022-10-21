#include <iostream>
#include <cmath>
#include "../vector.h"
/////////////////////////////////////////////////////////////////
// Libreria para el problema3 taller 1                         //
//                                                             //
// Define las clases cuerpo y colisionador los cuales hacen    //
// los cálculos del movimiento y fuerzas y el observador       //
//                                                             //
/////////////////////////////////////////////////////////////////
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
  double Getx(void){return r.x();}; //Inline
  double Gety(void){return r.y();}; //Inline
  double Getz(void){return r.z();}; //Inline
  friend class Colisionador; //Cualquier especimen de la clase colisionador le puede meter la mano a la clase cuerpo
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
//------------------ Class Colisionador ------------------------//
class Colisionador {
private:
public:
  void CalculeFuerza(Cuerpo * Planeta, int N, const double G);
  void CalculeFuerzaEntre(Cuerpo &  Planeta1, Cuerpo & Planeta2, const double G);
};
void Colisionador::CalculeFuerza(Cuerpo * Planeta, int N, const double G){
  int i,j;
  for(i=0;i<N;i++)
    Planeta[i].BorreFuerza();
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++)
      CalculeFuerzaEntre(Planeta[i],Planeta[j], G);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Planeta1, Cuerpo & Planeta2, const double G){
  vector3D r21,n,F1; double d21,F;
  r21 = Planeta2.r - Planeta1.r; d21=r21.norm(); n=r21/d21;
  F = G*Planeta1.m*Planeta2.m*pow(d21,-2.0);
  F1=F*n;
  Planeta1.SumeFuerza(F1); Planeta2.SumeFuerza(F1*(-1));
}
