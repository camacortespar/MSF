//OrbitPropagator.h
#pragma once
#include "master.h"


class OrbitPropagator{
 private:
  vector3D r, v, a;
  double Tmax=0, dt=0, masa=1, dmdt=0, max_alt=0, min_alt=0;
  bool coes, deg, stop;
  int tcuadro=1;

 public:
  std::string file="OP";
  
  template <typename T>
  void inicie(std::vector <double> &state0, double tspan, double DT,std::string name,const  T &CB, bool COES, bool DEG, perturbations &perts,double m,int resol, StopC &sc);
  template <typename T>
  void propagate_orbit(const T &CB,perturbations &perts);
  template <typename T>
  void CalculeAceleracion(const T &CB,perturbations &perts);
  void Mueva_r(double t, double coeficiente);
  void Mueva_v(double t, double coeficiente);
  void Mueva_m(double t, double coeficiente);
  void BorreAceleracion(void){a.cargue(0,0,0);};
  double Getx(void){return r.x();};                                                                      
  double Gety(void){return r.y();};
  double Getz(void){return r.z();};
  template <typename T>
  bool check_deorbit(T &CB);
  template <typename T>
  double alt(T &CB);
};

//-------------------------Implementar funciones------------------------------

template <typename T>
void OrbitPropagator::inicie(std::vector <double> &state0, double tspan, double DT,std::string name,const T &CB, bool COES, bool DEG,perturbations &perts,double m,int resol, StopC &sc)
{
  if(perts.isp==0)
    dmdt=0;
  else
    dmdt=perts.thrust/(perts.isp*CB.g0);
  
  max_alt=sc.max_alt; min_alt=sc.min_alt;
  coes=COES;
  deg=DEG;
  masa=m;
  tcuadro=resol;
  stop=sc.no;
  
  if (coes)
    {
      std::vector <double> state(6,0);
      state=coes2rv(state0,deg,CB);
      r.cargue(state[0],state[1],state[2]);
      v.cargue(state[3],state[4],state[5]);
    }
  else{
    r.cargue(state0[0],state0[1],state0[2]);
    v.cargue(state0[3],state0[4],state0[5]);
  }

  Tmax=tspan; dt=DT;
  file=name;
  propagate_orbit(CB,perts);
}
template <typename T>
void OrbitPropagator::CalculeAceleracion(const T &CB,perturbations &perts){
 BorreAceleracion();
 double aux=-CB.mu*std::pow(norma2(r),-1.5);
 a+=aux*r;
 if(perts.J2)
   {
     double z2=std::pow(r.z(),2);
     double r2=norma2(r);
     double norm_r=std::sqrt(r2);
     double tx=r.x()/norm_r*(5*z2/r2-1);
     double ty=r.y()/norm_r*(5*z2/r2-1);
     double tz=r.z()/norm_r*(5*z2/r2-3);
     aux=1.5*CB.J2*CB.mu*std::pow(CB.radius,2)/std::pow(norm_r,4);
     vector3D a_j2; a_j2.cargue(tx,ty,tz);
     a+=aux*a_j2;
   }
 if(perts.aero)
   {
     double r2=norma2(r);
     double norm_r=std::sqrt(r2);
     vector3D v_rel, drag, atmrot; double z, rho, vv;
     z=norm_r-CB.radius;
     rho=calc_atm_density(z,CB);
     atmrot.cargue(CB.atmrotvector[0],CB.atmrotvector[1],CB.atmrotvector[2]);

     v_rel=(v^atmrot); vv=std::sqrt(norma2(v_rel));

     drag=-1*v_rel*0.5*rho*perts.Cd*perts.A*vv/masa;
   }

     if(max_alt>alt(CB) and  min_alt<alt(CB) || stop) //condiciones para detenerse
       {
	 a+=perts.thrust_direction*v*perts.thrust/(masa*1000); //km/s² 
       }
}

void OrbitPropagator::Mueva_r(double t, double coeficiente){
  r+=v*t*coeficiente;
}

void OrbitPropagator::Mueva_v(double t, double coeficiente){
  v+=a*t*coeficiente;
}
void OrbitPropagator::Mueva_m(double t, double coeficiente){
  if(masa>1)
  masa+=dmdt*t*coeficiente;
}
template <typename T>
void OrbitPropagator::propagate_orbit(const T &CB,perturbations &perts)
{
  double E=0.1786178958448091e00;
  double L=-0.2123418310626054e0;
  double X=-0.6626458266981849e-1;
  double coeficiente1=(1-2*L)/2;
  double coeficiente2=(1-2*(X+E));
  int cuenta=0;
 
 std::ofstream fout;   //Salida de la orbita.                                                                                                                   
 fout.open(file+".dat");

 for(double t=0; t<=Tmax; t+=dt, cuenta++)
    {
    Mueva_r(dt,E);
    Mueva_m(dt,E);
    CalculeAceleracion(CB,perts);   Mueva_v(dt,coeficiente1);
    Mueva_r(dt,X);
    Mueva_m(dt,X);
    CalculeAceleracion(CB,perts);   Mueva_v(dt,L);
    Mueva_r(dt,coeficiente2);
    CalculeAceleracion(CB,perts);   Mueva_v(dt,L);
    Mueva_r(dt,X);
    Mueva_m(dt,X);
    CalculeAceleracion(CB,perts);   Mueva_v(dt,coeficiente1);
    Mueva_r(dt,E);
    Mueva_m(dt,E);
    
     if(check_deorbit(CB))
      {
	std::cout <<" after "<<t<<"seconds"<<std::endl;
	break;
	}
    
    if(cuenta%tcuadro==0){
      fout << r.x() <<"\t"<< r.y() <<"\t"<< r.z()<<"\t"<< v.x() <<"\t"<< v.y() <<"\t"<< v.z()<<"\t"<< t << std::endl;
      cuenta=0;}
    }
  fout.close();
}
template <typename T>
bool OrbitPropagator::check_deorbit(T &CB)
{
  double z=alt(CB);
  if (z<CB.deorbit_altitude){
    std::cout<<"Spacecraft "<<file<<" has deorbited";
    return true;
  }
  return false;
}
template <typename T>
double OrbitPropagator::alt(T &CB) //calcula la altitud
{
  double r2=norma2(r);
  double norm_r=std::sqrt(r2);
  double z=norm_r-CB.radius;
  return z;
}
