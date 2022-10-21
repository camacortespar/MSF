//main.cpp

#include "master.h"



int main(void)
{
  //----------------parametros del propagador--------------------
  
  int N=1; //numero de orbitas
  earth cb; //cuerpo central
  
  double r0=cb.radius+414;
  double omega=std::sqrt(cb.mu/(r0*r0*r0)), T=2*M_PI/omega, V0=r0*omega*0.9;
  
  std::vector <OrbitPropagator> OP(N);

  //parametros cuando se tienen r y v
  // std::vector <double> state0{r0,0,500,0,V0,10}; //condiciones iniciales

  //parametros cuando se tienen coes
  std::vector <double> state0{r0,0.0006189,51.6393,0.0,234.1955,105.6372};
  std::vector <double> state1(6,0);

  bool coes=true; //dice si los datos iniaciales están en coes o en velocidad y posición
  bool deg=true; //dice si los datos estan en grados o radianes
  
  double tspan=T*2;
  double dt=1; 

 //----------------parametros de Plot_orbit--------------------
  
  bool save=false; //guarda la gráfica
  
 //--------------------------------------------------------------

  state1= tlecoes("iss.txt",cb);


  
//OP[N](condiciones iniciales, tmax, dt, "donde se guardan los datos", cuerpo central, coes, deg)
  OP[0].inicie(state1,tspan,dt,"OP",cb,coes,deg);

//Plot_orbit(el vector de propagadores,cuerpo central ,"título", save);
  Plot_orbit(OP,cb,"N orbits", save);
  
  return 0;
}
