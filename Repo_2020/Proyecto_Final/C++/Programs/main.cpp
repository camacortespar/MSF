//main.cpp

#include "master.h"


int main(void)
{
  //----------------------------------SPICE parameters-----------
  
  furnsh_c("/home/wind/Escritorio/Simulacion/Proyecto-Simulaci-n/Carlos/SPICE_DATA/solar_system_kernel.mk"); //load metakernel

  std::string FRAME="ECLIPTJ2000";
  std::string OBSERVER="SUN";

  //get_objects("SPICE_DATA/de432s.bsp",true);
  
  //----------------parametros del propagador--------------------
  StopC sc;                    //stop conditions
  int N=1;                    //numero de orbitas
  earth cb;                  //cuerpo central
  perturbations perts;      //diccionario de perturbaciones
  double tspan=1;  //tmax en segundos
  double dt=0.1;           // paso de tiempo en segundos
  bool coes=true;        //si condiciones iniciales son coes
  bool deg=true;        //si están en grados
  double masa=100;     //kg
  int tcuadro=10;     //cada cuanto imprime datos
 
  //------------------definir las perturbaciones----------------
  
  //perts.J2=true;
  // perts.aero=true; perts.Cd=2.2; perts.A=std::pow(0.001,2); //Cd=coeficiente de fricción A=area en km²
  // perts.thrust=0.327; perts.isp=4300; perts.thrust_direction=1;
  
  //------------------definir el propagador---------------------
  
  std::vector <OrbitPropagator> OP(N);

  //-----------------definir condiciones iniciales-------------
  
  // si coes=false state0{rx,ry,rx,vx,vy,vz}
  //si coes=true   state0{a,e,i,ta,aop,raan}

  double r0=cb.radius+1; //{r0,0.0,20,0.0,234.1955,105.6372};
  std::vector <double> state0=tlecoes("iss.txt",cb);
  // std::vector <double>state1= tlecoes("HJ-2A.txt",cb); //sacar condiciones iniciales de los TLE
  // std::vector <double>state2= tlecoes("COSMOS2251.txt",cb);

 //------------------definir las stop conditions----------------

 sc.max_alt=10000; sc.min_alt=-10; sc.no=true; //mas condiciones es necesario implementarlas, no=true es para desactivarlas
  
 //------------------parametros de Plot_orbit--------------------
  
  bool save_orbit=false; //guarda la gráfica
  bool show_orbit=true; //muestra la grafica al final de correr
  std::vector<std::string> labels {"ISS"}; //labels de los propagadores
  std::string title_orbit ="Órbita del satélite HJ-2A"; //titulo de la gráfica

  //----------------parametros de Plot_coes---------------------

  bool show_coes=true;
  bool save_coes=false;
  bool hours=true;        //define las unidades de tiempo. Si las dos son false usa segundos
  bool days=false;      //tener cuidado de que las dos no esten en true
  
  
 //---------------------- Propagar orbitas ----------------------

  //si hay varios propagadores el nombre "OP" debe ser distinto
  //es el archivo .dat donde se guardan los datos
  
  OP[0].inicie(state0,tspan,dt,"OP",cb,coes,deg,perts,masa,tcuadro,sc);
  // OP[1].inicie(state1,tspan,dt,"OP1",cb,coes,deg,perts,masa,tcuadro,sc);
  //OP[2].inicie(state2,tspan,dt,"OP2",cb,coes,deg,perts,masa,tcuadro,sc);

  //-----------------------pintar las orbitas--------------------

  
  //Plot_orbit_gnuplot(el vector de propagadores,cuerpo central ,"título", save);
  //Plot_orbit_gnuplot(OP,cb,"N orbits", save_orbit); //grafica bonita pero sin ejes ni labels
  Plot_orbit(OP,cb,show_orbit,save_orbit,title_orbit,labels); //grafica en python con labeles y ejes (lento)
  
 Plot_coes(OP,cb,show_coes,save_coes,labels,hours,days);

  
  return 0;
}
