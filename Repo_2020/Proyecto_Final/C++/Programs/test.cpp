#include "master.h"






int main(void)
{
 furnsh_c("/home/wind/Escritorio/Simulacion/Proyecto-Simulaci-n/Carlos/SPICE_DATA/solar_system_kernel.mk"); //load metakernel

  std::string FRAME="ECLIPTJ2000";
  std::string OBSERVER="SUN";

  //get_objects("SPICE_DATA/de432s.bsp",true);
  
  //----------------parametros del propagador--------------------
  StopC sc;                    //stop conditions
  int N=1;                    //numero de orbitas
  earth cb;                  //cuerpo central
  perturbations perts;      //diccionario de perturbaciones
  double tspan=3600*24*2;  //tmax en segundos
  double dt=1;           // paso de tiempo en segundos
  bool coes=true;        //si condiciones iniciales son coes
  bool deg=true;        //si están en grados
  double masa=100;     //kg
  int tcuadro=10;     //cada cuanto imprime datos


  spkobj_c("de432.spk");





  return 0;
}
