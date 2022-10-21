#include <stdio.h>
#include "SpiceUsr.h"


#define METAKR "/home/wind/Escritorio/Simulacion/Proyecto-Simulaci-n/Carlos/SPICE_DATA/solar_system_kernel.mk"
#define SCLKID -82
#define STRLEN 50

int main(void)
{

  char calet  [STRLEN];
  char sclkst [STRLEN];
  char utctim [STRLEN]="2020-04-10";
  double et;

 furnsh_c ( METAKR );
 //prompt_c ( "Input UTC Time: ", STRLEN, utctim ); //lee una fecha de la consola
 
 printf ( "Converting UTC Time: %s\n", utctim ); //imprime la fecha que leyo

 str2et_c ( utctim, &et );//convierte las fecha a segundos J2000
 
 printf ( "   ET Seconds Past J2000: %16.3f\n", et ); //imprime los segundos J2000

 etcal_c ( et, STRLEN, calet ); //recive los segundo J2000 (et) y los devuelve a fecha calendario (calet)

 printf ( "   Calendar ET (etcal_c): %s\n", calet ); // imprime el tiempo en calendario (calet)

 timout_c ( et,     "YYYY-MON-DDTHR:MN:SC ::TDB",STRLEN, calet); //le da formato al tiempo calendario
 
 printf ( "Calendar ET (timout_c): %s\n", calet );

 


  return 0;
}
