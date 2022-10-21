#include <stdio.h>
#include "SpiceUsr.h"


#define METAKR "/home/wind/Escritorio/Simulacion/Proyecto-Simulaci-n/Carlos/SPICE_DATA/solar_system_kernel.mk"
#define STRLEN 50

int main(void)
{
  char calet  [STRLEN];
  char sclkst [STRLEN];
  char utctim [STRLEN]="2020-04-10";
  double et;
  SpiceDouble dist;
  SpiceDouble ltime;
  SpiceDouble pos   [3];
  SpiceDouble state [6];
  

 furnsh_c ( METAKR );
 //prompt_c ( "Input UTC Time: ", STRLEN, utctim ); //lee una fecha de la consola
 
 printf ( "Converting UTC Time: %s\n", utctim ); //imprime la fecha que leyo

 str2et_c ( utctim, &et );//convierte las fecha a segundos J2000
 
 printf ( "   ET Seconds Past J2000: %16.3f\n", et ); //imprime los segundos J2000

 spkezr_c ( "JUPITER BARYCENTER", et, "J2000", "LT+S","EARTH BARYCENTER", state, &ltime); //los datos de jupiter vistos desde la tierra

 printf ( "   Apparent state of Jupiter as seen "
               "from Earth in the J2000\n"        );
      printf ( "      frame (km, km/s):\n"          );
      printf ( "      X = %16.3f\n", state[0]       );
      printf ( "      Y = %16.3f\n", state[1]       );
      printf ( "      Z = %16.3f\n", state[2]       );
      printf ( "     VX = %16.3f\n", state[3]       );
      printf ( "     VY = %16.3f\n", state[4]       );
      printf ( "     VZ = %16.3f\n", state[5]       );

spkpos_c ( "EARTH", et,  "J2000", "LT+S","NEPTUNE BARYCENTER",   pos, &ltime ); //hace lo mismo que la funcion pasada pero no devuelve velocidad
 
      printf ( "   Apparent position of Earth as "
               "seen from NEPTUNE in the J2000\n"     );
      printf ( "      frame (km): \n"                 );
      printf ( "      X = %16.3f\n", pos[0]           );
      printf ( "      Y = %16.3f\n", pos[1]           );
      printf ( "      Z = %16.3f\n", pos[2]           );



//distancia entre la tierra y el sol
 spkpos_c ( "SUN",  et,  "J2000", "NONE","EARTH BARYCENTER", pos, &ltime);
 
      /*
      Compute the distance between the body centers in
      kilometers.
      */
      dist = vnorm_c ( pos );
 
      /*
      Convert this value to AU using convrt_c.
      */
      convrt_c ( dist, "KM", "AU", &dist );
 
      printf ( "   Actual distance between Sun and "
               "earth body centers:\n"               );
      printf ( "      (AU): %16.3f\n", dist           );

 return 0;
}
