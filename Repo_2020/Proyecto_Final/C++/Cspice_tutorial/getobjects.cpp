 #include <stdio.h>
#include <iostream>
#include "SpiceUsr.h"

#define  FILSIZ         256
#define  MAXIV          1000
#define  WINSIZ         ( 2 * MAXIV )
#define  TIMLEN         51
#define  MAXOBJ         1000

int main()
{

  SPICEDOUBLE_CELL        ( cover, WINSIZ );
  SPICEINT_CELL           ( ids,   MAXOBJ );

  SpiceChar               lsk     [ FILSIZ ];
  SpiceChar               spk     [ FILSIZ ]="/home/wind/Escritorio/Simulacion/Proyecto-Simulaci-n/Carlos/SPICE_DATA/de432s.bsp";
  SpiceChar               timstr  [ TIMLEN ];

  SpiceDouble             b;
  SpiceDouble             e;

  SpiceInt                i;
  SpiceInt                j;
  SpiceInt                niv;
  SpiceInt                obj;

  furnsh_c ("/home/wind/Escritorio/Simulacion/Proyecto-Simulaci-n/Carlos/SPICE_DATA/solar_system_kernel.mk" ); //metakernel

  spkobj_c ( spk, &ids ); //spk file

  /*
    We want to display the coverage for each object. Loop over
    the contents of the ID code set, find the coverage for
    each item in the set, and display the coverage.
  */
  for ( i = 0;  i < card_c( &ids );  i++  )
    {
      /*
	Find the coverage window for the current object.
	Empty the coverage window each time so we don't
	include data for the previous object.
      */                
      obj  =  SPICE_CELL_ELEM_I( &ids, i );

      scard_c  ( 0,        &cover );
      spkcov_c ( spk, obj, &cover );
      
      /*
	Get the number of intervals in the coverage window.
      */
      niv = wncard_c ( &cover );
	  
      /*
	Display a simple banner.
      */
      printf ( "%s\n", "========================================" );

      char body[1000];
      SpiceBoolean *found ;
      bodc2n_c((int)obj,10000000,body,found );

      printf ( "Coverage for object");
      std::cout<<body<<std::endl;

      /*
	Convert the coverage interval start and stop times to TDB
	calendar strings.
      */
      for ( j = 0;  j < niv;  j++  )
	{
	  /*
	  Get the endpoints of the jth interval.
	  */
	  wnfetd_c ( &cover, j, &b, &e );
	  
	  /*
	    Convert the endpoints to TDB calendar
	    format time strings and display them.
	  */

	  timout_c ( b,"YYYY MON DD HR:MN:SC.### (TDB) ::TDB",TIMLEN,timstr);
	    
	  printf ( "\n" "Interval:  %d\n" "Start: %s\n",j, timstr);
	  
	  timout_c ( e,"YYYY MON DD HR:MN:SC.### (TDB) ::TDB",TIMLEN,timstr);
	  printf ( "Stop: %s\n", timstr );
	}
    }
  return  0;
}
