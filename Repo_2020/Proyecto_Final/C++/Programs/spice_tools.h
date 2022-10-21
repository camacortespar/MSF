#include "master.h"



void get_objects(std::string c,bool display);


void get_objects(std::string c, bool display)
{
  // char a=c.c_str();
#define  FILSIZ 256
#define  MAXIV 1000
#define  WINSIZ ( 2 * MAXIV )
#define  TIMLEN 51
#define  MAXOBJ 1000
  
  SPICEDOUBLE_CELL ( cover, WINSIZ );
  SPICEINT_CELL ( ids,   MAXOBJ );
  
  SpiceChar lsk     [ FILSIZ ];
  SpiceChar spk     [ FILSIZ ];
  SpiceChar timstr  [ TIMLEN ];
  
  SpiceDouble b;
  SpiceDouble e;
  
  SpiceInt i;
  SpiceInt j;
  SpiceInt niv;
  SpiceInt obj;
  
  spkobj_c ("/home/wind/Escritorio/Simulacion/Proyecto-Simulaci-n/Carlos/SPICE_DATA/de432s.bsp", &ids);
  for ( i = 0;  i < card_c( &ids );  i++  )
    {
      obj  =  SPICE_CELL_ELEM_I( &ids, i );
      scard_c  ( 0,        &cover );
      spkcov_c ( spk, obj, &cover );
      niv = wncard_c ( &cover );
      
      printf ( "%s\n", "========================================" );
      
      printf ( "Coverage for object %d\n", (int)obj );
      
      for ( j = 0;  j < niv;  j++  )
	{
	  wnfetd_c ( &cover, j, &b, &e );
	  timout_c ( b,
		     "YYYY MON DD HR:MN:SC.### (TDB) ::TDB",
		     TIMLEN,
		     timstr                                  );
	  
	  printf ( "\n"
		   "Interval:  %d\n"
		   "Start:     %s\n",
		   j,
		   timstr            );
	  
	  timout_c ( e,
		     "YYYY MON DD HR:MN:SC.### (TDB) ::TDB",
		     TIMLEN,
		     timstr                                  );
	  printf ( "Stop:      %s\n", timstr );
	  
	}
      
    }
}
		



    
    
