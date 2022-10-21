#include <iostream>
#include <cmath>
#include "Vector.h"       //Libreria de los vectores



int main (void)
{
  double d=0;
  vector3D a, b, c;        //Declarar vector
  a.cargue(1,2,3);        //Colocar los valores
  a.show();              //Imprimir vector

  b.cargue(1,0,1);

  c=3*a+b;               //Suma de vectores y multiplicación por escalar
  c=a-b;                //Resta de vectores
  d=a*b;               //Producto punto de vectores
  c=a^b;              //Producto cruz de vecctores
  c=a+(a^b);         //El producto cruz tiene menos prioridad que la suam


  c.show();

  return 0;
}
