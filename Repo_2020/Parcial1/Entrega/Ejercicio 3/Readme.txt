
Por: Carlos Andrés del Valle

Compilación: el programa compila normal pero se usó una bandera de optimización para que corra más rápido.

g++ -O3 1.cpp


Ejecución: ./a.out > datos.dat

El programa imprime a la pantalla los datos de y contra t para hacer una gráfica.

el dt utilizado es dt=0.0001

Para hayar el Tc se hizo una busqueda manual con el siguiente algoritmo.

Suponga que el Tc se sabe que esta entre 0 y 1. en 0 es caotico y en 1 no. Así que se revisa en 0.5. Si no es caótico entonces se revisa 0.25 si sí, se revisa 0.75. y así sucesivamente.

Como la búsqueda fue manual, aprovechando el corto tiempo de ejecución del programa, se hizo un Makefile para hacer más rápido el proceso de generar la gráfica.

Se hizo la regresión en Exel y el n encontrado fue de

n=0.43
