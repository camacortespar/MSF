
Por: Carlos Andrés del Valle

Compilación: el programa compila normal pero se usó una bandera de optimización para que corra más rápido.

g++ -O3 1.cpp


Ejecución: ./a.out | gnuplot.

El programa crea el gif animado. Las funciones para graficar la altura contra el tiempo están comentadas.


Tiempos de rebote:

Para la pelota con Gamma=0 el tiempo de rebote es (12.0771-2.3814)/2=4.84785 u de tiempo

Se escogieron dos rebotes porque el 2 rebote presenta una pequeña irregularidad

Para la pelota con Gamma=10 el tiempo de rebote cada vez se reduce más pero el del primer rebote es 6.804-2.38140=4.4226 u de tiempo.

el dt utilizado es dt=0.0001
