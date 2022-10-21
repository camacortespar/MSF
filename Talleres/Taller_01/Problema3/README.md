# Sistema Sol-Júpiter

Se encuentran los archivos específicos de cada punto para mayor entendimiento, debido a ello creamos una librería ```planetas_lib.h``` la cual contiene las clases a utilizar en este problema: "Cuerpo" y "Colisionadaor". Cada uno de estos archivos genera un archivo ```.dat``` el cual contiene los datos generados y llevan el siguiente orden : "Sol(x,y) Jupiter(x,y) Troyano(x,y)", esto se aclara ya que los datos no tienen títulos.

Para correr cada uno de estos archivos y generar sus gráficas se utilizara ```make + comando```. Para mayor rapidez de ejecución del código se establece en los programas un tiempo máximo de 1.1T, si desean cambiarlo a 20 periodos aconsejamos hacer ```t_max=21.0*T``` en el programa requerido.

En la carpeta ```img``` se encuentran las imágenes generadas, estas se actualizarán en su escritorio cada vez que corran el ```make```. La carpeta ```plt``` contiene los comandos que se le pasan a gnuplot para generar las gráficas. Sientanse libres de editar si lo consideran necesario.

## Run


Para correr el primer item (Sistema Sol-Júpiter) deben hacer
```
make SJ
```
Para correr el segundo item (Sistema Sol-Júpiter coordenadas rotadas) deben hacer
```
make SJ_r
```
Para correr el tercer item (Sistema Sol-Júpiter+Troyano en L4 de Júpiter) deben hacer
```
make SJT
```
Para correr el cuarto item (Sistema Sol-Júpiter perturbando las velocidades iniciales del troyano) deben hacer
```
make SJT_p
```

Es aconsejable siempre limpiar al final de correr usando
```
make clean
```
