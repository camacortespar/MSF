#Rebote con oscilacion - solucion:
#Ley de potencias entre el periodo critico y la amplitud de oscilacion de la raqueta (pared inferior).
#
#Autor: Camilo Cortés Parra
import matplotlib.pyplot as plt
import numpy as np

#Listas con las mediciones
A = [2, 5, 10, 20, 50]                          #valores de amplitud
T_c = [3.5625, 5.3125, 7.5, 9.875, 15.5]        #periodos criticos

#Listas con el log10 de las mediciones
log_A = np.log10(A)
log_T = np.log10(T_c)

#Ajuste de potencias
m, b = np.polyfit(log_A,log_T, deg = 1)     #en escala log-log es una recta
n = m; C = 10**(b)
y = C*A**(n)                                #fit
print("Valores del ajuste -> n =",n,"y C =",C)

#Plot datos y fit
f = plt.figure(1)
plt.plot(A, T_c, 'yo', label=r"Datos")
plt.plot(A, y, 'k--', label=r"Ajuste")
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$log(A)$")
plt.ylabel(r"$log(T_{c})$")
plt.legend(loc='best')
f.set_size_inches(8, 6, forward = True)
plt.savefig("TcvsA_loglog.jpg")