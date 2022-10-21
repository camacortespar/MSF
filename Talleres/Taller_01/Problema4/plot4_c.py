#Ejercicio 4.c
#Código para plotear la tabla del inciso b y realizar el ajuste de potencias del inciso c

import matplotlib.pyplot as plt
import numpy as np

#Listas con las mediciones del inciso b
k = [0.1e12, 0.2e12, 0.5e12, 1e12, 2e12, 5e12, 10e12]                    #valores de K
tau_max = [1.84e8, 2.46e8, 3.47e8, 4.58e8, 6.04e8, 8.71e8, 1.15e9]       #torques maximos
t_max = [2.44e-4, 1.87e-4, 1.27e-4, 0.96e-4, 0.73e-4, 0.51e-4, 0.39e-4]  #tiempos de media oscilacion

#Listas con el log10 de las mediciones del inciso b
log_k = np.log10(k)
log_tau = np.log10(tau_max)
log_t = np.log10(t_max)

#Ajuste de potencias
m1, b1 = np.polyfit(log_k,log_tau, deg = 1)  #en escala log-log es una recta
a = m1; A = 10**(b1)
y1 = A*k**(a)                                #fit
print("Para el torque -> a =",a,"y A =",A)
m2, b2 = np.polyfit(log_k,log_t, deg = 1)
b = m2; B = 10**(b2)
y2 = B*k**(b)
print("Para el tiempo -> b =",b,"y B =",B)

#plot datos y fit
f = plt.figure(1)
plt.plot(k, tau_max, 'yo', label=r"$\tau_{max}$")
plt.plot(k, y1, 'k--')
plt.plot(k, t_max, 'bo', label=r"$t_{max}$")
plt.plot(k, y2, 'k--')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$K\,(\frac{g}{cm^{1/2}\cdot s^{2}})$")
plt.legend(loc='best')
f.set_size_inches(6, 6, forward = True)
plt.savefig("plot.jpg")
