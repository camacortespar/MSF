#Código para plotear datos guardados en .txt

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["text.usetex"] = True

#cargar datos en .txt
epidemia = np.loadtxt(
    "/home/camilo/Documentos/Maestria/2022-II/metodos/taller-01/01-ejercicio/p1a.txt"
)

tiempo = epidemia[:, 0]
susceptibles = epidemia[:, 1]
infectados = epidemia[:, 2]
recuperados = epidemia[:, 3]

slim = np.loadtxt(
    "/home/camilo/Documentos/Maestria/2022-II/metodos/taller-01/01-ejercicio/p1c.txt"
)

R0 = slim[:, 0]
sinf = slim[:, 1]

#plot de SIR
f = plt.figure(1)
plt.plot(tiempo, susceptibles, "b", label="Susceptibles")
plt.plot(tiempo, infectados, "r", label="Infectados")
plt.plot(tiempo, recuperados, "g", label="Recuperados")
plt.title("Modelo SIR sin muertes o nacimientos naturales")
plt.xlabel("Tiempo")
plt.ylabel("Susceptibles, Infectados, Recuperados")
plt.legend(loc='upper right')
plt.xlim(xmin=tiempo[0], xmax=tiempo[-1])
plt.ylim(ymin=0, ymax=1)
plt.savefig("sir.jpg")

#plot de S limite
g = plt.figure(2)
plt.plot(R0, sinf, "c", label=r"No. de Susceptibles para t $\rightarrow \infty$")
plt.title(r"$s_{\infty}$ vs $R_{0}$")
plt.legend()
plt.xlabel(r"R$_{0}=\beta / \gamma$", fontsize=15)
plt.ylabel(r"s$_{\infty}$", fontsize=15)
plt.xlim(xmin=1, xmax=R0[-1])
plt.ylim(ymin=0, ymax=sinf[0])
plt.savefig("slim.jpg")
