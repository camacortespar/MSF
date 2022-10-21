#Ejercicio 4.a
#Código para plotear datos guardados en .txt

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["text.usetex"] = True

#Nombre de los archivos a llamar
#filenumbers = ["1", "2", "3", "4", "5", "6", "7"]

#for i in range(8):
#    x = "k" + filenumbers[i] + "_p4a.txt"
#    filenumbers.append(x)
#filenames = filenumbers[7:-1]

#Carga datos en los .txt
data = np.loadtxt("C:\\Users\ASUS\OneDrive - Universidad Nacional de Colombia\Documentos\Física - Maestría\\2022-II\Mecánica Cuántica Avanzada\Taller 03\datos.dat")
data_gamma = np.loadtxt("C:\\Users\ASUS\OneDrive - Universidad Nacional de Colombia\Documentos\Física - Maestría\\2022-II\Mecánica Cuántica Avanzada\Taller 03\datos_gamma.dat")


t = data[:, 0]
normie = data[:, 1]
gamma = data_gamma[:, 1]

#tau = np.zeros((347639, 7))
#for i in range(7):
#    tau[:, i] = np.loadtxt("/home/camilo/Documentos/Maestria/2022-II/metodos/taller-01/04-ejercicio/" + filenames[i], usecols=1)
#    print("El torque máximo para K" + str(i+1) + " es:", "{:.2e}".format(max(tau[:,i])))
    
#Plot
f = plt.figure(1)
plt.plot(t, normie, color = "#FA8500", label=r"$\gamma = 0$")
plt.plot(t, gamma, color = "#FACB1F", label=r"$\gamma = 10$")
#plt.title("Torque del péndulo intermedio en función del tiempo")
plt.xlabel(r"$t$")
plt.ylabel(r"$y$")
plt.legend(loc='upper right')
#plt.xlim(xmin=0.1744, xmax=0.1752)
plt.savefig("eje_1.jpg")