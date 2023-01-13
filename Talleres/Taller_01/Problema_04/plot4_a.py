#Ejercicio 4.a
#Código para plotear datos guardados en .txt

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["text.usetex"] = True

#Nombre de los archivos a llamar
filenumbers = ["1", "2", "3", "4", "5", "6", "7"]

for i in range(8):
    x = "k" + filenumbers[i] + "_p4a.txt"
    filenumbers.append(x)
filenames = filenumbers[7:-1]

#Carga datos en los .txt
t = np.loadtxt("/home/camilo/Documentos/Maestria/2022-II/metodos/taller-01/04-ejercicio/" + filenames[0], usecols=0)

tau = np.zeros((347639, 7))
for i in range(7):
    tau[:, i] = np.loadtxt("/home/camilo/Documentos/Maestria/2022-II/metodos/taller-01/04-ejercicio/" + filenames[i], usecols=1)
    print("El torque máximo para K" + str(i+1) + " es:", "{:.2e}".format(max(tau[:,i])))
    
#Plot
f = plt.figure(1)
plt.plot(t, tau[:, 0], color = "#FA8500", label="K = 0.1e12")
plt.plot(t, tau[:, 1], color = "#FACB1F", label="K = 0.2e12")
plt.plot(t, tau[:, 2], color = "#6BFA48", label="K = 0.3e12")
plt.plot(t, tau[:, 3], color = "#63D5FA", label="K = 0.5e12")
plt.plot(t, tau[:, 4], color = "#5B41FA", label="K = 2e12")
plt.plot(t, tau[:, 5], color = "#F042FA", label="K = 5e12")
plt.plot(t, tau[:, 6], color = "#FA2E25", label="K = 10e12")
plt.title("Torque del péndulo intermedio en función del tiempo")
plt.xlabel(r"$t\,(s)$")
plt.ylabel(r"$\tau\,(\frac{g\cdot cm^{2}}{s^{2}})$")
plt.legend(loc='upper right')
plt.xlim(xmin=0.1744, xmax=0.1752)
plt.savefig("torque.jpg")

