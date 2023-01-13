#Ejercicio 4.d
#Reescalamiento del plot del inciso a

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["text.usetex"] = True

#Constantes globales
k = [0.1e12, 0.2e12, 0.5e12, 1e12, 2e12, 5e12, 10e12]  #valores de K
t0 = 0.174566                                          #instante del primer contacto (medido con xmgrace)
a = 0.39645076907099885; A = 8040.116032851302         #resultados del fit del inciso c
b = -0.3999314440182036; B = 6.108426217993056

#Nombre de los archivos a llamar
filenumbers = ["1", "2", "3", "4", "5", "6", "7"]

for i in range(8):
    x = "k" + filenumbers[i] + "_p4a.txt"
    filenumbers.append(x)
filenames = filenumbers[7:-1]

#Carga de datos en los .txt
t = np.loadtxt("/home/camilo/Documentos/Maestria/2022-II/metodos/taller-01/04-ejercicio/" + filenames[0], usecols=0)

tau = np.zeros((347639, 7))
for i in range(7):
    tau[:, i] = np.loadtxt("/home/camilo/Documentos/Maestria/2022-II/metodos/taller-01/04-ejercicio/" + filenames[i], usecols=1)
    
#Plot reescalado
f = plt.figure(1)
plt.plot((t-t0)*k[0]**-b, tau[:, 0]*k[0]**-a, color = "#FA8500", label="K = 0.1e12")
plt.plot((t-t0)*k[1]**-b, tau[:, 1]*k[1]**-a, color = "#FACB1F", label="K = 0.2e12")
plt.plot((t-t0)*k[2]**-b, tau[:, 2]*k[2]**-a, color = "#6BFA48", label="K = 0.3e12")
plt.plot((t-t0)*k[3]**-b, tau[:, 3]*k[3]**-a, color = "#63D5FA", label="K = 0.5e12")
plt.plot((t-t0)*k[4]**-b, tau[:, 4]*k[4]**-a, color = "#5B41FA", label="K = 2e12")
plt.plot((t-t0)*k[5]**-b, tau[:, 5]*k[5]**-a, color = "#F042FA", label="K = 5e12")
plt.plot((t-t0)*k[6]**-b, tau[:, 6]*k[6]**-a, color = "#FA2E25", label="K = 10e12")
plt.title("Torque del péndulo intermedio reescalado")
plt.xlabel(r"$(t-t_{0})\cdot K^{-b}$")
plt.ylabel(r"$\tau\cdot K^{-a}$")
plt.legend(loc='upper right')
plt.xlim(xmin=-10, xmax=20)
plt.savefig("torque_reescalado.jpg")

