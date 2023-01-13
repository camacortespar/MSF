import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["text.usetex"] = True

filenumbers = ["1", "2", "3", "4", "5", "6", "7"]

for i in range(8):
    x = "datosP4a_K" + filenumbers[i] + ".txt"
    filenumbers.append(x)
filenames = filenumbers[7:-1]

TodoslosK = np.zeros((695277, 7))
tiempo = np.loadtxt(
    "d:/Esteban/Uni Esteban/Simulaciones/RepositorioGrupoSimulaciones/Taller1/Problema4/DatosP1aV2/"
    + filenames[0],
    usecols=0,
)

for i in range(7):
    TodoslosK[:, i] = np.loadtxt(
        "d:/Esteban/Uni Esteban/Simulaciones/RepositorioGrupoSimulaciones/Taller1/Problema4/DatosP1aV2/"
        + filenames[i],
        usecols=1,
    )

print(TodoslosK[-1, :])
print(filenames)

plt.plot(tiempo, TodoslosK)
plt.xlim(xmin=0.174, xmax=0.178)
plt.savefig(
    "d:/Esteban/Uni Esteban/Simulaciones/RepositorioGrupoSimulaciones/Taller1/Problema4/TorqueVStiempo_mejorado.png"
)
