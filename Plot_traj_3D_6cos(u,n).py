import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Abre el archivo
with open('Tr(3D_6C(ura,neptu)_1dia).txt') as fichero:
    # Lee todas las líneas del archivo
    lineas = fichero.readlines()

# Inicializa una lista de listas para almacenar las columnas
A = [[] for _ in range(len(lineas[0].split()))]

# Procesa cada línea y divide en columnas
for linea in lineas:
    valores = [float(valor) for valor in linea.split()]
    
    

 # Almacena cada valor en la lista correspondiente
    for i, valor in enumerate(valores):
        A[i].append(valor)
"""A partir d'aqui tenim separades les dades de les trajectories per columnes, per exemple A[0] és el temps i A[2] la ysol"""



fig= plt.figure()
ax= Axes3D(fig)


ax.plot(A[1],A[2],A[3],"-", color="orange", lw=2, label="Sol")
ax.plot(A[7],A[8],A[9],"-", color="forestgreen",lw=1, label="Terra")
ax.plot(A[13],A[14],A[15], "-",color="red",lw=1, label="Mart")
ax.plot(A[19],A[20],A[21], "-",color="darkgoldenrod",lw=1, label="Jupiter")
ax.plot(A[25],A[26],A[27], "-",color="darkcyan",lw=1, label="Urà")
ax.plot(A[31],A[32],A[33],"-", color="darkorchid", lw=1,label="Neptú")
plt.legend(loc='best', facecolor='w', fontsize=12)

ax.set_xlabel('x[1]')
ax.set_ylabel('y[1]')
ax.set_zlabel('z[1]')

plt.show()



