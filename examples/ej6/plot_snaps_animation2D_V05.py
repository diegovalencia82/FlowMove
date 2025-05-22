import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def read_snapshot(filepath):
    # Leer datos del archivo snapshot
    with open(filepath, 'r') as file:
        lines = file.readlines()
        # Extraer el tiempo del segundo dato de la primera línea que contiene metadatos
        time = float(lines[0].split()[1])
        nfluid = int(lines[0].split()[3])
        # Saltar la primera línea que contiene metadatos
        lines = lines[1:]
        # Extraer las posiciones x e y, el tipo y las velocidades de cada partícula de cada línea
        data = np.array([[float(val) for val in line.split()[0:]] for line in lines])

        """
        # Extraer cada columna en vectores con sus nombres específicos
        ids = data[:, 0].astype(int)
        posicion_x = data[:, 1]
        posicion_y = data[:, 2]
        velocidad_x = data[:, 3]
        velocidad_y = data[:, 4]
        masa = data[:, 5]
        rho = data[:, 6]
        pressure = data[:, 7]
        internal_energy = data[:, 8]
        itype = data[:, 9].astype(int)
        hsml = data[:, 10]
        viscosities_particles = data[:, 11]
        """
        
    return time, nfluid, data

def update(frame):
    ax.clear()
    ax.set_xlim(-0.1, 1.7)  # Ajustar límites x según tus necesidades
    ax.set_ylim(-0.1, 1.1)  # Ajustar límites y según tus necesidades

    # Extraer las posiciones x, y y viscosidad de cada partícula
    positions_x = frame[1][:, 1]
    positions_y = frame[1][:, 2]
    viscosities = frame[1][:, 11] 
    itype = frame[1][:,9].astype(int)
    
    
    # Determinar el umbral de viscosidad (puedes ajustarlo según tus datos)
    viscosity_threshold = 0.01
    
    # Seleccionar índices de partículas con menor y mayor viscosidad
    low_viscosity_indices = np.where((viscosities < viscosity_threshold) & (itype == 1))[0]
    high_viscosity_indices = np.where((viscosities >= viscosity_threshold) & (itype == 1))[0]
    wall_particles_indices = np.where((itype == -1 ))[0]
    
    # Dibujar partículas de menor viscosidad (en azul)
    ax.scatter(positions_x[low_viscosity_indices], positions_y[low_viscosity_indices], color='blue', s=3, label='Water')
    
    # Dibujar partículas de mayor viscosidad (en café)
    ax.scatter(positions_x[high_viscosity_indices], positions_y[high_viscosity_indices], color='brown', s=3, label='Mud')

    # Dibujar partículas de pared (en gris)
    ax.scatter(positions_x[wall_particles_indices], positions_y[wall_particles_indices], color='gray', s=3)
    
    # Mostrar el tiempo en la gráfica
    time_factor = 1.0
    adjusted_time = frame[0] * time_factor
    ax.text(0.05, 0.95, f't = {adjusted_time:.2f}', transform=ax.transAxes, ha='left', va='top', fontsize=12)
    ax.legend(loc='upper right')

# Definir la carpeta donde se encuentran los archivos
folder = './'

# Crear una figura y un eje para la animación
fig, ax = plt.subplots()

# Lista para almacenar los datos de todas las instantáneas
data_frames = []

# Iterar sobre los nombres de archivo del tipo 'snapshot_0001', 'snapshot_0002', ..., 'snapshot_1600'
for i in range(1583):
    # Formatear el nombre de archivo con el número de snapshot
    filename = f'snapshot_{i:04d}'
    # Crear la ruta completa al archivo
    filepath = os.path.join(folder, filename)
    # Verificar si el archivo existe
    if os.path.exists(filepath):
        # Leer datos del archivo snapshot y almacenarlos en la lista
        time, nfluid, data = read_snapshot(filepath)
        data_frames.append((time, data))
    else:
        print(f'El archivo {filename} no existe en la carpeta especificada.')

# Crear la animación
anim = FuncAnimation(fig, update, frames=data_frames, interval=10)  # Intervalo de 10 milisegundos entre cada frame

# Guardar la animación como un archivo .mp4
anim.save('animacion.mp4', writer='ffmpeg')

plt.tight_layout()
plt.show()
