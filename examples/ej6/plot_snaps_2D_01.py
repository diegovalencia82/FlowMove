import os
import numpy as np
import matplotlib.pyplot as plt

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
        
    return time, nfluid, data

# Definir la carpeta donde se encuentran los archivos
folder = './'

# Configurar la cuadrícula de subfiguras (e.g., 2 filas x 2 columnas)
nrows, ncols = 4, 4
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 8))

# Iterar sobre los nombres de archivo del tipo 'snapshot_0001', 'snapshot_0002', ..., 'snapshot_4000'
# y las subfiguras en el grid
j = 0
for idx, ax in enumerate(axes.flat):
    i = idx * 125 + 1  # Ajustar para seleccionar snapshots espaciados cada 125, por ejemplo 1, 125,...
    # Formatear el nombre de archivo con el número de snapshot
    filename = f'snapshot_{i:04d}'
    # Crear la ruta completa al archivo
    filepath = os.path.join(folder, filename)
    j = j + 1
    
    # Verificar si el archivo existe
    if os.path.exists(filepath):
        # Leer datos del archivo snapshot
        time, nfluid, data = read_snapshot(filepath)

        # Limpiar el gráfico antes de dibujar el nuevo (aunque no es necesario en este contexto)
        ax.clear()
        ax.set_xlim(-0.1, 1.7)  # Ajustar límites x según tus necesidades
        ax.set_ylim(-0.1, 0.7)  # Ajustar límites y según tus necesidades

        # Extraer las posiciones x, y y viscosidad de cada partícula
        positions_x = data[:, 1]
        positions_y = data[:, 2]
        viscosities = data[:, 11]
        itype = data[:, 9].astype(int)
        
        # Determinar el umbral de viscosidad (puedes ajustarlo según tus datos)
        viscosity_threshold = 0.01
        
        # Seleccionar índices de partículas con menor y mayor viscosidad
        low_viscosity_indices = np.where((viscosities < viscosity_threshold) & (itype == 1))[0]
        high_viscosity_indices = np.where((viscosities >= viscosity_threshold) & (itype == 1))[0]
        wall_particles_indices = np.where(itype == -1)[0]

        if j == 1:
            label_water = 'Water'
            label_mud = 'Mud'
        else:
            label_water = None
            label_mud = None
            
        # Dibujar partículas de menor viscosidad (en azul)
        ax.scatter(positions_x[low_viscosity_indices], positions_y[low_viscosity_indices], color='blue', s=1, label=label_water)
        
        # Dibujar partículas de mayor viscosidad (en café)
        ax.scatter(positions_x[high_viscosity_indices], positions_y[high_viscosity_indices], color='brown', s=1, label=label_mud)

        # Dibujar partículas de pared (en gris)
        ax.scatter(positions_x[wall_particles_indices], positions_y[wall_particles_indices], color='gray', s=1)

        # Mostrar el tiempo en la gráfica
        time_factor = 1.0
        adjusted_time = time * time_factor
        #ax.set_title(f't = {adjusted_time:.2f}', fontsize=10)
        #ax.legend(loc='upper right', fontsize=8)
        ax.text(0.1, 0.96, f't = {adjusted_time:.2f}', transform=ax.transAxes, ha='left', va='top', fontsize=10)

        if (j==1):
            ax.set_ylabel('Y [m]')
            ax.set_xticks([])  # Ocultar etiquetas del eje X
            ax.legend(loc='upper right', fontsize=8)
        if (j==5):
            ax.set_ylabel('Y [m]')
            ax.set_xticks([])  # Ocultar etiquetas del eje X
#            ax.tick_params(axis='x', labelbottom=False)  # Ocultar etiquetas del eje X
        if (j==9):
            ax.set_ylabel('Y [m]')
            ax.set_xticks([])  # Ocultar etiquetas del eje X
#            ax.tick_params(axis='x', labelbottom=False)  # Ocultar etiquetas del eje X
        if (j==13):
            ax.set_ylabel('Y [m]')
            ax.set_xlabel('X [m]')
        if (j==14 or j==15 or j==16):
            ax.set_xlabel('X [m]')
            ax.set_yticks([])  # Ocultar etiquetas del eje Y
        if (j==2 or j==3 or j==4 or j==6 or j==7 or j==8 or j==10 or j==11 or j==12 ):
            ax.set_xticks([])  # Ocultar etiquetas del eje X
            ax.set_yticks([])  # Ocultar etiquetas del eje Y
        

    else:
        ax.set_title(f'Archivo {filename} no encontrado', fontsize=10)
        ax.axis('off')

# Ajustar el espacio entre subplots para evitar superposición
plt.tight_layout()

# Guardar la figura como archivo .png
plt.savefig('fig_snapshots_MultFase_water001_mud1.png', dpi=300, bbox_inches='tight')
print('Guardado grid_snapshots.png')

# Mostrar la figura
plt.show()

# Cerrar la figura al finalizar
plt.close(fig)
