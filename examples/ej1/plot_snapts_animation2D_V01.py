import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

def calculate_velocity(data):
    """
    Calcula la magnitud de la velocidad para cada partícula (ignorando componente Y).
    """
    velocity = np.sqrt(data[:, 3]**2 + data[:, 4]**2)  # Usar componentes X y Z
    return velocity

def read_snapshot(filepath):
    """
    Lee los datos de un archivo snapshot.
    """
    with open(filepath, 'r') as file:
        lines = file.readlines()
        time = float(lines[0].split()[1])
        lines = lines[1:]
        data = np.array([[float(val) for val in line.split()] for line in lines])
    return time, data

def update(frame):
    """
    Actualiza el gráfico 2D para cada frame de la animación.
    """
    ax_xz.clear()
    
    # Configurar límites y etiquetas
    ax_xz.set_xlim(x_min, x_max)
    ax_xz.set_ylim(z_min, z_max)
    ax_xz.set_xlabel('X')
    ax_xz.set_ylabel('Z')
    ax_xz.set_aspect('equal')
    
    # Filtrar datos según selección
    if sel == 0:
        data = frame[1][frame[1][:, 9] != -1]
    else:
        data = frame[1]
        
    fluid_particles = data[data[:, 9] == 1]
    mountain_particles = data[data[:, 9] == -1][::5]  # Submuestreo
    
    # Calcular velocidad y colores
    velocity = calculate_velocity(fluid_particles)
    colors = cmap(norm(velocity))

    # Dibujar partículas
    if sel == 1:
        ax_xz.scatter(mountain_particles[:, 1], mountain_particles[:, 2], 
                     color='green', s=0.1, alpha=0.1, label='Montaña')


    if mespejo == 1:
        ax_xz.scatter(x_max-fluid_particles[:, 1], fluid_particles[:, 2], 
                      color=colors, s=tamanop, alpha=0.8, label='Fluido')
    else: 
        ax_xz.scatter(fluid_particles[:, 1], fluid_particles[:, 2], 
                      color=colors, s=tamanop, alpha=0.8, label='Fluido')
    
    # Texto con tiempo
    ax_xz.text(0.05, 0.95, f't = {frame[0]:.2f}', 
              transform=ax_xz.transAxes, ha='left', va='top', fontsize=12)

def main():
    global sel, tamanop, ht, x_min, x_max, z_min, z_max, ax_xz, norm, cmap, mespejo

    # Configuración de usuario
    sel = int(input("Mostrar partículas virtuales (1 para sí, 0 para no): ") or "0")
    tamanop = float(input("Tamaño de partículas (presiona Enter para usar 1): ") or "1")
    mespejo = float(input("Modo normal 0 espejo 1 (presiona Enter para usar 0): ") or "0")
    ht = float(input("Altura total (ht) (Enter para usar 0.25): ") or "0.25")
    x_min = float(input("Límite mínimo en X (Enter para usar 0.0 ): ") or "0.0")
    x_max = float(input("Límite máximo en X (Enter para usar 1.61): ") or "1.61")
    z_min = float(input("Límite mínimo en Z (Enter para usar 0.0 ): ") or "0.0")
    z_max = float(input("Límite máximo en Z (Enter para usar 0.6 ): ") or "0.6")
    
    folder = input("Carpeta de snapshots (Enter para usar ./): ") or "./"
    num_snapshots = int(input("Número de snapshots (Enter para usar 2000): ") or "2000")

    norm = Normalize(vmin=0, vmax=np.sqrt(2*9.8*ht))
    cmap = plt.get_cmap('winter')

    # Configurar figura
    fig = plt.figure(figsize=(8, 6))
    ax_xz = fig.add_subplot(111)
    
    # Cargar datos
    data_frames = []
    for i in range(num_snapshots+1):
        filepath = os.path.join(folder, f'snapshot_{i:04d}')
        print('Loading ',f'snapshot_{i:04d}')
        if os.path.exists(filepath):
            time, data = read_snapshot(filepath)
            data_frames.append((time, data))
    
    # Configurar animación
    anim = FuncAnimation(fig, update, frames=data_frames, interval=25)
    
    # Barra de color
    cbar_ax = fig.add_axes([0.15, 0.92, 0.7, 0.02])
    sm = ScalarMappable(norm=Normalize(0, np.sqrt(2*9.8*ht)), cmap='winter')
    plt.colorbar(sm, cax=cbar_ax, orientation='horizontal').set_label('Magnitud de Velocidad')
    
    # Guardar y mostrar
    anim.save('animacion_2D.mp4', writer='ffmpeg')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

if __name__ == "__main__":
    main()
