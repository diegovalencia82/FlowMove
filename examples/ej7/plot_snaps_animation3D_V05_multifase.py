import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.gridspec import GridSpec
from scipy.interpolate import griddata, RBFInterpolator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def interpolar_puntos(puntos, num_puntos=100, metodo='cubic'):
    """
    Interpola los puntos de la matriz en una cuadrícula regular.

    Parámetros:
    - puntos: numpy array con las coordenadas (X, Y, Z) de los puntos.
    - num_puntos: número de puntos en cada eje de la cuadrícula (opcional, por defecto 100).
    - metodo: método de interpolación (opcional, por defecto 'cubic').

    Retorna:
    - Xh_grid: cuadrícula de coordenadas X.
    - Yh_grid: cuadrícula de coordenadas Y.
    - Zh_grid: valores interpolados de Z en la cuadrícula.
    """
    xh_grid = np.linspace(np.min(puntos[:, 0]), np.max(puntos[:, 0]), num_puntos)
    yh_grid = np.linspace(np.min(puntos[:, 1]), np.max(puntos[:, 1]), num_puntos)
    Xh_grid, Yh_grid = np.meshgrid(xh_grid, yh_grid)
    Zh_grid = griddata((puntos[:, 0], puntos[:, 1]), puntos[:, 2], (Xh_grid, Yh_grid), method=metodo)
    
    if np.isnan(Zh_grid).any():
        print("La interpolación contiene valores NaN.")
    
    return Xh_grid, Yh_grid, Zh_grid

def interpolar_puntos_3d(puntos, num_puntos=100, kernel='thin_plate_spline'):
    """
    Interpola los puntos tridimensionales en una cuadrícula 3D regular utilizando RBFInterpolator.

    Parámetros:
    - puntos: numpy array con las coordenadas (X, Y, Z) de los puntos.
    - num_puntos: número de puntos en cada eje de la cuadrícula (opcional, por defecto 100).
    - kernel: tipo de kernel para la interpolación radial (opcional, por defecto 'thin_plate_spline').

    Retorna:
    - Xh_grid: cuadrícula de coordenadas X.
    - Yh_grid: cuadrícula de coordenadas Y.
    - Zh_grid: cuadrícula de coordenadas Z.
    - Zh_interpolados: valores interpolados de Z en la cuadrícula tridimensional.
    """
    xh_grid = np.linspace(np.min(puntos[:, 0]), np.max(puntos[:, 0]), num_puntos)
    yh_grid = np.linspace(np.min(puntos[:, 1]), np.max(puntos[:, 1]), num_puntos)
    zh_grid = np.linspace(np.min(puntos[:, 2]), np.max(puntos[:, 2]), num_puntos)
    Xh_grid, Yh_grid, Zh_grid = np.meshgrid(xh_grid, yh_grid, zh_grid)
    
    rbf_interpolator = RBFInterpolator(puntos[:, :3], puntos[:, 2], kernel=kernel)
    puntos_grid = np.vstack([Xh_grid.ravel(), Yh_grid.ravel(), Zh_grid.ravel()]).T
    Zh_interpolados = rbf_interpolator(puntos_grid).reshape(Xh_grid.shape)
    
    return Xh_grid, Yh_grid, Zh_grid, Zh_interpolados

def crear_caja(x_min, x_max, y_min, y_max, z_min, z_max):
    """
    Crea las caras de una caja 3D.

    Parámetros:
    - x_min, x_max: límites en el eje X.
    - y_min, y_max: límites en el eje Y.
    - z_min, z_max: límites en el eje Z.

    Retorna:
    - vertices: lista de caras de la caja.
    """
    vertices = [
        [(x_min, y_min, z_min), (x_max, y_min, z_min), (x_max, y_max, z_min), (x_min, y_max, z_min)],
        [(x_min, y_min, z_max), (x_max, y_min, z_max), (x_max, y_max, z_max), (x_min, y_max, z_max)],
        [(x_min, y_min, z_min), (x_max, y_min, z_min), (x_max, y_min, z_max), (x_min, y_min, z_max)],
        [(x_min, y_max, z_min), (x_max, y_max, z_min), (x_max, y_max, z_max), (x_min, y_max, z_max)],
        [(x_min, y_min, z_min), (x_min, y_max, z_min), (x_min, y_max, z_max), (x_min, y_min, z_max)],
        [(x_max, y_min, z_min), (x_max, y_max, z_min), (x_max, y_max, z_max), (x_max, y_min, z_max)],
    ]
    return vertices

def dibujar_caja(ax, vertices):
    """
    Dibuja una caja 3D en el gráfico.

    Parámetros:
    - ax: eje 3D donde se dibujará la caja.
    - vertices: lista de caras de la caja.
    """
    for cara in vertices:
        ax.add_collection3d(Poly3DCollection([cara], facecolors='cyan', edgecolors='black', alpha=0.2))

def read_snapshot(filepath):
    """
    Lee los datos de un archivo snapshot.

    Parámetros:
    - filepath: ruta del archivo snapshot.

    Retorna:
    - time: tiempo del snapshot.
    - data: array con los datos de las partículas.
    """
    with open(filepath, 'r') as file:
        lines = file.readlines()
        time = float(lines[0].split()[1])
        lines = lines[1:]
        data = np.array([[float(val) for val in line.split()] for line in lines])
    return time, data

def calculate_velocity(data):
    """
    Calcula la magnitud de la velocidad para cada partícula.

    Parámetros:
    - data: array con los datos de las partículas.

    Retorna:
    - velocity: array con la magnitud de la velocidad de cada partícula.
    """
    velocity = np.sqrt(data[:, 4]**2 + data[:, 5]**2 + data[:, 6]**2)
    return velocity

def update(frame):
    """
    Actualiza el gráfico para cada frame de la animación.
    """
    ax_3d.clear()
    ax_xz.clear()
    
    # Configuración de ejes (igual que antes)
    ax_3d.set_xlabel('X')
    ax_3d.set_ylabel('Y')
    ax_3d.set_zlabel('Z')
    
    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min
    max_range = max(x_range, y_range, z_range)

    ax_3d.set_xlim(x_min, x_min + max_range)
    ax_3d.set_ylim(y_min, y_min + max_range)
    ax_3d.set_zlim(z_min, z_min + max_range)
    ax_3d.set_box_aspect([1, 1, 1])

    ax_xz.set_xlim(x_min, x_max)
    ax_xz.set_ylim(z_min, z_max)
    ax_xz.set_xlabel('X')
    ax_xz.set_ylabel('Z')
    ax_xz.set_aspect('equal', 'box')
    
    vertices_caja = crear_caja(x_min, x_max, y_min, y_max, z_min, z_max)
    dibujar_caja(ax_3d, vertices_caja)
    
    if sel == 0:
        data = frame[1][frame[1][:, 11] != -1]
    else:
        data = frame[1]
        
    fluid_particles = data[data[:, 11] == 1]
    mountain_particles1 = data[data[:, 11] == -1]
    mountain_particles = mountain_particles1[::5]
    
    # Separar partículas por viscosidad
    umbral_viscosidad = 0.1  # Ajustar este valor según necesidad
    mask_agua = fluid_particles[:, 13] < umbral_viscosidad
    mask_lodo = fluid_particles[:, 13] >= umbral_viscosidad
    
    agua_particles = fluid_particles[mask_agua]
    lodo_particles = fluid_particles[mask_lodo]
    
    # Calcular velocidades
    velocity_agua = calculate_velocity(agua_particles)
    velocity_lodo = calculate_velocity(lodo_particles)
    
    # Crear normas y mapas de color separados
    norm = Normalize(vmin=0, vmax=np.sqrt(2*9.8*ht))
    cmap_agua = plt.get_cmap('winter')
    cmap_lodo = plt.get_cmap('copper')
    
    # Generar colores
    colors_agua = cmap_agua(norm(velocity_agua))
    colors_lodo = cmap_lodo(norm(velocity_lodo))
    
    # Graficar partículas de agua
    ax_3d.scatter(agua_particles[:, 1], agua_particles[:, 2], agua_particles[:, 3], 
                 color=colors_agua, s=tamanop, alpha=0.8, label='Agua')
    
    # Graficar partículas de lodo
    ax_3d.scatter(lodo_particles[:, 1], lodo_particles[:, 2], lodo_particles[:, 3], 
                 color=colors_lodo, s=tamanop, alpha=0.8, label='Lodo')
    
    if sel == 1:
        ax_xz.scatter(mountain_particles[:, 1], mountain_particles[:, 3], color='green', s=0.1, alpha=0.1)
    
    # Graficar proyección X-Z
    ax_xz.scatter(agua_particles[:, 1], agua_particles[:, 3], color=colors_agua, s=tamanop, alpha=0.8)
    ax_xz.scatter(lodo_particles[:, 1], lodo_particles[:, 3], color=colors_lodo, s=tamanop, alpha=0.8)

    if sel == 1:
        Xh_grid, Yh_grid, Zh_grid = interpolar_puntos(mountain_particles[:, 1:4], num_puntos=100, metodo='cubic')
        ax_3d.plot_surface(Xh_grid, Yh_grid, Zh_grid, cmap='Greens', alpha=0.7, edgecolor='none')
    
    ax_3d.text2D(0.05, 0.95, f't = {frame[0]:.2f}', transform=ax_3d.transAxes, ha='left', va='top', fontsize=12)
    ax_3d.legend(loc='upper right')

def update1(frame):
    """
    Actualiza el gráfico para cada frame de la animación.

    Parámetros:
    - frame: tupla con el tiempo y los datos de las partículas.
    """
    ax_3d.clear()
    ax_xz.clear()
    
    ax_3d.set_xlabel('X')
    ax_3d.set_ylabel('Y')
    ax_3d.set_zlabel('Z')
    
    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min
    max_range = max(x_range, y_range, z_range)

    ax_3d.set_xlim(x_min, x_min + max_range)
    ax_3d.set_ylim(y_min, y_min + max_range)
    ax_3d.set_zlim(z_min, z_min + max_range)
    ax_3d.set_box_aspect([1, 1, 1])

    ax_xz.set_xlim(x_min, x_max)
    ax_xz.set_ylim(z_min, z_max)
    ax_xz.set_xlabel('X')
    ax_xz.set_ylabel('Z')
    ax_xz.set_aspect('equal', 'box')
    
    vertices_caja = crear_caja(x_min, x_max, y_min, y_max, z_min, z_max)
    dibujar_caja(ax_3d, vertices_caja)
    
    if sel == 0:
        data = frame[1][frame[1][:, 11] != -1]
    else:
        data = frame[1]
        
    fluid_particles = data[data[:, 11] == 1]
    mountain_particles1 = data[data[:, 11] == -1]
    mountain_particles = mountain_particles1[::5]
    
    velocity = calculate_velocity(fluid_particles)
    
    norm = Normalize(vmin=0, vmax=np.sqrt(2*9.8*ht))
    cmap = plt.get_cmap('copper')
    sm = ScalarMappable(norm=norm, cmap=cmap)
    colors = sm.to_rgba(velocity)
    
    ax_3d.scatter(fluid_particles[:, 1], fluid_particles[:, 2], fluid_particles[:, 3], color=colors, s=tamanop, alpha=0.8)
    
    if sel == 1:
        ax_xz.scatter(mountain_particles[:, 1], mountain_particles[:, 3], color='green', s=0.1, alpha=0.1)
    
    ax_xz.scatter(fluid_particles[:, 1], fluid_particles[:, 3], color=colors, s=tamanop, alpha=0.8)

    if sel == 1:
        Xh_grid, Yh_grid, Zh_grid = interpolar_puntos(mountain_particles[:, 1:4], num_puntos=100, metodo='cubic')
        ax_3d.plot_surface(Xh_grid, Yh_grid, Zh_grid, cmap='Greens', alpha=0.7, edgecolor='none')
    
    ax_3d.text2D(0.05, 0.95, f't = {frame[0]:.2f}', transform=ax_3d.transAxes, ha='left', va='top', fontsize=12)


def main():
    global sel, tamanop, ht, x_min, x_max, y_min, y_max, z_min, z_max, ax_3d, ax_xz

    # Entradas del usuario
    sel = int(input("Mostrar partículas de virtuales (1 para sí, 0 para no) (presiona Enter para usar 0): ") or "0")
    tamanop = float(input("Tamaño de las partículas (presiona Enter para usar 1): ") or "1")
    ht = float(input("Altura total (ht) (presiona Enter para usar 0.5): ") or "0.5")
    x_min = float(input("Límite mínimo en X (presiona Enter para usar 0.0): ") or "0.0")
    x_max = float(input("Límite máximo en X (presiona Enter para usar 1.2): ") or "1.2")
    y_min = float(input("Límite mínimo en Y (presiona Enter para usar 0.0): ") or "0.0")
    y_max = float(input("Límite máximo en Y (presiona Enter para usar 0.5): ") or "0.5")
    z_min = float(input("Límite mínimo en Z (presiona Enter para usar 0.0): ") or "0.0")
    z_max = float(input("Límite máximo en Z (presiona Enter para usar 0.5): ") or "0.5")

    folder = input("Carpeta de los archivos snapshot (Presiona Enter para usar ./): " or "./")
    num_snapshots = int(input("Número de snapshots a leer (presiona Enter para usar 300): ") or "301")

    fig = plt.figure(figsize=(12, 6))
    gs = GridSpec(1, 2, width_ratios=[2, 1])
    ax_3d = fig.add_subplot(gs[0], projection='3d')
    ax_xz = fig.add_subplot(gs[1])

    data_frames = []
    for i in range(0, num_snapshots):
        filename = f'snapshot_{i:04d}'
        filepath = os.path.join(folder, filename)
        if os.path.exists(filepath):
            print('Loading = ', filename)
            time, data = read_snapshot(filepath)
            data_frames.append((time, data))
        else:
            print(f'El archivo {filename} no existe en la carpeta especificada.')

    anim = FuncAnimation(fig, update, frames=data_frames, interval=25)
    
    # Crear ejes para las barras (ajustar posiciones según necesidad)
    cbar_ax_agua = fig.add_axes([0.15, 0.92, 0.3, 0.02])  # Barra para agua (arriba-izquierda)
    cbar_ax_lodo = fig.add_axes([0.55, 0.92, 0.3, 0.02])  # Barra para lodo (arriba-derecha)
    
    norm = Normalize(vmin=0, vmax=np.sqrt(2*9.8*ht))
    
    # Barra para agua (winter)
    cmap_agua = plt.get_cmap('winter')
    sm_agua = ScalarMappable(norm=norm, cmap=cmap_agua)
    cbar_agua = plt.colorbar(sm_agua, cax=cbar_ax_agua, orientation='horizontal')
    cbar_ax_agua.text(0.5, 1.5, 'Agua', transform=cbar_ax_agua.transAxes, 
                     ha='center', va='bottom', fontsize=10)

    # Barra para lodo (copper)
    cmap_lodo = plt.get_cmap('copper')
    sm_lodo = ScalarMappable(norm=norm, cmap=cmap_lodo)
    cbar_lodo = plt.colorbar(sm_lodo, cax=cbar_ax_lodo, orientation='horizontal')
    cbar_ax_lodo.text(0.5, 1.5, 'Lodo', transform=cbar_ax_lodo.transAxes, 
                     ha='center', va='bottom', fontsize=10)

   # Ajustar el layout para dejar espacio (modificar el rect)
    plt.tight_layout(rect=[0, 0, 1, 0.85])  # 15% de espacio superior

    anim.save('animacion_3D.mp4', writer='ffmpeg')
    #plt.tight_layout(rect=[0, 0, 1, 0.9])
    plt.show()

if __name__ == "__main__":
    main()
