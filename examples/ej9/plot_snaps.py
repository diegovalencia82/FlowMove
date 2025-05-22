import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def interpolar_puntos(puntos, num_puntos=100, metodo='cubic'):
    xh_grid = np.linspace(np.min(puntos[:, 0]), np.max(puntos[:, 0]), num_puntos)
    yh_grid = np.linspace(np.min(puntos[:, 1]), np.max(puntos[:, 1]), num_puntos)
    Xh_grid, Yh_grid = np.meshgrid(xh_grid, yh_grid)
    Zh_grid = griddata((puntos[:, 0], puntos[:, 1]), puntos[:, 2], (Xh_grid, Yh_grid), method=metodo)
    return Xh_grid, Yh_grid, Zh_grid

def crear_caja(x_min, x_max, y_min, y_max, z_min, z_max):
    vertices = np.array([
        [x_min, y_min, z_min], [x_max, y_min, z_min],
        [x_max, y_max, z_min], [x_min, y_max, z_min],
        [x_min, y_min, z_max], [x_max, y_min, z_max],
        [x_max, y_max, z_max], [x_min, y_max, z_max]
    ])
    caras = [
        [0, 1, 2, 3], [4, 5, 6, 7],
        [0, 1, 5, 4], [2, 3, 7, 6],
        [0, 3, 7, 4], [1, 2, 6, 5]
    ]
    return vertices, caras

def dibujar_caja(ax, vertices_transformados, caras):
    caras_coords = [vertices_transformados[cara] for cara in caras]
    ax.add_collection3d(Poly3DCollection(
        caras_coords, facecolors='cyan', edgecolors='black', alpha=0.2
    ))

def read_snapshot(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()
        time = float(lines[0].split()[1])
        data = np.array([[float(val) for val in line.split()] for line in lines[1:]])
    return time, data

def calculate_velocity(data):
    return np.sqrt(data[:, 4]**2 + data[:, 5]**2 + data[:, 6]**2)

def transformar_puntos(puntos, theta=45, phi=0, traslacion=[0, 0, 1.8]):
    theta, phi = np.radians(theta), np.radians(phi)
    Rz = np.array([[np.cos(phi), -np.sin(phi), 0],
                   [np.sin(phi), np.cos(phi), 0],
                   [0, 0, 1]])
    Ry = np.array([[np.cos(theta), 0, np.sin(theta)],
                   [0, 1, 0],
                   [-np.sin(theta), 0, np.cos(theta)]])
    return np.dot(puntos, (Ry @ Rz).T) + traslacion

def plot_snapshot(ax_3d, ax_xz, frame, sel, tamanop, ht, x_min, x_max, y_min, y_max, z_min, z_max):
    for ax in [ax_3d, ax_xz]:
        ax.clear()
    
    # Configuración de ejes 3D
    max_range = max(x_max-x_min, y_max-y_min, z_max-z_min)
    ax_3d.set_xlim(x_min, x_min + max_range)
    ax_3d.set_ylim(y_min, y_min + max_range)
    ax_3d.set_zlim(z_min, z_min + max_range)
    ax_3d.set_box_aspect([1, 1, 1])
    ax_3d.set_xlabel('X')
    ax_3d.set_ylabel('Y')
    ax_3d.set_zlabel('Z')
    
    # Configuración de ejes XZ
    ax_xz.set_xlim(x_min, x_min + max_range)
    ax_xz.set_ylim(z_min, z_min + max_range)
    ax_xz.set_aspect('equal', 'box')
    ax_xz.set_xlabel('X')
    ax_xz.set_ylabel('Z')
    
    # Dibujar caja
    vertices, caras = crear_caja(x_min, x_max, y_min, y_max, z_min, z_max)
    dibujar_caja(ax_3d, transformar_puntos(vertices), caras)
    
    # Procesar datos
    data = frame[1] if sel else frame[1][frame[1][:, 11] != -1]
    fluid_particles = data[data[:, 11] == 1]
    mountain_particles = data[data[:, 11] == -1][::5]
    
    # Calcular velocidades
    velocity = calculate_velocity(fluid_particles)
    norm = Normalize(vmin=0, vmax=np.sqrt(2*9.8*ht))
    colors = plt.cm.copper(norm(velocity))
    
    # Graficar partículas
    ax_3d.scatter(fluid_particles[:,1], fluid_particles[:,2], fluid_particles[:,3], 
                 c=colors, s=tamanop, alpha=0.8)
    ax_xz.scatter(fluid_particles[:,1], fluid_particles[:,3], 
                 c=colors, s=tamanop, alpha=0.8)
    
    if sel:
        ax_xz.scatter(mountain_particles[:,1], mountain_particles[:,3], 
                     color='green', s=0.1, alpha=0.1)
        Xh_grid, Yh_grid, Zh_grid = interpolar_puntos(mountain_particles[:,1:4])
        ax_3d.plot_surface(Xh_grid, Yh_grid, Zh_grid, cmap='Greens', alpha=0.7)
    
    ax_3d.text2D(0.05, 0.95, f't = {frame[0]:.2f}', transform=ax_3d.transAxes,
                fontsize=10, ha='left', va='top')

def main():
    # Configuración de usuario (sin cambios)
    sel = int(input("Mostrar partículas virtuales (1 para sí, 0 para no) [0]: ") or "0")
    tamanop = float(input("Tamaño de las partículas [1]: ") or "1")
    ht = float(input("Altura total (ht) [0.5]: ") or "0.5")
    x_min = float(input("Límite mínimo en X [0.0]: ") or "0.0")
    x_max = float(input("Límite máximo en X [2.5]: ") or "2.5")
    y_min = float(input("Límite mínimo en Y [0.0]: ") or "0.0")
    y_max = float(input("Límite máximo en Y [0.3]: ") or "0.3")
    z_min = float(input("Límite mínimo en Z [0.0]: ") or "0.0")
    z_max = float(input("Límite máximo en Z [0.5]: ") or "0.5")
    folder = input("Carpeta de snapshots [./]: ") or "./"
    total_snaps = int(input("Número total de snapshots [300]: ") or "300") + 1

    # Cargar snapshots
    data_frames = []
    for i in range(total_snaps):
        filepath = os.path.join(folder, f'snapshot_{i:04d}')
        if os.path.exists(filepath):
            data_frames.append(read_snapshot(filepath))
            print(f'Cargado: snapshot_{i:04d}')
    
    # Seleccionar 5 snapshots equidistantes
    step = max(1, (len(data_frames) - 1) // 4)
    selected_frames = [data_frames[i] for i in range(0, len(data_frames), step)][:5]

    # Configurar figura
    fig = plt.figure(figsize=(10, 20))
    gs = GridSpec(5, 1, hspace=0.4)
    
    # Crear subplots
    for i, frame in enumerate(selected_frames):
        row_gs = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[i], width_ratios=[2, 1], wspace=0.1)
        ax3d = fig.add_subplot(row_gs[0], projection='3d')
        axxz = fig.add_subplot(row_gs[1])
        plot_snapshot(ax3d, axxz, frame, sel, tamanop, ht, x_min, x_max, y_min, y_max, z_min, z_max)
    
    # Barra de color
    cbar_ax = fig.add_axes([0.3, 0.94, 0.4, 0.015])
    sm = ScalarMappable(norm=Normalize(0, np.sqrt(2*9.8*ht)), cmap='copper')
    plt.colorbar(sm, cax=cbar_ax, orientation='horizontal', label='Magnitud de Velocidad')
    
    plt.savefig('fig_snapshots_incli.png', dpi=300, bbox_inches='tight')
    print("Figura guardada como 'fig_snapshots_incli.png'")
    plt.close()

if __name__ == "__main__":
    main()
