import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Parámetros fijos
sel = 0
tamanop = 1
ht = 0.5
x_min, x_max = 0.0, 1.2
y_min, y_max = 0.0, 0.5
z_min, z_max = 0.0, 0.5
folder = "./"
snap_ids = [0, 75, 150, 225, 300]

def read_snapshot(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()
        time = float(lines[0].split()[1])
        lines = lines[1:]
        data = np.array([[float(val) for val in line.split()] for line in lines])
    return time, data

def calculate_velocity(data):
    return np.sqrt(data[:, 4]**2 + data[:, 5]**2 + data[:, 6]**2)

def crear_caja(x_min, x_max, y_min, y_max, z_min, z_max):
    return [
        [(x_min, y_min, z_min), (x_max, y_min, z_min), (x_max, y_max, z_min), (x_min, y_max, z_min)],
        [(x_min, y_min, z_max), (x_max, y_min, z_max), (x_max, y_max, z_max), (x_min, y_max, z_max)],
        [(x_min, y_min, z_min), (x_max, y_min, z_min), (x_max, y_min, z_max), (x_min, y_min, z_max)],
        [(x_min, y_max, z_min), (x_max, y_max, z_min), (x_max, y_max, z_max), (x_min, y_max, z_max)],
        [(x_min, y_min, z_min), (x_min, y_max, z_min), (x_min, y_max, z_max), (x_min, y_min, z_max)],
        [(x_max, y_min, z_min), (x_max, y_max, z_min), (x_max, y_max, z_max), (x_max, y_min, z_max)],
    ]

def dibujar_caja(ax, vertices):
    for cara in vertices:
        ax.add_collection3d(Poly3DCollection([cara], facecolors='cyan', edgecolors='black', alpha=0.2))

# Crear figura
fig = plt.figure(figsize=(10, 18))
gs = GridSpec(len(snap_ids)*2, 2, width_ratios=[2, 1], height_ratios=[1]*len(snap_ids)*2)

norm = Normalize(vmin=0, vmax=np.sqrt(2 * 9.8 * ht))
cmap_agua = plt.get_cmap('winter')
cmap_lodo = plt.get_cmap('copper')

# Barras de color
cbar_ax_agua = fig.add_axes([0.15, 0.96, 0.3, 0.02])
cbar_ax_lodo = fig.add_axes([0.55, 0.96, 0.3, 0.02])
sm_agua = ScalarMappable(norm=norm, cmap=cmap_agua)
sm_lodo = ScalarMappable(norm=norm, cmap=cmap_lodo)
plt.colorbar(sm_agua, cax=cbar_ax_agua, orientation='horizontal')
plt.colorbar(sm_lodo, cax=cbar_ax_lodo, orientation='horizontal')
cbar_ax_agua.set_title('Agua', fontsize=10)
cbar_ax_lodo.set_title('Lodo', fontsize=10)

# Graficar cada snapshot
for i, snap in enumerate(snap_ids):
    filename = f'snapshot_{snap:04d}'
    filepath = os.path.join(folder, filename)
    if not os.path.exists(filepath):
        continue

    time, data = read_snapshot(filepath)
    if sel == 0:
        data = data[data[:, 11] != -1]

    fluid = data[data[:, 11] == 1]
    mountain = data[data[:, 11] == -1][::5]

    umbral_viscosidad = 0.1
    agua = fluid[fluid[:, 13] < umbral_viscosidad]
    lodo = fluid[fluid[:, 13] >= umbral_viscosidad]

    vel_agua = calculate_velocity(agua)
    vel_lodo = calculate_velocity(lodo)

    colors_agua = cmap_agua(norm(vel_agua))
    colors_lodo = cmap_lodo(norm(vel_lodo))

    ax3d = fig.add_subplot(gs[2*i, 0], projection='3d')
    axxz = fig.add_subplot(gs[2*i+1, 1])

    # Vista 3D
    ax3d.set_title(f"t = {time:.2f} s", fontsize=9)
    ax3d.scatter(agua[:, 1], agua[:, 2], agua[:, 3], color=colors_agua, s=tamanop, alpha=0.8)
    ax3d.scatter(lodo[:, 1], lodo[:, 2], lodo[:, 3], color=colors_lodo, s=tamanop, alpha=0.8)
    ax3d.set_xlim(x_min, x_max)
    ax3d.set_ylim(y_min, y_max)
    ax3d.set_zlim(z_min, z_max)
    ax3d.set_xlabel("X")
    ax3d.set_ylabel("Y")
    ax3d.set_zlabel("Z")
    ax3d.view_init(elev=20, azim=-60)
    ax3d.set_box_aspect([1, 1, 1])
    dibujar_caja(ax3d, crear_caja(x_min, x_max, y_min, y_max, z_min, z_max))

    # Proyección X-Z
    axxz.scatter(agua[:, 1], agua[:, 3], color=colors_agua, s=tamanop, alpha=0.8)
    axxz.scatter(lodo[:, 1], lodo[:, 3], color=colors_lodo, s=tamanop, alpha=0.8)
    axxz.set_xlim(x_min, x_max)
    axxz.set_ylim(z_min, z_max)
    axxz.set_xlabel("X")
    axxz.set_ylabel("Z")
    axxz.set_aspect("equal", "box")

# Guardar figura final
plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig("fig_snapshots_multifase_5frames_vertical.png", dpi=300)
plt.show()

