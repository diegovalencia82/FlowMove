import os
import numpy as np
import matplotlib.pyplot as plt

def read_snapshot(filepath):
    """Lee archivos snapshot ignorando coordenada Y y velocidad Vy"""
    with open(filepath, 'r') as file:
        lines = file.readlines()
        time = float(lines[0].split()[1])
        lines = lines[1:]
        # Leer solo 13 columnas relevantes (ignorando Y y Vy)
        data = np.array([[float(val) for val in line.split()[:15]] for line in lines])
    return time, data

def read_data_file(filepath):
    """Lee datos experimentales"""
    with open(filepath, 'r') as file:
        lines = file.readlines()
        if lines[0].startswith("#"):
            lines = lines[1:]
        data = np.array([[float(val) for val in line.split()] for line in lines])
    return data[:, 0], data[:, 1]

def analyze_snapshot(data):
    """Analiza datos en 2D (X-Z) para flujo 3D"""
    # Nuevos índices para datos 3D filtrados
    ids = data[:, 0]
    xs = data[:, 1]   # X
    zs = data[:, 3]   # Z (índice 3 en 3D)
    vxs = data[:, 4]  # Vx
    vzs = data[:, 6]  # Vz
    masses = data[:, 7]
    rhos = data[:, 8]
    ps = data[:, 9]
    us = data[:, 10]
    itypes = data[:, 11]
    hsmls = data[:, 12]
    eta_c = data[:, 13]
    

    # Filtrar partículas de fluido
    fluid_mask = itypes != -1
    fluid_particles = data[fluid_mask]
    
    # Obtener coordenadas X-Z
    x_fluid = fluid_particles[:, 1]
    z_fluid = fluid_particles[:, 3]
    max_x = np.max(x_fluid)
    max_z = np.max(z_fluid)
    
    # Calcular distancias a máximo X
    distances = np.abs(x_fluid - max_x)
    closest_indices = np.argsort(distances)[:2]
    closest_particles = fluid_particles[closest_indices]
    
    return np.mean(closest_particles[:, 1]), max_x, max_z

def analizar_residuales(residuals):
    """Calcula estadísticas de residuales (igual que en 2D)"""
    resultados = {
        'promedio': np.mean(residuals),
        'desviacion_std': np.std(residuals, ddof=1),
        'varianza': np.var(residuals, ddof=1),
        'rango': np.max(residuals) - np.min(residuals),
        'mad': np.mean(np.abs(residuals - np.mean(residuals))),
        'iqr': np.percentile(residuals, 75) - np.percentile(residuals, 25)
    }
    return resultados

# Configuración principal
folder = './'
H, B = 0.3, 0.6  # Valores de referencia
avg_x_values = []
time_values = []

# Procesar snapshots
for i in range(140):
    filename = f'snapshot_{i:04d}'
    filepath = os.path.join(folder, filename)
    
    if os.path.exists(filepath):
        time, data = read_snapshot(filepath)
        avg_x, max_x, max_z = analyze_snapshot(data)
        avg_x_values.append(avg_x)
        time_values.append(time)

# Convertir a arrays y normalizar
time_array = np.array(time_values) * (9.82 / H)**0.5
x_array = (np.array(avg_x_values) - B) / H

# Leer datos experimentales
experiment_time, experiment_x = read_data_file('Fig12_01_ETSIN_Front_Dressler_Martin_2.dat')

# Configuración de estilo mejorada
plt.rcParams.update({
    'font.size': 18,
    'axes.labelsize': 20,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 18,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white'
})

# Crear figura única
fig, ax = plt.subplots(figsize=(10, 7))

# Graficar datos principales
ax.plot(time_array, x_array, 'k-', linewidth=2, label='Model 3D')
ax.plot(experiment_time, experiment_x, 'o', color='gray', 
        markersize=5, markerfacecolor='none', markeredgewidth=1.5,
        label='Experimental data')

# Configurar ejes y cuadrícula
ax.set(xlabel='$t(g/H)^{1/2}$', ylabel='$x/H$', 
       xlim=(0, 2.7), ylim=(0, 3.6))
ax.grid(True, linestyle=':', alpha=0.7, color='gray')

# Calcular estadísticas de residuales
time_window = 0.1  # Ventana temporal para filtrado
filtered_x = [np.mean(x_array[(time_array >= t-time_window) & (time_array <= t+time_window)]) 
              for t in experiment_time]
residuals = filtered_x - experiment_x
stats = analizar_residuales(residuals)

# Añadir leyenda principal
legend = ax.legend(loc='upper left', frameon=True, edgecolor='black')
legend.get_frame().set_facecolor('white')

# Añadir cuadro estadístico
stats_text = (f'Residual statistics:\n'
              f'• Range: {stats["rango"]:.3f}\n'
              f'• σ: {stats["desviacion_std"]:.3f}\n'
              f'• Mean: {stats["promedio"]:.3f}')

ax.text(0.97, 0.05, stats_text, transform=ax.transAxes,
        ha='right', va='bottom', fontsize=16,
        bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.5'))

# Guardar y mostrar
plt.tight_layout()
plt.savefig('comparacion_3D.png', dpi=300, bbox_inches='tight')
plt.show()
