import numpy as np
import matplotlib.pyplot as plt

# ─── OPCIÓN: Incluir o no el bloque de agua ───────────────────────────────
INCLUIR_FLUIDO = False   # Pon False para que solo aparezcan las partículas frontera

# ======================= FUNCIONES MODIFICADAS PARA 2D =======================
def generar_lamina_puntos_2D(punto_inicial, ancho, altura, espacio, angulo=0):
    """Genera una lámina de puntos en 2D (plano X-Z)"""
    angulo_rad = np.radians(angulo)
    rot = np.array([
        [np.cos(angulo_rad), -np.sin(angulo_rad)],
        [np.sin(angulo_rad),  np.cos(angulo_rad)]
    ])
    num_x = int(ancho/espacio) + 1
    num_z = int(altura/espacio) + 1
    puntos = []
    for i in range(num_x):
        for j in range(num_z):
            x, z = i*espacio, j*espacio
            x_rot, z_rot = rot @ np.array([x, z])
            puntos.append([punto_inicial[0]+x_rot, 0.0, punto_inicial[2]+z_rot])
    return np.array(puntos)

def generar_bloque_2D(punto_inicial, ancho, altura, espacio):
    """Genera un bloque fluido en 2D (X-Z)"""
    num_x = int(ancho/espacio) + 1
    num_z = int(altura/espacio) + 1
    puntos = []
    for i in range(num_x):
        for j in range(num_z):
            puntos.append([punto_inicial[0]+i*espacio, 0.0, punto_inicial[2]+j*espacio])
    return np.array(puntos)

# ==================== ROTACION Y TRANSLACION DE PUNTOS ======================
def transformar_puntos(puntos, theta, phi, traslacion):
    theta = np.radians(theta)
    phi   = np.radians(phi)
    Rz = np.array([[ np.cos(phi), -np.sin(phi), 0],
                   [ np.sin(phi),  np.cos(phi), 0],
                   [          0,            0,  1]])
    Ry = np.array([[ np.cos(theta), 0, np.sin(theta)],
                   [             0, 1,            0],
                   [-np.sin(theta), 0, np.cos(theta)]])
    rot = Ry @ Rz
    return puntos @ rot.T + traslacion

# ======================= VISUALIZACIÓN 2D ===================================
def graficar_configuracion(paredes=None, fluido=None, titulo="Configuración 2D",
                           tamano_puntos=10):
    fig, ax = plt.subplots()
    if paredes is not None:
        ax.scatter(paredes[:,0], paredes[:,2], s=tamano_puntos, label='Paredes')
    if INCLUIR_FLUIDO and fluido is not None:
        ax.scatter(fluido[:,0], fluido[:,2], s=tamano_puntos, label='Fluido')
    ax.set_xlabel('X'); ax.set_ylabel('Z')
    ax.set_title(titulo); ax.grid(True); ax.set_aspect('equal')
    ax.legend()
    plt.show()

# ======================= GENERACIÓN DEL CANAL 2D ============================
delta = 0.006
pared_delta = delta/2
pared_altura = 0.61
largo = 1.61

# Paredes
pared_izq = generar_lamina_puntos_2D([0,0,pared_delta], pared_delta/2, pared_altura-pared_delta, pared_delta)
pared_der = generar_lamina_puntos_2D([largo,0,pared_delta], pared_delta/2, pared_altura-pared_delta, pared_delta)
base     = generar_lamina_puntos_2D([0,0,0], largo, pared_delta/2, pared_delta)
puntos_paredes = np.vstack([pared_izq, pared_der, base])

# Fluido (opcional)
if INCLUIR_FLUIDO:
    fluid_alto = 0.3   # ejemplo de altura de fluido
    fluid_largo = 0.5  # ejemplo de largo de fluido
    fluido = generar_bloque_2D([delta*1.4,0,delta*1.4],
                               fluid_largo-delta, fluid_alto-delta, delta)
    puntos_fluido = fluido
else:
    puntos_fluido = np.empty((0,3))

# Rotar y trasladar
theta = 0; traslacion = np.array([0,0,0])
puntos_paredes = transformar_puntos(puntos_paredes, theta, 0, traslacion)
if INCLUIR_FLUIDO:
    puntos_fluido = transformar_puntos(puntos_fluido, theta, 0, traslacion)

# Gráficas
graficar_configuracion(paredes=puntos_paredes, fluido=puntos_fluido,
                       titulo="Sólo Paredes" if not INCLUIR_FLUIDO else "Paredes y Fluido")

# ======================= GENERAR snapshot_000 ================================
def generar_sph_snapshot_2D(X_topo, Z_topo, X_fluido, Z_fluido, filename='snapshot_000'):
    num_wall  = len(X_topo)
    num_fluid = len(X_fluido)
    total     = num_wall + num_fluid
    rho0 = 1000.0; g = 9.82; gamma = 7.0
    b = rho0 * (25*np.sqrt(2*g*0.3))**2 / gamma  # ejemplo c, ht fijo a 0.3
    mass = np.full(total, rho0 * delta**2)
    rho  = np.zeros(total); pressure = np.zeros(total)
    # fluid
    for i in range(num_fluid):
        rho[i]      = rho0 * (1 + (rho0*g*(0.3 - abs(Z_fluido[i]))/b))**(1.0/gamma)
        pressure[i] = rho0*g*(0.3 - Z_fluido[i])
    with open(filename, 'w') as f:
        f.write(f"0 0.0 {total} {num_fluid} {num_wall}\n")
        for i in range(num_fluid):
            f.write(f"{i+1} {X_fluido[i]:.6f} {Z_fluido[i]:.6f} 0 0 {mass[i]:.6f} "
                    f"{rho[i]:.6f} {pressure[i]:.6f} 357.1 1 {delta:.6f} 0.001 0 0\n")
        for i in range(num_wall):
            idx = num_fluid + i + 1
            f.write(f"{idx} {X_topo[i]:.6f} {Z_topo[i]:.6f} 0 0 {mass[idx-1]:.6f} "
                    f"{rho0/2:.6f} 0 357.1 -1 {delta:.6f} 0 0 0\n")

X_topo   = puntos_paredes[:,0]; Z_topo   = puntos_paredes[:,2]
X_fluido = puntos_fluido[:,0]; Z_fluido = puntos_fluido[:,2]
generar_sph_snapshot_2D(X_topo, Z_topo, X_fluido, Z_fluido, filename='snapshot_000')

