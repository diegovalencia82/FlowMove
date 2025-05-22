import numpy as np
import matplotlib.pyplot as plt

# Parámetros para el espacio de partículas
particle_spacing = 0.008  # Espaciado entre partículas del fluido
delta = particle_spacing  # Pequeño delta para no sobreponerse con las paredes

# Parámetros de la caja
box_length = 1.61
box_height = 0.61

# Constantes físicas
g = 9.82
gamma = 7.0
water_height = 0.2 + particle_spacing
ht = water_height
beta0 = 25.0
beta = beta0 * 1.0
c = beta * np.sqrt(2. * g * ht)
rho0 = 1000.
b = rho0 * c * c / gamma

def create_box(particle_spacing, box_length, box_height, vel_x, vel_z, viscosity, n_parameter):
    """Crear partículas para las paredes de la caja y asignar velocidades uniformes."""
    x_box = np.arange(0, box_length, particle_spacing / 2)
    z_box = np.arange(particle_spacing / 2., box_height, particle_spacing / 2)
    
    wall_particles = []
    # Paredes laterales
    for x in [0, box_length]:
        for z in z_box:
            wall_particles.append([x, z])
    # Piso
    for x in x_box:
        wall_particles.append([x, 0])
    
    wall_particles = np.array(wall_particles)

    # Crear un array de velocidades, todas con el mismo valor
    velocities = np.full((wall_particles.shape[0], 2), [vel_x, vel_z])
    # Crea viscosidad para la pared, se puede transformar en rugosidad
    viscosities =  np.full(wall_particles.shape[0], viscosity)
    n_parameters = np.full(wall_particles.shape[0], n_parameter)

    return wall_particles, velocities, viscosities, n_parameters

def create_fluid_block(particle_spacing, delta, water_length, water_height, x_inicio, y_inicio, vel_x, vel_z, density, viscosity, n_parameter):
    """
    Crear partículas para el bloque de agua y asignar velocidades, densidades y viscosidades uniformes.
    El bloque comienza desde el punto (x_inicio, y_inicio).
    
    Parámetros:
    - particle_spacing: Espaciado entre partículas.
    - delta: Parámetro de ajuste para el inicio.
    - water_length: Longitud del bloque de agua.
    - water_height: Altura del bloque de agua.
    - x_inicio, y_inicio: Coordenadas de inicio del bloque de agua.
    - vel_x, vel_z: Velocidades en x y z.
    - density: Valor de la densidad de las partículas.
    - viscosity: Valor de la viscosidad de las partículas.
    
    Retorna:
    - particles: Coordenadas de las partículas.
    - velocities: Velocidades de cada partícula.
    - densities: Vector de densidades de cada partícula.
    - viscosities: Vector de viscosidades de cada partícula.
    """
    # Generar coordenadas de las partículas
    x_water = np.arange(x_inicio, x_inicio + water_length, particle_spacing)
    z_water = np.arange(y_inicio, y_inicio + water_height + particle_spacing, particle_spacing)
    xx_water, zz_water = np.meshgrid(x_water, z_water)
    particles = np.vstack([xx_water.ravel(), zz_water.ravel()]).T

    # Crear un array de velocidades, densidades, viscosidades y exponente n
    velocities = np.full((particles.shape[0], 2), [vel_x, vel_z])
    densities = np.full(particles.shape[0], density)
    viscosities = np.full(particles.shape[0], viscosity)
    n_parameters = np.full(particles.shape[0], n_parameter)

    return particles, velocities, densities, viscosities, n_parameters


def crear_circulo_particulas(x_centro, y_centro, ancho, alto, espaciado, radio, vel_x, vel_y):
    """
    Crea un bloque de partículas rectangulares centrado en (x_centro, y_centro) 
    y filtra las partículas que estén dentro de un radio determinado, asignando
    una velocidad constante a cada partícula.

    Parámetros:
    - x_centro, y_centro: Coordenadas del centro del bloque.
    - ancho, alto: Dimensiones del bloque.
    - espaciado: Espaciado entre las partículas.
    - radio: Radio máximo para filtrar partículas desde el centro del bloque.
    - vel_x, vel_y: Componentes de la velocidad para cada partícula.

    Retorna:
    - particulas_filtradas: Array con las coordenadas de las partículas dentro del radio.
    - velocidades: Array con las velocidades de cada partícula.
    """
    # Generar las coordenadas de las partículas
    x = np.arange(x_centro - ancho/2, x_centro + ancho/2 + espaciado, espaciado)
    y = np.arange(y_centro - alto/2, y_centro + alto/2 + espaciado, espaciado)
    xx, yy = np.meshgrid(x, y)
    particulas = np.vstack([xx.ravel(), yy.ravel()]).T

    # Filtrar partículas dentro del radio
    distancias = np.sqrt((particulas[:, 0] - x_centro)**2 + (particulas[:, 1] - y_centro)**2)
    particulas_filtradas = particulas[distancias <= radio]

    # Crear un array de velocidades, todas con el mismo valor
    velocidades = np.full((particulas_filtradas.shape[0], 2), [vel_x, vel_y])

    return particulas_filtradas, velocidades


def setup_particles(particle_spacing, water_particles, wall_particles):
    """Inicializar y calcular propiedades de las partículas."""
    num_fluid_particles = water_particles.shape[0]
    num_box_particles = wall_particles.shape[0]
    num_particles = num_fluid_particles + num_box_particles
    
    ids_fluid = np.arange(1, num_fluid_particles + 1)
    ids_wall = np.arange(num_fluid_particles + 1, num_particles + 1)
#    velocities = np.zeros((num_particles, 2))
    
    dx = particle_spacing
    rho = np.zeros(num_particles)
    pressure = np.zeros(num_particles)
    mass = np.zeros(num_particles)
    
    for i in range(num_fluid_particles):
        rho[i] = rho0 * (1 + (rho0 * g * (ht - abs(water_particles[i, 1])) / b)) ** (1.0 / gamma)
        pressure[i] = b * ((rho[i] / rho0) ** gamma - 1.0)
        mass[i] = dx * dx * rho[i]

    # Para las partículas de pared, mantener rho0
    rho[num_fluid_particles:] = rho0
    pressure[num_fluid_particles:] = 0.0
    mass[num_fluid_particles:] = dx * dx * rho0

    internal_energy = np.full(num_particles, 357.1)
    itype = np.full(num_particles, -1)
    itype[:num_fluid_particles] = 1  # Partículas del fluido
    hsml = np.full(num_particles, dx)
    
    return ids_fluid, ids_wall, rho, pressure, mass, internal_energy, itype, hsml

def write_snapshot(filename, ids_fluid, ids_wall, water_particles, wall_particles, velocities, mass, rho, pressure, internal_energy, itype, hsml, viscosities_particles, n_exponent):
    """Guardar los datos de las partículas en un archivo."""
    num_fluid_particles = water_particles.shape[0]
    num_box_particles = wall_particles.shape[0]
    num_particles = num_fluid_particles + num_box_particles
    
    with open(filename, 'w') as f:
        f.write(f"0   0.0   {num_particles}   {num_fluid_particles}   {num_box_particles}\n")

        # Escribir partículas del fluido
        for i in range(num_fluid_particles):
            f.write(f"{ids_fluid[i]}   {water_particles[i, 0]:.6f}   {water_particles[i, 1]:.6f}   "
                    f"{velocities[i, 0]:.6f}   {velocities[i, 1]:.6f}   {mass[i]:.6f}   {rho[i]:.6f}   "
                    f"{pressure[i]:.6f}   {internal_energy[i]:.6f}   {itype[i]}   {hsml[i]:.6f}   {viscosities_particles[i]:.6f}   0.0   {n_exponent[i]:.6}\n")

        # Escribir partículas de las paredes
        for i in range(num_box_particles):
            j = num_fluid_particles + i
            f.write(f"{ids_wall[i]}   {wall_particles[i, 0]:.6f}   {wall_particles[i, 1]:.6f}   "
                    f"{velocities[j, 0]:.6f}   {velocities[j, 1]:.6f}   {mass[j]:.6f}   {rho[j]:.6f}   "
                    f"{pressure[j]:.6f}   {internal_energy[j]:.6f}   {itype[j]}   {hsml[j]:.6f}   {viscosities_particles[j]:.6f}  0.0  {n_exponent[j]:.6}\n")

def plot_particles(water_particles, wall_particles):
    """Graficar las partículas en 2D."""
    fig, ax = plt.subplots()

    # Paredes de la caja (puntos más pequeños y transparentes)
    ax.scatter(wall_particles[:, 0], wall_particles[:, 1], c='r', s=10, alpha=0.1, label='Box')

    # Bloque de agua (puntos más grandes)
    ax.scatter(water_particles[:, 0], water_particles[:, 1], c='b', s=10, label='Water')

    ax.set_xlabel('X')
    ax.set_ylabel('Z')
    ax.legend()
    ax.set_aspect('equal')  # Asegurar que los ejes tengan la misma escala
    ax.grid(True)
    plt.show()

# Uso de las funciones

vel_x_wall = 0.0
vel_z_wall = 0.0
viscosity = 0.0
n_wall = 1.0
wall_particles, wall_velocities, wall_viscosities, n_walls = create_box(particle_spacing, box_length, box_height, vel_x_wall, vel_z_wall, viscosity, n_wall)

# Parámetros del bloque de agua
x_inicio = 0.95
y_inicio = 0. + particle_spacing*1.5
water_length = 0.3
water_height = 0.6
vel_x = 0.0
vel_z = 0.0
density = rho0
viscosity = 0.001002
n_parameter = 1.0
water_particles, velocities_fluid, densities, viscosities, n_parameters = create_fluid_block(particle_spacing, delta, water_length, water_height, x_inicio, y_inicio, vel_x, vel_z, density, viscosity, n_parameter)

# Parámetros del bloque de lodo
x_inicio = 0.35
y_inicio = 0. + particle_spacing*1.5
water_length = 0.3
water_height = 0.4
vel_x = 0.0
vel_z = 0.0
density = rho0
viscosity = 0.5
n_parameter_lodo = 1.5
lodo_particles, velocities_lodo, densities_lodo, viscosities_lodo, n_parameters_lodo = create_fluid_block(particle_spacing, delta, water_length, water_height, x_inicio, y_inicio, vel_x, vel_z, density, viscosity, n_parameter_lodo)


'''
# Crea el circulo de particulas
x_centro = 1.1
y_centro = 0.4
radio = 0.1
ancho = radio * 5
alto = radio * 5
vel_x_circulo = 2.0
vel_z_circulo = 0.0
circulo_particulas, velocities_circulo = crear_circulo_particulas(x_centro, y_centro, ancho, alto, particle_spacing, radio, vel_x_circulo, vel_z_circulo)
'''

#todas_las_particulas = np.vstack([water_particles, circulo_particulas])
todas_las_particulas = np.vstack([water_particles,lodo_particles])
#velocities_particulas = np.vstack([velocities_fluid, velocities_circulo, wall_velocities])
velocities_particulas = np.vstack([velocities_fluid, velocities_lodo, wall_velocities])
viscosities_particles = np.concatenate([viscosities, viscosities_lodo, wall_viscosities])
n_exponent = np.concatenate([n_parameters, n_parameters_lodo, n_walls])

# Da parametros a las partículas
ids_fluid, ids_wall, rho, pressure, mass, internal_energy, itype, hsml = setup_particles(particle_spacing, todas_las_particulas, wall_particles)
write_snapshot('snapshot_000', ids_fluid, ids_wall, todas_las_particulas, wall_particles, velocities_particulas, mass, rho, pressure, internal_energy, itype, hsml, viscosities_particles, n_exponent)

plot_particles(todas_las_particulas, wall_particles)


