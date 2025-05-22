#!/bin/bash

# Asegurar que se recibe un argumento
if [ -z "$1" ]; then
    echo "Error: Debes especificar el nombre del directorio de destino."
    echo "Uso: $0 <ruta_completa_del_directorio>"
    exit 1
fi

EJEMPLO=$1  # Asignar el argumento a la variable

# Verificar si el directorio existe
if [ ! -d "$EJEMPLO" ]; then
    echo "Error: El directorio examples/$EJEMPLO no existe."
    exit 1
fi

# Copiar el archivo de parámetros
cp "$EJEMPLO/param.inc" include/.

# Limpiar compilación anterior y compilar el código
make clean
make sph

# Mover el ejecutable "sph" al directorio correspondiente
mv sph "$EJEMPLO"

echo "El ejecutable se ha movido a $EJEMPLO/"
