# Atom Simulator 3D

Simulación interactiva de un átomo en 3D usando OpenGL.

## Requisitos

- CMake 3.10+
- Compilador C++17
- GLFW3
- OpenGL 3.3+

## Instalación de dependencias (Ubuntu)

```bash
sudo apt install build-essential cmake libglfw3-dev libgl1-mesa-dev
```

## Compilación y ejecución

```bash
mkdir build
cd build
cmake ..
make
./opengl_app
```

## Controles

- **Mouse**: Arrastrar para rotar la vista
- **Scroll**: Acercar/Alejar
- **↑/↓**: Cambiar elemento atómico

## Estructura

El proyecto simula un átomo con:
- Núcleo (protones y neutrones)
- Nube electrónica con orbitales (s, p, d)
- Física básica de interacciones coulombianas