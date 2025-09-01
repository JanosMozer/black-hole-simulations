# Black Hole Simulations

This project provides C++ simulations of gravitational lensing and black hole visualization using OpenGL. It includes both a 3D visualization and a simpler 2D simulation.

## Features

*   **3D Black Hole Simulation (`BlackHole3D`)**:
    *   **Primary Visualization (Compute Shader)**: An advanced ray-tracing simulation of gravitational lensing around a black hole. This is the intended visualization but requires OpenGL 4.3+.
    *   **Fallback 3D Visualization**: For systems that do not support compute shaders (like macOS), a fallback renderer is used. This displays a 3D scene with the central black hole, nearby stars, and a grid that visually deforms to represent the curvature of spacetime. The camera is fully interactive.
    *   **Physics**: Simulates the gravitational attraction between celestial bodies.

*   **2D Lensing Simulation (`BlackHole2D`)**:
    *   A simplified 2D visualization of gravitational lensing effects.
    *   Note: This simulation uses deprecated OpenGL functions and may not work correctly on modern OpenGL core profiles.

## macOS Compatibility Notice

This project has known compatibility issues on macOS systems. The primary 3D visualization relies on OpenGL compute shaders, which require **OpenGL 4.3** or higher. macOS officially supports up to **OpenGL 4.1**.

Due to this limitation, the project will automatically use a fallback 3D rendering mode on Mac systems. This fallback uses a CPU-based ray tracer to accurately simulate the gravitational lensing effect. However, because this is computationally intensive and does not use the GPU, the rendering resolution has been lowered to maintain usable performance. Even with this adjustment, the simulation was quite slow during testing on a MacBook, so the 3D visualization may not be 100% fluid or accurate.

## Dependencies

The project's dependencies are managed using [vcpkg](https://github.com/microsoft/vcpkg).
*   **GLEW**: The OpenGL Extension Wrangler Library.
*   **GLFW3**: A multi-platform library for OpenGL, window, and input.
*   **GLM**: OpenGL Mathematics, a header-only C++ mathematics library for graphics software.

## Build Instructions

1.  **Clone the repository**
    ```bash
    git clone https://github.com/JanosMozer/black-hole-simulations
    cd black-hole-simulations
    ```

2.  **Install vcpkg and dependencies**
    If you don't have vcpkg, clone it into the project directory.
    ```bash
    git clone https://github.com/Microsoft/vcpkg.git
    ./vcpkg/bootstrap-vcpkg.sh
    ```
    Then, install the required libraries.
    ```bash
    ./vcpkg/vcpkg install
    ```

3.  **Configure and Build with CMake**
    Create a build directory and run CMake to configure the project.
    ```bash
    cmake -B build -S . -DCMAKE_TOOLCHAIN_FILE=./vcpkg/scripts/buildsystems/vcpkg.cmake
    ```
    Then, build the executables.
    ```bash
    cmake --build build
    ```

## Usage

The executables will be located in the `build` directory.

*   To run the 3D simulation:
    ```bash
    ./build/BlackHole3D
    ```

*   To run the 2D simulation:
    ```bash
    ./build/BlackHole2D
    ```

### 3D Simulation Controls
*   **Mouse Drag**: Orbit the camera around the central black hole.
*   **Scroll Wheel**: Zoom in and out.

## Acknowledgements

Special thanks to **kavan010**, from whom I copied the OpenGL configurations that helped get this project running.
