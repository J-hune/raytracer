# Raytracer in C++

This project implements a **Raytracer** in C++ using **OpenGL** and **GLUT** libraries to generate realistic renders of 3D scenes.
The raytracer supports various shapes and materials, along with advanced lighting techniques such as global illumination and caustics via **photon mapping**.
This allows the simulation of realistic lighting effects, such as reflections, refractions, and light diffusion on surfaces.

<img src="/img/github/preview.png" alt="Raytracer Preview" width="500"/>

## Features and Architecture

### Light Management and Visual Effects

- **Direct Lighting**: Computed using the **Phong** lighting model, which simulates the impact of light based on the position, orientation, and properties of each material. This model includes:
    - **Diffuse Reflection**: For light distribution on matte surfaces.
    - **Specular Reflection**: For shiny highlights on smooth surfaces.
- **Indirect Lighting**: Simulated using **photon mapping** to represent light bouncing between surfaces. This creates a more realistic rendering, as areas in direct shadow can still receive reflected light.
- **Caustics**: Complex light effects caused by transparent (glass) or reflective (mirror) materials are accurately simulated using **photon mapping**, enabling light concentration effects (like light passing through a glass of water).


### Material Models and Textures

Object materials are defined to replicate realistic physical properties:
- **Mirror**: Perfectly specular reflective surfaces, ideal for reflective objects like mirrors.
- **Glass**: A material with both reflection and refraction, calculated using the **Fresnel effect** to simulate intensity variation based on the angle of incidence.
- **Diffuse Material**: A standard material model supporting diffuse and specular reflection (via the **Blinn-Phong** model).


#### Textures and Skybox

- **Textures**: The raytracer supports applying textures to objects.
- **Skybox**: A skybox can be included to simulate background environments, creating a realistic ambiance for the 3D scene.

### Photon Mapping

Photon mapping is used to simulate **indirect lighting** and **caustics**:
- **Photon Emission**: Photons are emitted from light sources and interact with scene materials.
- **Photon Map Storage**: Photon interactions with objects are stored and organized in a KDTree, enabling efficient illumination estimation.
- **Illumination Estimation**: Using photons near each intersection point, the raytracer calculates indirect lighting and caustics, with a **Cone filter** applied for smoothing.

---

## Render Examples
### Render Breakdown
Final renders result from combining multiple passes:

| Direct Lighting (Phong)                    | Indirect Lighting                              |
|--------------------------------------------|------------------------------------------------|
| ![Direct Lighting](/img/github/direct.png) | ![Indirect Lighting](/img/github/indirect.png) |

| Caustics                              | Final Render                                 |
|---------------------------------------|----------------------------------------------|
| ![Caustics](/img/github/caustics.png) | ![Final Render](/img/github/final_small.png) |

### Render Examples with Skybox

| Skybox with Reflection                  | Skybox with Refraction                  |
|-----------------------------------------|-----------------------------------------|
| ![Skybox 1](/img/github/reflection.png) | ![Skybox 2](/img/github/refraction.png) |

Additional render examples are available in the `/img/github` folder.

## Compilation and Execution

### Prerequisites
- `CMake` (version 3.11 or higher)
- `OpenGL` and `GLUT`: For graphical display and real-time visualization.

### Compilation Instructions

1. Clone the repository and create a build directory:
   ```bash
   git clone https://github.com/J-hune/raytracer.git
   cd raytracer
   mkdir build
   cd build
   ```

2. Generate the project using **CMake**:
   ```bash
   cmake ..
   ```

3. Compile the code:
   ```bash
   make
   ```

4. Run the program:
   ```bash
   ./raytracer
   ```

Generated images will be saved in the `/build/renders` folder.

---

## Licence
This project is licensed under the **MIT License**. You are free to use, modify, and redistribute it under the terms of this license. For more details, see the `LICENSE` file.