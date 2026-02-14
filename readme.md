# Raytracer

GPU path tracer built around glTF 2.0 scenes. Blender is the authoring tool; the
application only loads, displays and renders the exported scene.

Geometry traversal, path integration, temporal accumulation and tone mapping
run on the GPU through CUDA and OptiX. The display buffer is shared directly
with OpenGL, without copying rendered pixels through the CPU.

## Scene format

Only `.gltf` and `.glb` are supported. The loader preserves:

- indexed triangle meshes, smooth normals, tangents and UVs;
- shared meshes as instances instead of duplicating geometry;
- metallic-roughness materials, normal maps, UV0/UV1 and texture transforms;
- transmission, IOR, dispersion, volume attenuation and emissive strength;
- perspective cameras and punctual lights;
- `raytracer_aperture` and `raytracer_focus_distance` camera extras.

In Blender, export glTF 2.0 with cameras, punctual lights, materials and custom
properties enabled. Prefer GLB for self-contained scenes.

HDRI settings live in the Blender scene custom properties:

- `raytracer_hdri`: path relative to the exported glTF file;
- `raytracer_hdri_rotation`: rotation in radians;
- `raytracer_hdri_strength`: lighting multiplier;
- `raytracer_exposure`: display exposure in stops.

For colored absorption, add Blender's `glTF Material Output` node to the
material and set a non-zero thickness. This exports standard
`KHR_materials_volume`; the renderer does not use a private material format.

The optional Blender extension in
`tools/blender/raytracer_tools` exposes these settings in the 3D View sidebar.
It can aim the active camera at the selected object, set its focus distance and
export the current scene as GLB. Blender remains the scene editor; the viewer
does not duplicate its placement or material UI.

## Build

Requirements:

- CMake 3.24 or newer;
- a C++20 compiler;
- CUDA Toolkit 12 or newer;
- an NVIDIA Turing-or-newer GPU with an R570-or-newer driver;
- GLFW 3.4, OpenGL, libpng and zlib development packages;
- Git and an internet connection during the first configure.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/bin/raytracer scene.glb
```

The window accumulates one sample per pixel and frame:

- hold right mouse to look around;
- use `WASD` or `ZQSD`, `E` and `C` to move;
- use the wheel to change movement speed and Shift to accelerate;
- press `F` for a 256 spp final render, OptiX denoising and PNG + EXR export;
- press Escape to close.

The final files are written to `renders/<scene>.png` and `.exr`. Offline
renders use the same profiles:

```sh
./build/bin/raytracer scene.glb --output render.png
./build/bin/raytracer scene.glb --profile final --samples 1024 --output render.exr
```

The PNG is denoised and tone mapped. The EXR keeps denoised scene-linear HDR
values. fastgltf `v0.9.0`, OptiX `v9.0.0`, stb and TinyEXR are fetched at
pinned revisions.

## Previous renders

Images produced by the original CPU renderer remain in [`img/github`](img/github).
The CPU renderer, OFF importer and legacy PPM assets have been removed.

## License

MIT. See [`LICENSE`](LICENSE).
