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
- metallic-roughness materials and texture transforms;
- transmission, IOR, dispersion, volume attenuation and emissive strength;
- perspective cameras and punctual lights;
- `raytracer_aperture` and `raytracer_focus_distance` camera extras.

In Blender, export glTF 2.0 with cameras, punctual lights, materials and custom
properties enabled. Prefer GLB for self-contained scenes.

## Build

Requirements:

- CMake 3.24 or newer;
- a C++20 compiler;
- CUDA Toolkit 12 or newer;
- an NVIDIA Turing-or-newer GPU with an R570-or-newer driver;
- GLFW 3.4 and OpenGL development packages;
- Git and an internet connection during the first configure.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/bin/raytracer scene.glb
```

The window accumulates one sample per pixel and frame. Press `Escape` to close
it. fastgltf `v0.9.0` and OptiX `v9.0.0` are fetched at pinned tags.

## Previous renders

Images produced by the original CPU renderer remain in [`img/github`](img/github).
The CPU renderer, OFF importer and legacy PPM assets have been removed.

## License

MIT. See [`LICENSE`](LICENSE).
