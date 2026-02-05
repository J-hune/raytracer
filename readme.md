# Raytracer

GPU path tracer built around glTF 2.0 scenes. Blender is the authoring tool; the
application only loads, displays and renders the exported scene.

The renderer is being migrated to CUDA/OptiX. The current executable validates
the new scene pipeline and reports the resources prepared for GPU upload.

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
- Git and an internet connection during the first configure.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/bin/raytracer scene.glb
```

fastgltf is fetched at the pinned `v0.9.0` tag. CUDA and OptiX become required
in the next renderer milestone.

## Previous renders

Images produced by the original CPU renderer remain in [`img/github`](img/github).
The CPU renderer, OFF importer and legacy PPM assets have been removed.

## License

MIT. See [`LICENSE`](LICENSE).
