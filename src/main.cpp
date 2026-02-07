#include "display.hpp"
#include "renderer.hpp"
#include "scene.hpp"

#include <exception>
#include <filesystem>
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: raytracer <scene.gltf|scene.glb>\n";
        return 2;
    }

    try {
        constexpr std::uint32_t width = 1280;
        constexpr std::uint32_t height = 720;
        const auto scene = rt::loadScene(std::filesystem::path(argv[1]));
        std::cout << "Loaded " << argv[1] << '\n'
                  << "  " << scene.geometries.size() << " geometries, "
                  << scene.instances.size() << " instances\n"
                  << "  " << scene.materials.size() << " materials, "
                  << scene.textures.size() << " textures, "
                  << scene.images.size() << " images\n"
                  << "  " << scene.cameras.size() << " cameras, "
                  << scene.lights.size() << " lights\n";
        if (!scene.cameras.empty())
            std::cout << "  camera aperture " << scene.cameras.front().aperture
                      << ", focus " << scene.cameras.front().focusDistance << " m\n";

        rt::Display display(width, height);
        rt::Renderer renderer(scene, width, height);
        while (display.open()) {
            renderer.render(display.map());
            display.present(renderer.samples());
        }
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
