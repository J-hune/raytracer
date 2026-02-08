#include "display.hpp"
#include "image.hpp"
#include "renderer.hpp"
#include "scene.hpp"

#include <charconv>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string_view>

namespace {

struct Options {
    std::filesystem::path scene;
    std::optional<std::filesystem::path> output;
    std::uint32_t samples = 256;
};

Options options(int argc, char** argv) {
    if (argc < 2)
        throw std::runtime_error(
            "usage: raytracer <scene.gltf|scene.glb> [--output image.png] [--samples count]");

    Options result;
    result.scene = argv[1];
    for (int index = 2; index < argc; index += 2) {
        if (index + 1 >= argc)
            throw std::runtime_error("Missing option value");
        const std::string_view name = argv[index];
        if (name == "--output")
            result.output = argv[index + 1];
        else if (name == "--samples") {
            const std::string_view value = argv[index + 1];
            const auto parsed = std::from_chars(value.data(), value.data() + value.size(),
                                                result.samples);
            if (parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size() ||
                result.samples == 0)
                throw std::runtime_error("Invalid sample count");
        } else {
            throw std::runtime_error("Unknown option: " + std::string(name));
        }
    }
    return result;
}

}

int main(int argc, char** argv) {
    try {
        constexpr std::uint32_t width = 1280;
        constexpr std::uint32_t height = 720;
        const auto arguments = options(argc, argv);
        const auto scene = rt::loadScene(arguments.scene);
        std::cout << "Loaded " << arguments.scene << '\n'
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

        rt::Renderer renderer(scene, width, height);
        if (arguments.output) {
            while (renderer.samples() < arguments.samples)
                renderer.render();
            rt::writePng(*arguments.output, width, height, renderer.pixels());
            std::cout << "Wrote " << *arguments.output << " at "
                      << renderer.samples() << " spp\n";
            return 0;
        }

        rt::Display display(width, height);
        while (display.open()) {
            renderer.render(display.map());
            display.present(renderer.samples());
        }
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
