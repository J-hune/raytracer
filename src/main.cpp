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
    rt::Profile profile = rt::Profile::Preview;
    std::uint32_t samples = 32;
    bool profileSet = false;
    bool samplesSet = false;
    bool denoise = true;
};

Options options(int argc, char** argv) {
    if (argc < 2)
        throw std::runtime_error(
            "usage: raytracer <scene.gltf|scene.glb> [--profile preview|final] "
            "[--output image.png] [--samples count] [--denoise on|off]");

    Options result;
    result.scene = argv[1];
    for (int index = 2; index < argc; index += 2) {
        if (index + 1 >= argc)
            throw std::runtime_error("Missing option value");
        const std::string_view name = argv[index];
        if (name == "--output")
            result.output = argv[index + 1];
        else if (name == "--profile") {
            const std::string_view value = argv[index + 1];
            if (value == "preview")
                result.profile = rt::Profile::Preview;
            else if (value == "final")
                result.profile = rt::Profile::Final;
            else
                throw std::runtime_error("Expected preview or final profile");
            result.profileSet = true;
        } else if (name == "--samples") {
            const std::string_view value = argv[index + 1];
            const auto parsed = std::from_chars(value.data(), value.data() + value.size(),
                                                result.samples);
            if (parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size() ||
                result.samples == 0)
                throw std::runtime_error("Invalid sample count");
            result.samplesSet = true;
        } else if (name == "--denoise") {
            const std::string_view value = argv[index + 1];
            if (value == "on")
                result.denoise = true;
            else if (value == "off")
                result.denoise = false;
            else
                throw std::runtime_error("Expected on or off for denoise");
        } else {
            throw std::runtime_error("Unknown option: " + std::string(name));
        }
    }
    if (result.output && !result.profileSet)
        result.profile = rt::Profile::Final;
    if (!result.samplesSet)
        result.samples = result.profile == rt::Profile::Final ? 256U : 32U;
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
        if (!scene.environment.pixels.empty())
            std::cout << "  HDRI " << scene.environment.width << 'x'
                      << scene.environment.height << ", strength "
                      << scene.environment.strength << '\n';

        const char* profile =
            arguments.profile == rt::Profile::Final ? "final" : "preview";
        std::cout << "  profile " << profile << '\n';
        rt::Renderer renderer(scene, width, height, arguments.profile);
        if (arguments.output) {
            while (renderer.samples() < arguments.samples)
                renderer.render();
            if (arguments.profile == rt::Profile::Final && arguments.denoise)
                renderer.denoise();
            if (arguments.output->extension() == ".exr")
                rt::writeExr(*arguments.output, width, height,
                             renderer.linearPixels());
            else if (arguments.output->extension() == ".png")
                rt::writePng(*arguments.output, width, height, renderer.pixels());
            else
                throw std::runtime_error("Output must use .png or .exr");
            std::cout << "Wrote " << *arguments.output << " at "
                      << renderer.samples() << " spp\n";
            return 0;
        }

        rt::Display display(width, height);
        auto camera = scene.cameras.front();
        bool finalRendering = arguments.profile == rt::Profile::Final;
        bool finalReady = false;
        std::string status = finalRendering ? "final rendering" : "preview";
        const auto finalSamples =
            arguments.samplesSet ? arguments.samples : 256U;
        const auto capture = std::filesystem::path("renders") /
                             arguments.scene.stem();
        auto capturePng = capture;
        auto captureExr = capture;
        capturePng += ".png";
        captureExr += ".exr";
        while (display.open()) {
            if (display.update(camera)) {
                renderer.setProfile(rt::Profile::Preview);
                renderer.setCamera(camera);
                finalRendering = false;
                finalReady = false;
                status = "preview";
            }
            if (display.finalRequested()) {
                renderer.setProfile(rt::Profile::Final);
                finalRendering = true;
                finalReady = false;
                status = "final rendering";
            }

            auto* output = display.map();
            if (finalReady)
                renderer.copyOutput(output);
            else
                renderer.render(output);
            display.present(renderer.samples(), status);

            if (finalRendering && renderer.samples() >= finalSamples) {
                renderer.denoise();
                std::filesystem::create_directories(capture.parent_path());
                rt::writePng(capturePng, width, height, renderer.pixels());
                rt::writeExr(captureExr, width, height, renderer.linearPixels());
                std::cout << "Wrote " << capturePng << " and " << captureExr
                          << " at "
                          << renderer.samples() << " spp\n";
                finalRendering = false;
                finalReady = true;
                status = "final saved";
            }
        }
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
