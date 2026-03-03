#include "display.hpp"
#include "image.hpp"
#include "renderer.hpp"
#include "scene.hpp"

#include <algorithm>
#include <chrono>
#include <charconv>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string_view>
#include <unistd.h>

namespace {

using Clock = std::chrono::steady_clock;

double secondsSince(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

std::string duration(double seconds) {
    const auto value = static_cast<std::uint64_t>(seconds);
    std::ostringstream output;
    if (value >= 3600)
        output << value / 3600 << 'h' << std::setw(2) << std::setfill('0')
               << value / 60 % 60 << 'm';
    else if (value >= 60)
        output << value / 60 << 'm' << std::setw(2) << std::setfill('0')
               << value % 60 << 's';
    else
        output << value << 's';
    return output.str();
}

class Progress {
public:
    explicit Progress(std::uint32_t total)
        : total_(total), interactive_(isatty(STDOUT_FILENO)), start_(Clock::now()),
          update_(start_) {
        show(0);
    }

    void show(std::uint32_t current) {
        const auto now = Clock::now();
        const auto interval = interactive_ ? std::chrono::milliseconds(100)
                                           : std::chrono::seconds(5);
        if (current != total_ && current && now - update_ < interval)
            return;
        update_ = now;

        const double elapsed = std::chrono::duration<double>(now - start_).count();
        const double rate = current / std::max(elapsed, 1e-6);
        const auto percent = 100U * current / total_;
        std::ostringstream line;
        line << "Rendering ";
        if (interactive_) {
            constexpr std::uint32_t width = 24;
            const auto filled = width * current / total_;
            line << '[' << std::string(filled, '#')
                 << std::string(width - filled, '.') << "] ";
        }
        line << std::setw(3) << percent << "% | " << current << '/' << total_
             << " spp | " << duration(elapsed);
        if (current)
            line << " | " << std::fixed << std::setprecision(2) << rate << " spp/s";
        if (current && current != total_)
            line << " | ETA " << duration((total_ - current) / rate);

        if (interactive_)
            std::cout << "\r\033[2K" << line.str() << std::flush;
        else
            std::cout << line.str() << '\n' << std::flush;
    }

    void finish() {
        if (interactive_)
            std::cout << '\n';
    }

private:
    std::uint32_t total_;
    bool interactive_;
    Clock::time_point start_;
    Clock::time_point update_;
};

template <typename Function>
void phase(std::string_view label, Function&& function) {
    std::cout << label << "..." << std::flush;
    const auto start = Clock::now();
    function();
    std::cout << " done in " << duration(secondsSince(start)) << '\n';
}

struct Options {
    std::filesystem::path scene;
    std::optional<std::filesystem::path> output;
    rt::Profile profile = rt::Profile::Final;
    std::uint32_t width = 1280;
    std::uint32_t height = 720;
    std::uint32_t samples = 32;
    bool profileSet = false;
    bool samplesSet = false;
    bool denoise = true;
};

std::uint32_t positiveInteger(std::string_view value, std::string_view name) {
    std::uint32_t result = 0;
    const auto parsed =
        std::from_chars(value.data(), value.data() + value.size(), result);
    if (parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size() ||
        result == 0)
        throw std::runtime_error("Invalid " + std::string(name));
    return result;
}

Options options(int argc, char** argv) {
    if (argc < 2)
        throw std::runtime_error(
            "usage: raytracer <scene.gltf|scene.glb> [--profile preview|final] "
            "[--width pixels] [--height pixels] [--samples count] "
            "[--denoise on|off] [--output image.png]");

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
        } else if (name == "--width") {
            result.width = positiveInteger(argv[index + 1], "width");
        } else if (name == "--height") {
            result.height = positiveInteger(argv[index + 1], "height");
        } else if (name == "--samples") {
            result.samples = positiveInteger(argv[index + 1], "sample count");
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
        const auto arguments = options(argc, argv);
        const auto width = arguments.width;
        const auto height = arguments.height;
        const auto loadStart = Clock::now();
        if (arguments.output)
            std::cout << "Loading " << arguments.scene << "..." << std::flush;
        const auto scene = rt::loadScene(arguments.scene);
        if (arguments.output)
            std::cout << " done in " << duration(secondsSince(loadStart)) << '\n';
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
        std::cout << "  profile " << profile << ", " << width << 'x' << height
                  << '\n';
        const auto initializeStart = Clock::now();
        if (arguments.output)
            std::cout << "Initializing renderer..." << std::flush;
        rt::Renderer renderer(scene, width, height, arguments.profile);
        if (arguments.output)
            std::cout << " done in " << duration(secondsSince(initializeStart)) << '\n';
        if (arguments.output) {
            Progress progress(arguments.samples);
            while (renderer.samples() < arguments.samples) {
                renderer.render();
                progress.show(renderer.samples());
            }
            progress.finish();
            if (arguments.profile == rt::Profile::Final && arguments.denoise)
                phase("Denoising with OptiX", [&] { renderer.denoise(); });
            phase("Writing " + arguments.output->string(), [&] {
                if (arguments.output->extension() == ".exr")
                    rt::writeExr(*arguments.output, width, height,
                                 renderer.linearPixels());
                else if (arguments.output->extension() == ".png")
                    rt::writePng(*arguments.output, width, height, renderer.pixels());
                else
                    throw std::runtime_error("Output must use .png or .exr");
            });
            std::cout << "Done at " << renderer.samples() << " spp\n";
            return 0;
        }

        rt::Display display(width, height);
        auto camera = scene.cameras.front();
        std::string status = "accumulating";
        const auto capture = std::filesystem::path("renders") /
                             arguments.scene.stem();
        auto capturePng = capture;
        auto captureExr = capture;
        capturePng += ".png";
        captureExr += ".exr";
        while (display.open()) {
            if (display.update(camera)) {
                renderer.setCamera(camera);
                status = "accumulating";
            }

            auto* output = display.map();
            renderer.render(output);
            if (display.finalRequested()) {
                renderer.denoise();
                renderer.copyOutput(output);
                std::filesystem::create_directories(capture.parent_path());
                rt::writePng(capturePng, width, height, renderer.pixels());
                rt::writeExr(captureExr, width, height, renderer.linearPixels());
                std::cout << "Wrote " << capturePng << " and " << captureExr
                          << " at " << renderer.samples() << " spp\n";
                status = "captured";
            } else {
                status = "accumulating";
            }
            display.present(renderer.samples(), status);
        }
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
