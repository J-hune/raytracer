#include "image.hpp"

#include <png.h>
#include <tinyexr.h>

#include <cstring>
#include <stdexcept>
#include <vector>

namespace rt {

void writePng(const std::filesystem::path& path, std::uint32_t width,
              std::uint32_t height, std::span<const std::uint8_t> pixels) {
    const auto stride = static_cast<std::size_t>(width) * 4;
    if (pixels.size() != stride * height)
        throw std::runtime_error("Invalid RGBA image size");

    std::vector<std::uint8_t> flipped(pixels.size());
    for (std::uint32_t row = 0; row < height; ++row)
        std::memcpy(flipped.data() + row * stride,
                    pixels.data() + (height - row - 1) * stride, stride);

    png_image image{};
    image.version = PNG_IMAGE_VERSION;
    image.width = width;
    image.height = height;
    image.format = PNG_FORMAT_RGBA;
    if (!png_image_write_to_file(&image, path.c_str(), false, flipped.data(), 0, nullptr))
        throw std::runtime_error("Unable to write PNG: " + std::string(image.message));
}

void writeExr(const std::filesystem::path& path, std::uint32_t width,
              std::uint32_t height, std::span<const float> pixels) {
    const auto stride = static_cast<std::size_t>(width) * 4;
    if (pixels.size() != stride * height)
        throw std::runtime_error("Invalid RGBA image size");

    std::vector<float> flipped(pixels.size());
    for (std::uint32_t row = 0; row < height; ++row)
        std::memcpy(flipped.data() + row * stride,
                    pixels.data() + (height - row - 1) * stride,
                    stride * sizeof(float));

    const char* message = nullptr;
    const auto status = SaveEXR(flipped.data(), static_cast<int>(width),
                                static_cast<int>(height), 4, true,
                                path.c_str(), &message);
    if (status != TINYEXR_SUCCESS) {
        const std::string error = message ? message : "unknown error";
        FreeEXRErrorMessage(message);
        throw std::runtime_error("Unable to write EXR: " + error);
    }
}

}
