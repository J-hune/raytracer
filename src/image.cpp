#include "image.hpp"

#include <png.h>

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

}
