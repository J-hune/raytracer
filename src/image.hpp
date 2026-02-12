#pragma once

#include <cstdint>
#include <filesystem>
#include <span>

namespace rt {

void writePng(const std::filesystem::path& path, std::uint32_t width,
              std::uint32_t height, std::span<const std::uint8_t> pixels);
void writeExr(const std::filesystem::path& path, std::uint32_t width,
              std::uint32_t height, std::span<const float> pixels);

}
