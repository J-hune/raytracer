#pragma once

#include "scene.hpp"

#include <cstdint>
#include <memory>
#include <vector>

namespace rt {

enum class Profile {
    Preview,
    Final
};

class Renderer {
public:
    Renderer(const Scene& scene, std::uint32_t width, std::uint32_t height,
             Profile profile = Profile::Preview);
    ~Renderer();

    Renderer(const Renderer&) = delete;
    Renderer& operator=(const Renderer&) = delete;

    void render(void* output = nullptr);
    void denoise();
    void setCamera(const Camera& camera);
    void setProfile(Profile profile);
    void copyOutput(void* output) const;
    std::vector<std::uint8_t> pixels() const;
    std::vector<float> linearPixels() const;
    std::uint32_t samples() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}
