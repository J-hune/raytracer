#pragma once

#include "scene.hpp"

#include <cstdint>
#include <memory>

namespace rt {

class Renderer {
public:
    Renderer(const Scene& scene, std::uint32_t width, std::uint32_t height);
    ~Renderer();

    Renderer(const Renderer&) = delete;
    Renderer& operator=(const Renderer&) = delete;

    void render(void* output);
    std::uint32_t samples() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}
