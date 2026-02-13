#pragma once

#include "scene.hpp"

#include <cstdint>
#include <string_view>

struct GLFWwindow;
struct cudaGraphicsResource;

namespace rt {

class Display {
public:
    Display(std::uint32_t width, std::uint32_t height);
    ~Display();

    Display(const Display&) = delete;
    Display& operator=(const Display&) = delete;

    bool open() const;
    bool update(Camera& camera);
    bool finalRequested();
    void* map();
    void present(std::uint32_t samples, std::string_view status);

private:
    GLFWwindow* window_ = nullptr;
    cudaGraphicsResource* resource_ = nullptr;
    std::uint32_t width_;
    std::uint32_t height_;
    unsigned int pixelBuffer_ = 0;
    unsigned int texture_ = 0;
    double lastTime_ = 0.0;
    double lastMouseX_ = 0.0;
    double lastMouseY_ = 0.0;
    float speed_ = 3.0f;
    bool looking_ = false;
    bool finalDown_ = false;
    bool finalRequested_ = false;
};

}
