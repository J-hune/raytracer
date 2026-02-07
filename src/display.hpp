#pragma once

#include <cstdint>

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
    void* map();
    void present(std::uint32_t samples);

private:
    GLFWwindow* window_ = nullptr;
    cudaGraphicsResource* resource_ = nullptr;
    std::uint32_t width_;
    std::uint32_t height_;
    unsigned int pixelBuffer_ = 0;
    unsigned int texture_ = 0;
};

}
