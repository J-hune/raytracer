#define GL_GLEXT_PROTOTYPES
#include "display.hpp"

#include <GLFW/glfw3.h>
#include <cuda_gl_interop.h>

#include <stdexcept>
#include <string>

namespace rt {
namespace {

void check(cudaError_t result) {
    if (result != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(result));
}

}

Display::Display(std::uint32_t width, std::uint32_t height)
    : width_(width), height_(height) {
    if (!glfwInit())
        throw std::runtime_error("Unable to initialize GLFW");
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);
    window_ = glfwCreateWindow(static_cast<int>(width), static_cast<int>(height),
                               "Raytracer", nullptr, nullptr);
    if (!window_)
        throw std::runtime_error("Unable to create the render window");
    glfwMakeContextCurrent(window_);
    glfwSwapInterval(0);

    glGenBuffers(1, &pixelBuffer_);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pixelBuffer_);
    glBufferData(GL_PIXEL_UNPACK_BUFFER, width * height * 4, nullptr, GL_STREAM_DRAW);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
    check(cudaGraphicsGLRegisterBuffer(&resource_, pixelBuffer_,
                                       cudaGraphicsRegisterFlagsWriteDiscard));

    glGenTextures(1, &texture_);
    glBindTexture(GL_TEXTURE_2D, texture_);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, static_cast<int>(width),
                 static_cast<int>(height), 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
}

Display::~Display() {
    if (resource_)
        cudaGraphicsUnregisterResource(resource_);
    if (texture_)
        glDeleteTextures(1, &texture_);
    if (pixelBuffer_)
        glDeleteBuffers(1, &pixelBuffer_);
    if (window_)
        glfwDestroyWindow(window_);
    glfwTerminate();
}

bool Display::open() const {
    return window_ && !glfwWindowShouldClose(window_);
}

void* Display::map() {
    check(cudaGraphicsMapResources(1, &resource_));
    void* data = nullptr;
    std::size_t bytes = 0;
    check(cudaGraphicsResourceGetMappedPointer(&data, &bytes, resource_));
    return data;
}

void Display::present(std::uint32_t samples) {
    check(cudaGraphicsUnmapResources(1, &resource_));
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pixelBuffer_);
    glBindTexture(GL_TEXTURE_2D, texture_);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, static_cast<int>(width_),
                    static_cast<int>(height_), GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);

    glClear(GL_COLOR_BUFFER_BIT);
    glEnable(GL_TEXTURE_2D);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f); glVertex2f(-1.0f, -1.0f);
    glTexCoord2f(1.0f, 0.0f); glVertex2f(1.0f, -1.0f);
    glTexCoord2f(1.0f, 1.0f); glVertex2f(1.0f, 1.0f);
    glTexCoord2f(0.0f, 1.0f); glVertex2f(-1.0f, 1.0f);
    glEnd();

    const auto title = "Raytracer - " + std::to_string(samples) + " spp";
    glfwSetWindowTitle(window_, title.c_str());
    glfwSwapBuffers(window_);
    glfwPollEvents();
    if (glfwGetKey(window_, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window_, GLFW_TRUE);
}

}
