#define GL_GLEXT_PROTOTYPES
#include "display.hpp"

#include <GLFW/glfw3.h>
#include <cuda_gl_interop.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace rt {
namespace {

void check(cudaError_t result) {
    if (result != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(result));
}

Vec3 operator+(Vec3 a, Vec3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 operator-(Vec3 a, Vec3 b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 operator-(Vec3 value) {
    return {-value.x, -value.y, -value.z};
}

Vec3 operator*(Vec3 value, float scale) {
    return {value.x * scale, value.y * scale, value.z * scale};
}

float dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}

Vec3 normalized(Vec3 value) {
    const float inverse = 1.0f / std::sqrt(dot(value, value));
    return value * inverse;
}

Vec3 rotate(Vec3 value, Vec3 axis, float angle) {
    const float cosine = std::cos(angle);
    const float sine = std::sin(angle);
    return value * cosine + cross(axis, value) * sine +
           axis * (dot(axis, value) * (1.0f - cosine));
}

Vec3 column(const Mat4& matrix, std::size_t index) {
    return {matrix.values[index], matrix.values[index + 1],
            matrix.values[index + 2]};
}

void column(Mat4& matrix, std::size_t index, Vec3 value) {
    matrix.values[index] = value.x;
    matrix.values[index + 1] = value.y;
    matrix.values[index + 2] = value.z;
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
    lastTime_ = glfwGetTime();
    glfwSetWindowUserPointer(window_, this);
    glfwSetScrollCallback(window_, [](GLFWwindow* window, double, double offset) {
        auto* display = static_cast<Display*>(glfwGetWindowUserPointer(window));
        display->speed_ = std::clamp(
            display->speed_ * std::pow(1.25f, static_cast<float>(offset)),
            0.1f, 50.0f);
    });

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

bool Display::update(Camera& camera) {
    glfwPollEvents();
    const double now = glfwGetTime();
    const float elapsed =
        static_cast<float>(std::min(now - lastTime_, 0.1));
    lastTime_ = now;
    const auto down = [&](int key) {
        return glfwGetKey(window_, key) == GLFW_PRESS;
    };

    const bool finalDown = down(GLFW_KEY_F);
    finalRequested_ |= finalDown && !finalDown_;
    finalDown_ = finalDown;

    auto right = normalized(column(camera.transform, 0));
    auto up = normalized(column(camera.transform, 4));
    auto forward = -normalized(column(camera.transform, 8));
    bool changed = false;

    const bool look = glfwGetMouseButton(window_, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
    if (look && !looking_) {
        looking_ = true;
        glfwGetCursorPos(window_, &lastMouseX_, &lastMouseY_);
        glfwSetInputMode(window_, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    } else if (!look && looking_) {
        looking_ = false;
        glfwSetInputMode(window_, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    }
    if (looking_) {
        double x;
        double y;
        glfwGetCursorPos(window_, &x, &y);
        const float yaw = static_cast<float>(lastMouseX_ - x) * 0.002f;
        const float pitch = static_cast<float>(lastMouseY_ - y) * 0.002f;
        lastMouseX_ = x;
        lastMouseY_ = y;
        if (yaw != 0.0f || pitch != 0.0f) {
            constexpr Vec3 worldUp{0.0f, 1.0f, 0.0f};
            forward = normalized(rotate(forward, worldUp, yaw));
            right = normalized(cross(forward, worldUp));
            const auto pitched = normalized(rotate(forward, right, pitch));
            if (std::abs(dot(pitched, worldUp)) < 0.995f)
                forward = pitched;
            right = normalized(cross(forward, worldUp));
            up = normalized(cross(right, forward));
            column(camera.transform, 0, right);
            column(camera.transform, 4, up);
            column(camera.transform, 8, -forward);
            changed = true;
        }
    }

    Vec3 movement{};
    if (down(GLFW_KEY_W) || down(GLFW_KEY_Z))
        movement = movement + forward;
    if (down(GLFW_KEY_S))
        movement = movement - forward;
    if (down(GLFW_KEY_A) || down(GLFW_KEY_Q))
        movement = movement - right;
    if (down(GLFW_KEY_D))
        movement = movement + right;
    if (down(GLFW_KEY_E))
        movement = movement + up;
    if (down(GLFW_KEY_C))
        movement = movement - up;
    if (dot(movement, movement) > 0.0f) {
        const float boost = down(GLFW_KEY_LEFT_SHIFT) ? 4.0f : 1.0f;
        const auto position = column(camera.transform, 12) +
            normalized(movement) * (speed_ * boost * elapsed);
        column(camera.transform, 12, position);
        changed = true;
    }
    return changed;
}

bool Display::finalRequested() {
    return std::exchange(finalRequested_, false);
}

void* Display::map() {
    check(cudaGraphicsMapResources(1, &resource_));
    void* data = nullptr;
    std::size_t bytes = 0;
    check(cudaGraphicsResourceGetMappedPointer(&data, &bytes, resource_));
    return data;
}

void Display::present(std::uint32_t samples, std::string_view status) {
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

    const auto title = "Raytracer | " + std::string(status) + " | " +
        std::to_string(samples) +
        " spp | RMB look  WASD/ZQSD move  wheel speed  F final";
    glfwSetWindowTitle(window_, title.c_str());
    glfwSwapBuffers(window_);
    if (glfwGetKey(window_, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window_, GLFW_TRUE);
}

}
