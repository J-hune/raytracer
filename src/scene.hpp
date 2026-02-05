#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

namespace rt {

struct Vec2 {
    float x;
    float y;
};

struct Vec3 {
    float x;
    float y;
    float z;
};

struct Vec4 {
    float x;
    float y;
    float z;
    float w;
};

struct Mat4 {
    std::array<float, 16> values;
};

struct Vertex {
    Vec3 position;
    Vec3 normal;
    Vec4 tangent;
    Vec2 uv;
};

struct TextureRef {
    std::int32_t texture = -1;
    std::uint32_t texCoord = 0;
    Vec2 offset{0.0f, 0.0f};
    Vec2 scale{1.0f, 1.0f};
    float rotation = 0.0f;
    float strength = 1.0f;
};

struct Material {
    std::string name;
    Vec4 baseColor{1.0f, 1.0f, 1.0f, 1.0f};
    Vec3 emissive{0.0f, 0.0f, 0.0f};
    Vec3 attenuationColor{1.0f, 1.0f, 1.0f};
    float metallic = 1.0f;
    float roughness = 1.0f;
    float transmission = 0.0f;
    float ior = 1.5f;
    float thickness = 0.0f;
    float attenuationDistance = std::numeric_limits<float>::infinity();
    float emissiveStrength = 1.0f;
    float dispersion = 0.0f;
    float alphaCutoff = 0.5f;
    std::uint32_t flags = 0;
    TextureRef baseColorTexture;
    TextureRef metallicRoughnessTexture;
    TextureRef normalTexture;
    TextureRef emissiveTexture;
    TextureRef transmissionTexture;
    TextureRef thicknessTexture;
};

struct Geometry {
    std::string name;
    std::vector<Vertex> vertices;
    std::vector<std::uint32_t> indices;
    std::uint32_t material = 0;
};

struct Instance {
    std::string name;
    Mat4 transform;
    std::uint32_t geometry;
};

struct Image {
    std::string name;
    std::vector<std::byte> encoded;
    std::uint32_t mimeType = 0;
};

struct Texture {
    std::string name;
    std::int32_t image = -1;
    std::uint32_t magFilter = 0;
    std::uint32_t minFilter = 0;
    std::uint32_t wrapU = 10497;
    std::uint32_t wrapV = 10497;
};

struct Camera {
    std::string name;
    Mat4 transform;
    float verticalFov = 0.7853982f;
    float aspectRatio = 0.0f;
    float nearPlane = 0.01f;
    float farPlane = std::numeric_limits<float>::infinity();
    float aperture = 0.0f;
    float focusDistance = 10.0f;
};

struct Light {
    std::string name;
    Mat4 transform;
    Vec3 color;
    float intensity;
    float range;
    float innerCone;
    float outerCone;
    std::uint32_t type;
};

struct Scene {
    std::vector<Geometry> geometries;
    std::vector<Instance> instances;
    std::vector<Material> materials;
    std::vector<Image> images;
    std::vector<Texture> textures;
    std::vector<Camera> cameras;
    std::vector<Light> lights;
};

Scene loadScene(const std::filesystem::path& path);

}
