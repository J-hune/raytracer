#include "scene.hpp"

#include <fastgltf/core.hpp>
#include <fastgltf/tools.hpp>
#include <simdjson.h>
#include <stb_image.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <ranges>
#include <stdexcept>
#include <string_view>

namespace rt {
namespace {

constexpr std::uint32_t DoubleSided = 1U << 0U;
constexpr std::uint32_t Unlit = 1U << 1U;
constexpr std::uint32_t AlphaMask = 1U << 2U;
constexpr std::uint32_t AlphaBlend = 1U << 3U;

struct CameraExtras {
    float aperture = 0.0f;
    float focusDistance = 10.0f;
};

struct Extras {
    std::vector<CameraExtras> cameras;
    std::string environment;
    float rotation = 0.0f;
    float strength = 1.0f;
    float exposure = 0.0f;
};

void parseExtras(simdjson::dom::object* extras, std::size_t index,
                 fastgltf::Category category, void* userData) {
    auto& state = *static_cast<Extras*>(userData);
    if (category == fastgltf::Category::Cameras) {
        state.cameras.resize(std::max(state.cameras.size(), index + 1));
        auto& camera = state.cameras[index];
        double value = 0.0;
        if ((*extras)["raytracer_aperture"].get_double().get(value) == simdjson::SUCCESS)
            camera.aperture = static_cast<float>(value);
        if ((*extras)["raytracer_focus_distance"].get_double().get(value) ==
            simdjson::SUCCESS)
            camera.focusDistance = static_cast<float>(value);
        return;
    }
    if (category != fastgltf::Category::Scenes)
        return;

    std::string_view path;
    if ((*extras)["raytracer_hdri"].get_string().get(path) == simdjson::SUCCESS)
        state.environment = path;
    double value = 0.0;
    if ((*extras)["raytracer_hdri_rotation"].get_double().get(value) ==
        simdjson::SUCCESS)
        state.rotation = static_cast<float>(value);
    if ((*extras)["raytracer_hdri_strength"].get_double().get(value) ==
        simdjson::SUCCESS)
        state.strength = static_cast<float>(value);
    if ((*extras)["raytracer_exposure"].get_double().get(value) == simdjson::SUCCESS)
        state.exposure = static_cast<float>(value);
}

Vec2 vec2(const fastgltf::math::nvec2& value) {
    return {value[0], value[1]};
}

Vec3 vec3(const fastgltf::math::nvec3& value) {
    return {value[0], value[1], value[2]};
}

Vec4 vec4(const fastgltf::math::nvec4& value) {
    return {value[0], value[1], value[2], value[3]};
}

Mat4 mat4(const fastgltf::math::fmat4x4& value) {
    Mat4 result{};
    for (std::size_t column = 0; column < 4; ++column)
        for (std::size_t row = 0; row < 4; ++row)
            result.values[column * 4 + row] = value[column][row];
    return result;
}

std::runtime_error error(std::string_view message, std::string_view name = {}) {
    auto text = std::string(message);
    if (!name.empty())
        text += ": " + std::string(name);
    return std::runtime_error(text);
}

TextureRef textureRef(const fastgltf::TextureInfo& source) {
    TextureRef result;
    result.texture = static_cast<std::int32_t>(source.textureIndex);
    result.texCoord = static_cast<std::uint32_t>(source.texCoordIndex);
    if (source.transform) {
        result.offset = vec2(source.transform->uvOffset);
        result.scale = vec2(source.transform->uvScale);
        result.rotation = source.transform->rotation;
        if (source.transform->texCoordIndex)
            result.texCoord = static_cast<std::uint32_t>(*source.transform->texCoordIndex);
    }
    return result;
}

template<typename TextureInfo>
TextureRef textureRef(const std::optional<TextureInfo>& source) {
    if (!source)
        return {};
    auto result = textureRef(static_cast<const fastgltf::TextureInfo&>(*source));
    if constexpr (requires { source->scale; })
        result.strength = source->scale;
    return result;
}

Material material(const fastgltf::Material& source) {
    Material result;
    result.name = source.name;
    result.baseColor = vec4(source.pbrData.baseColorFactor);
    result.emissive = vec3(source.emissiveFactor);
    result.metallic = source.pbrData.metallicFactor;
    result.roughness = source.pbrData.roughnessFactor;
    result.ior = source.ior;
    result.dispersion = source.dispersion;
    result.emissiveStrength = source.emissiveStrength;
    result.alphaCutoff = source.alphaCutoff;
    result.flags = (source.doubleSided ? DoubleSided : 0U) |
                   (source.unlit ? Unlit : 0U) |
                   (source.alphaMode == fastgltf::AlphaMode::Mask ? AlphaMask : 0U) |
                   (source.alphaMode == fastgltf::AlphaMode::Blend ? AlphaBlend : 0U);

    if (source.transmission) {
        result.transmission = source.transmission->transmissionFactor;
        result.transmissionTexture = textureRef(source.transmission->transmissionTexture);
    }
    if (source.volume) {
        result.thickness = source.volume->thicknessFactor;
        result.attenuationColor = vec3(source.volume->attenuationColor);
        result.attenuationDistance = source.volume->attenuationDistance;
        result.thicknessTexture = textureRef(source.volume->thicknessTexture);
    }

    result.baseColorTexture = textureRef(source.pbrData.baseColorTexture);
    result.metallicRoughnessTexture = textureRef(source.pbrData.metallicRoughnessTexture);
    result.normalTexture = textureRef(source.normalTexture);
    result.emissiveTexture = textureRef(source.emissiveTexture);
    return result;
}

Geometry geometry(const fastgltf::Asset& asset, const fastgltf::Mesh& mesh,
                  const fastgltf::Primitive& primitive, std::size_t primitiveIndex) {
    if (primitive.type != fastgltf::PrimitiveType::Triangles)
        throw error("Only triangle primitives are supported", mesh.name);

    const auto positions = primitive.findAttribute("POSITION");
    if (positions == primitive.attributes.end() || !primitive.indicesAccessor)
        throw error("Malformed mesh primitive", mesh.name);

    const auto& positionAccessor = asset.accessors[positions->accessorIndex];
    Geometry result;
    result.name = mesh.name.empty() ? "mesh." + std::to_string(primitiveIndex)
                                    : std::string(mesh.name) + "." + std::to_string(primitiveIndex);
    result.vertices.resize(positionAccessor.count);

    fastgltf::iterateAccessorWithIndex<fastgltf::math::fvec3>(
        asset, positionAccessor, [&](const auto& value, std::size_t index) {
            result.vertices[index].position = vec3(value);
        });

    const auto normals = primitive.findAttribute("NORMAL");
    if (normals != primitive.attributes.end()) {
        fastgltf::iterateAccessorWithIndex<fastgltf::math::fvec3>(
            asset, asset.accessors[normals->accessorIndex], [&](const auto& value, std::size_t index) {
                result.vertices.at(index).normal = vec3(value);
            });
    }

    const auto tangents = primitive.findAttribute("TANGENT");
    if (tangents != primitive.attributes.end()) {
        fastgltf::iterateAccessorWithIndex<fastgltf::math::fvec4>(
            asset, asset.accessors[tangents->accessorIndex], [&](const auto& value, std::size_t index) {
                result.vertices.at(index).tangent = {value[0], value[1], value[2], value[3]};
            });
    }

    const auto texcoords = primitive.findAttribute("TEXCOORD_0");
    if (texcoords != primitive.attributes.end()) {
        fastgltf::iterateAccessorWithIndex<fastgltf::math::fvec2>(
            asset, asset.accessors[texcoords->accessorIndex], [&](const auto& value, std::size_t index) {
                result.vertices.at(index).uv = {value[0], value[1]};
            });
    }
    const auto texcoords1 = primitive.findAttribute("TEXCOORD_1");
    if (texcoords1 != primitive.attributes.end()) {
        fastgltf::iterateAccessorWithIndex<fastgltf::math::fvec2>(
            asset, asset.accessors[texcoords1->accessorIndex],
            [&](const auto& value, std::size_t index) {
                result.vertices.at(index).uv1 = {value[0], value[1]};
            });
    }

    const auto& indexAccessor = asset.accessors[*primitive.indicesAccessor];
    result.indices.resize(indexAccessor.count);
    fastgltf::copyFromAccessor<std::uint32_t>(asset, indexAccessor, result.indices.data());
    if (result.indices.size() % 3 != 0 ||
        std::ranges::any_of(result.indices, [&](auto index) { return index >= result.vertices.size(); }))
        throw error("Invalid triangle indices", result.name);

    result.material = primitive.materialIndex
        ? static_cast<std::uint32_t>(*primitive.materialIndex + 1)
        : 0U;
    return result;
}

Image image(const fastgltf::Asset& asset, const fastgltf::Image& source) {
    Image result;
    result.name = source.name;
    const auto copy = [&](const std::byte* bytes, std::size_t size, fastgltf::MimeType mime) {
        result.encoded.assign(bytes, bytes + size);
        result.mimeType = static_cast<std::uint32_t>(mime);
    };

    std::visit(fastgltf::visitor{
        [&](const fastgltf::sources::Array& data) {
            copy(data.bytes.data(), data.bytes.size_bytes(), data.mimeType);
        },
        [&](const fastgltf::sources::Vector& data) {
            copy(data.bytes.data(), data.bytes.size(), data.mimeType);
        },
        [&](const fastgltf::sources::ByteView& data) {
            copy(data.bytes.data(), data.bytes.size_bytes(), data.mimeType);
        },
        [&](const fastgltf::sources::BufferView& data) {
            const auto bytes = fastgltf::DefaultBufferDataAdapter{}(asset, data.bufferViewIndex);
            copy(bytes.data(), bytes.size_bytes(), data.mimeType);
        },
        [&](const auto&) {
            throw error("Image data was not loaded", source.name);
        }
    }, source.data);
    return result;
}

Texture texture(const fastgltf::Asset& asset, const fastgltf::Texture& source) {
    Texture result;
    result.name = source.name;
    if (!source.imageIndex)
        throw error("Only core glTF image sources are supported", source.name);
    result.image = static_cast<std::int32_t>(*source.imageIndex);
    if (source.samplerIndex) {
        const auto& sampler = asset.samplers[*source.samplerIndex];
        result.magFilter = sampler.magFilter ? fastgltf::to_underlying(*sampler.magFilter) : 0U;
        result.minFilter = sampler.minFilter ? fastgltf::to_underlying(*sampler.minFilter) : 0U;
        result.wrapU = fastgltf::to_underlying(sampler.wrapS);
        result.wrapV = fastgltf::to_underlying(sampler.wrapT);
    }
    return result;
}

Camera camera(const fastgltf::Camera& source, const fastgltf::math::fmat4x4& transform,
              const CameraExtras& extras) {
    const auto* perspective = std::get_if<fastgltf::Camera::Perspective>(&source.camera);
    if (!perspective)
        throw error("Orthographic cameras are not supported", source.name);

    Camera result;
    result.name = source.name;
    result.transform = mat4(transform);
    result.verticalFov = perspective->yfov;
    result.aspectRatio = perspective->aspectRatio.value_or(0.0f);
    result.nearPlane = perspective->znear;
    result.farPlane = perspective->zfar.value_or(std::numeric_limits<float>::infinity());
    result.aperture = extras.aperture;
    result.focusDistance = extras.focusDistance;
    return result;
}

Light light(const fastgltf::Light& source, const fastgltf::math::fmat4x4& transform) {
    return {
        std::string(source.name),
        mat4(transform),
        vec3(source.color),
        source.intensity,
        source.range.value_or(std::numeric_limits<float>::infinity()),
        source.innerConeAngle.value_or(0.0f),
        source.outerConeAngle.value_or(0.7853982f),
        static_cast<std::uint32_t>(source.type)
    };
}

Environment environment(const std::filesystem::path& scenePath, const Extras& extras) {
    Environment result;
    result.rotation = extras.rotation;
    result.strength = extras.strength;
    result.exposure = extras.exposure;
    if (extras.environment.empty())
        return result;

    const auto path = scenePath.parent_path() / extras.environment;
    int width;
    int height;
    int channels;
    float* source = stbi_loadf(path.c_str(), &width, &height, &channels, 4);
    if (!source)
        throw error(stbi_failure_reason(), path.string());
    result.width = static_cast<std::uint32_t>(width);
    result.height = static_cast<std::uint32_t>(height);
    result.pixels.resize(static_cast<std::size_t>(width) * height);
    std::memcpy(result.pixels.data(), source,
                result.pixels.size() * sizeof(result.pixels.front()));
    stbi_image_free(source);
    return result;
}

}

Scene loadScene(const std::filesystem::path& path) {
    if (path.extension() != ".gltf" && path.extension() != ".glb")
        throw error("Expected a glTF 2.0 .gltf or .glb scene", path.string());

    constexpr auto extensions =
        fastgltf::Extensions::KHR_lights_punctual |
        fastgltf::Extensions::KHR_materials_dispersion |
        fastgltf::Extensions::KHR_materials_emissive_strength |
        fastgltf::Extensions::KHR_materials_ior |
        fastgltf::Extensions::KHR_materials_transmission |
        fastgltf::Extensions::KHR_materials_unlit |
        fastgltf::Extensions::KHR_materials_volume |
        fastgltf::Extensions::KHR_mesh_quantization |
        fastgltf::Extensions::KHR_texture_transform;
    constexpr auto options =
        fastgltf::Options::LoadExternalBuffers |
        fastgltf::Options::LoadExternalImages |
        fastgltf::Options::GenerateMeshIndices;

    auto file = fastgltf::MappedGltfFile::FromPath(path);
    if (!file)
        throw error(fastgltf::getErrorMessage(file.error()), path.string());

    Extras extras;
    fastgltf::Parser parser(extensions);
    parser.setExtrasParseCallback(parseExtras);
    parser.setUserPointer(&extras);
    auto loaded = parser.loadGltf(file.get(), path.parent_path(), options);
    if (!loaded)
        throw error(fastgltf::getErrorMessage(loaded.error()), path.string());
    auto asset = std::move(loaded.get());
    if (const auto validation = fastgltf::validate(asset); validation != fastgltf::Error::None)
        throw error(fastgltf::getErrorMessage(validation), path.string());

    Scene scene;
    scene.environment = environment(path, extras);
    scene.materials.reserve(asset.materials.size() + 1);
    scene.materials.emplace_back();
    for (const auto& source : asset.materials)
        scene.materials.push_back(material(source));
    for (const auto& source : asset.images)
        scene.images.push_back(image(asset, source));
    for (const auto& source : asset.textures)
        scene.textures.push_back(texture(asset, source));

    std::vector<std::vector<std::uint32_t>> meshGeometries(asset.meshes.size());
    for (std::size_t meshIndex = 0; meshIndex < asset.meshes.size(); ++meshIndex) {
        const auto& mesh = asset.meshes[meshIndex];
        for (std::size_t primitiveIndex = 0; primitiveIndex < mesh.primitives.size(); ++primitiveIndex) {
            meshGeometries[meshIndex].push_back(static_cast<std::uint32_t>(scene.geometries.size()));
            scene.geometries.push_back(geometry(asset, mesh, mesh.primitives[primitiveIndex], primitiveIndex));
        }
    }

    if (asset.scenes.empty())
        throw error("glTF contains no scene", path.string());
    const auto sceneIndex = asset.defaultScene.value_or(0);
    fastgltf::iterateSceneNodes(asset, sceneIndex, fastgltf::math::fmat4x4{},
        [&](const fastgltf::Node& node, const fastgltf::math::fmat4x4& transform) {
            if (node.meshIndex) {
                for (const auto geometryIndex : meshGeometries[*node.meshIndex])
                    scene.instances.push_back({std::string(node.name), mat4(transform), geometryIndex});
            }
            if (node.cameraIndex) {
                const auto index = *node.cameraIndex;
                const auto values =
                    index < extras.cameras.size() ? extras.cameras[index] : CameraExtras{};
                scene.cameras.push_back(camera(asset.cameras[index], transform, values));
            }
            if (node.lightIndex)
                scene.lights.push_back(light(asset.lights[*node.lightIndex], transform));
        });
    return scene;
}

}
