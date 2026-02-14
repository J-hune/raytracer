#include "renderer.hpp"

#include "gpu_shared.hpp"

#include <cuda_runtime.h>
#include <optix_function_table_definition.h>
#include <optix_stack_size.h>
#include <optix_stubs.h>
#include <stb_image.h>

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

extern "C" const unsigned char deviceCode[];
extern "C" const unsigned long deviceCodeSize;

namespace rt {
namespace {

void cudaCheck(cudaError_t result) {
    if (result != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(result));
}

void optixCheck(OptixResult result, const char* operation) {
    if (result != OPTIX_SUCCESS)
        throw std::runtime_error(std::string(operation) + ": " +
                                 optixGetErrorName(result));
}

class Buffer {
public:
    Buffer() = default;
    explicit Buffer(std::size_t bytes) { resize(bytes); }
    ~Buffer() { reset(); }

    Buffer(Buffer&& other) noexcept : data_(std::exchange(other.data_, nullptr)),
                                      size_(std::exchange(other.size_, 0)) {}
    Buffer& operator=(Buffer&& other) noexcept {
        if (this != &other) {
            reset();
            data_ = std::exchange(other.data_, nullptr);
            size_ = std::exchange(other.size_, 0);
        }
        return *this;
    }

    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;

    void resize(std::size_t bytes) {
        if (bytes == size_)
            return;
        reset();
        if (bytes)
            cudaCheck(cudaMalloc(&data_, bytes));
        size_ = bytes;
    }

    void upload(const void* source, std::size_t bytes) {
        resize(bytes);
        if (bytes)
            cudaCheck(cudaMemcpy(data_, source, bytes, cudaMemcpyHostToDevice));
    }

    CUdeviceptr device() const {
        return reinterpret_cast<CUdeviceptr>(data_);
    }

    void* data() const { return data_; }
    std::size_t size() const { return size_; }

private:
    void reset() {
        if (data_)
            cudaFree(data_);
        data_ = nullptr;
        size_ = 0;
    }

    void* data_ = nullptr;
    std::size_t size_ = 0;
};

template<typename T>
struct alignas(OPTIX_SBT_RECORD_ALIGNMENT) Record {
    std::array<char, OPTIX_SBT_RECORD_HEADER_SIZE> header;
    T data;
};

struct Empty {};

float3 normalized(float x, float y, float z) {
    const float inverse = 1.0f / std::sqrt(x * x + y * y + z * z);
    return make_float3(x * inverse, y * inverse, z * inverse);
}

GpuTextureRef gpuTextureRef(const TextureRef& source) {
    if (source.texture >= 0 && source.texCoord > 1)
        throw std::runtime_error("Only glTF TEXCOORD_0 and TEXCOORD_1 are supported");
    return {
        source.texture,
        source.texCoord,
        make_float2(source.offset.x, source.offset.y),
        make_float2(source.scale.x, source.scale.y),
        source.rotation,
        source.strength
    };
}

GpuMaterial gpuMaterial(const Material& source) {
    return {
        make_float4(source.baseColor.x, source.baseColor.y, source.baseColor.z,
                    source.baseColor.w),
        make_float3(source.emissive.x, source.emissive.y, source.emissive.z),
        make_float3(source.attenuationColor.x, source.attenuationColor.y,
                    source.attenuationColor.z),
        source.metallic,
        source.roughness,
        source.transmission,
        source.ior,
        source.thickness,
        source.attenuationDistance,
        source.emissiveStrength,
        source.dispersion,
        gpuTextureRef(source.baseColorTexture),
        gpuTextureRef(source.metallicRoughnessTexture),
        gpuTextureRef(source.normalTexture),
        gpuTextureRef(source.emissiveTexture),
        gpuTextureRef(source.transmissionTexture),
        gpuTextureRef(source.thicknessTexture)
    };
}

float3 point(const Mat4& matrix, const Vec3& value) {
    const auto& m = matrix.values;
    return make_float3(
        m[0] * value.x + m[4] * value.y + m[8] * value.z + m[12],
        m[1] * value.x + m[5] * value.y + m[9] * value.z + m[13],
        m[2] * value.x + m[6] * value.y + m[10] * value.z + m[14]);
}

float3 difference(float3 a, float3 b) {
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

float3 crossProduct(float3 a, float3 b) {
    return make_float3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                       a.x * b.y - a.y * b.x);
}

float length(float3 value) {
    return std::sqrt(value.x * value.x + value.y * value.y + value.z * value.z);
}

float luminance(float3 value) {
    return 0.2126f * value.x + 0.7152f * value.y + 0.0722f * value.z;
}

}

struct Renderer::Impl {
    struct GeometryState {
        Buffer vertices;
        Buffer indices;
        Buffer acceleration;
        OptixTraversableHandle handle = 0;
        std::uint32_t material = 0;
    };

    struct ImageState {
        Buffer pixels;
        std::uint32_t width = 0;
        std::uint32_t height = 0;
    };

    ~Impl() {
        if (denoiser)
            optixDenoiserDestroy(denoiser);
        if (pipeline)
            optixPipelineDestroy(pipeline);
        if (photonRaygen)
            optixProgramGroupDestroy(photonRaygen);
        if (raygen)
            optixProgramGroupDestroy(raygen);
        if (miss)
            optixProgramGroupDestroy(miss);
        if (hit)
            optixProgramGroupDestroy(hit);
        if (module)
            optixModuleDestroy(module);
        if (context)
            optixDeviceContextDestroy(context);
    }

    void initialize(const Scene& scene, std::uint32_t renderWidth,
                    std::uint32_t renderHeight, Profile profile) {
        if (scene.instances.empty() || scene.cameras.empty())
            throw std::runtime_error("The scene needs geometry and a perspective camera");
        width = renderWidth;
        height = renderHeight;
        launch.maxDepth = profile == Profile::Final ? 12U : 5U;

        cudaCheck(cudaFree(nullptr));
        optixCheck(optixInit(), "optixInit");
        OptixDeviceContextOptions contextOptions{};
        optixCheck(optixDeviceContextCreate(nullptr, &contextOptions, &context),
                   "optixDeviceContextCreate");

        uploadTextures(scene);
        uploadMaterials(scene);
        uploadEnvironment(scene);
        buildGeometry(scene);
        buildInstances(scene);
        buildLights(scene);
        configureCaustics(scene);
        createPipeline();
        createBindingTable();
        accumulation.resize(static_cast<std::size_t>(width) * height * sizeof(float4));
        albedoGuide.resize(accumulation.size());
        normalGuide.resize(accumulation.size());
        output.resize(static_cast<std::size_t>(width) * height * sizeof(uchar4));
        configureCamera(scene.cameras.front());
        reset();
        if (profile == Profile::Final)
            buildCaustics();
    }

    void uploadMaterials(const Scene& scene) {
        std::vector<GpuMaterial> gpuMaterials;
        gpuMaterials.reserve(scene.materials.size());
        for (const auto& material : scene.materials)
            gpuMaterials.push_back(gpuMaterial(material));
        materials.upload(gpuMaterials.data(), gpuMaterials.size() * sizeof(GpuMaterial));
    }

    void uploadTextures(const Scene& scene) {
        images.reserve(scene.images.size());
        for (const auto& source : scene.images) {
            if (source.encoded.size() > static_cast<std::size_t>(INT_MAX))
                throw std::runtime_error("Texture image is too large");
            int width;
            int height;
            int channels;
            auto* decoded = stbi_load_from_memory(
                reinterpret_cast<const stbi_uc*>(source.encoded.data()),
                static_cast<int>(source.encoded.size()), &width, &height, &channels, 4);
            if (!decoded)
                throw std::runtime_error("Unable to decode texture " + source.name +
                                         ": " + stbi_failure_reason());
            auto& image = images.emplace_back();
            image.width = static_cast<std::uint32_t>(width);
            image.height = static_cast<std::uint32_t>(height);
            image.pixels.upload(decoded,
                static_cast<std::size_t>(width) * height * sizeof(uchar4));
            stbi_image_free(decoded);
        }

        std::vector<GpuTexture> records;
        records.reserve(scene.textures.size());
        for (const auto& source : scene.textures) {
            const auto& image = images.at(static_cast<std::size_t>(source.image));
            records.push_back({
                reinterpret_cast<const uchar4*>(image.pixels.device()),
                image.width, image.height, source.wrapU, source.wrapV
            });
        }
        if (!records.empty()) {
            textures.upload(records.data(), records.size() * sizeof(GpuTexture));
            launch.textures = reinterpret_cast<const GpuTexture*>(textures.device());
        }
    }

    void uploadEnvironment(const Scene& scene) {
        const auto& source = scene.environment;
        launch.environmentRotation = source.rotation;
        launch.environmentStrength = source.strength;
        launch.exposure = source.exposure;
        if (source.pixels.empty())
            return;

        std::vector<float4> pixels;
        std::vector<float> cdf;
        pixels.reserve(source.pixels.size());
        cdf.reserve(source.pixels.size());
        double total = 0.0;
        for (std::uint32_t y = 0; y < source.height; ++y) {
            const float sine = std::sin(3.14159265f *
                (static_cast<float>(y) + 0.5f) / source.height);
            for (std::uint32_t x = 0; x < source.width; ++x) {
                const auto& pixel = source.pixels[y * source.width + x];
                pixels.push_back(make_float4(pixel.x, pixel.y, pixel.z, pixel.w));
                total += luminance(make_float3(pixel.x, pixel.y, pixel.z)) * sine;
                cdf.push_back(static_cast<float>(total));
            }
        }
        if (total <= 0.0)
            return;
        for (auto& value : cdf)
            value = static_cast<float>(value / total);

        environment.upload(pixels.data(), pixels.size() * sizeof(float4));
        environmentCdf.upload(cdf.data(), cdf.size() * sizeof(float));
        launch.environment =
            reinterpret_cast<const float4*>(environment.device());
        launch.environmentCdf =
            reinterpret_cast<const float*>(environmentCdf.device());
        launch.environmentWidth = source.width;
        launch.environmentHeight = source.height;
        environmentPower = static_cast<float>(
            total * 19.7392088 / static_cast<double>(source.width * source.height) *
            source.strength);
        launch.environmentWeight = environmentPower;
    }

    void buildGeometry(const Scene& scene) {
        geometries.reserve(scene.geometries.size());
        for (const auto& geometry : scene.geometries) {
            auto& state = geometries.emplace_back();
            std::vector<GpuVertex> vertices;
            vertices.reserve(geometry.vertices.size());
            for (const auto& vertex : geometry.vertices) {
                vertices.push_back({
                    make_float3(vertex.position.x, vertex.position.y, vertex.position.z),
                    make_float3(vertex.normal.x, vertex.normal.y, vertex.normal.z),
                    make_float4(vertex.tangent.x, vertex.tangent.y, vertex.tangent.z,
                                vertex.tangent.w),
                    make_float2(vertex.uv.x, vertex.uv.y),
                    make_float2(vertex.uv1.x, vertex.uv1.y)
                });
            }
            state.vertices.upload(vertices.data(), vertices.size() * sizeof(GpuVertex));
            state.indices.upload(geometry.indices.data(),
                                 geometry.indices.size() * sizeof(std::uint32_t));
            state.material = geometry.material;

            const CUdeviceptr deviceVertices = state.vertices.device();
            const std::uint32_t flags = OPTIX_GEOMETRY_FLAG_NONE;
            OptixBuildInput input{};
            input.type = OPTIX_BUILD_INPUT_TYPE_TRIANGLES;
            input.triangleArray.vertexBuffers = &deviceVertices;
            input.triangleArray.numVertices =
                static_cast<std::uint32_t>(geometry.vertices.size());
            input.triangleArray.vertexFormat = OPTIX_VERTEX_FORMAT_FLOAT3;
            input.triangleArray.vertexStrideInBytes = sizeof(GpuVertex);
            input.triangleArray.indexBuffer = state.indices.device();
            input.triangleArray.numIndexTriplets =
                static_cast<std::uint32_t>(geometry.indices.size() / 3);
            input.triangleArray.indexFormat = OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
            input.triangleArray.indexStrideInBytes = sizeof(uint3);
            input.triangleArray.flags = &flags;
            input.triangleArray.numSbtRecords = 1;

            OptixAccelBuildOptions options{};
            options.buildFlags = OPTIX_BUILD_FLAG_PREFER_FAST_TRACE;
            options.operation = OPTIX_BUILD_OPERATION_BUILD;
            OptixAccelBufferSizes sizes{};
            optixCheck(optixAccelComputeMemoryUsage(context, &options, &input, 1, &sizes),
                       "optixAccelComputeMemoryUsage");
            Buffer scratch(sizes.tempSizeInBytes);
            state.acceleration.resize(sizes.outputSizeInBytes);
            optixCheck(optixAccelBuild(context, nullptr, &options, &input, 1,
                                      scratch.device(), scratch.size(),
                                      state.acceleration.device(), state.acceleration.size(),
                                      &state.handle, nullptr, 0),
                       "optixAccelBuild");
        }
        cudaCheck(cudaDeviceSynchronize());
    }

    void buildInstances(const Scene& scene) {
        std::vector<OptixInstance> instances(scene.instances.size());
        for (std::size_t index = 0; index < scene.instances.size(); ++index) {
            const auto& source = scene.instances[index];
            auto& target = instances[index];
            for (std::size_t row = 0; row < 3; ++row)
                for (std::size_t column = 0; column < 4; ++column)
                    target.transform[row * 4 + column] =
                        source.transform.values[column * 4 + row];
            target.instanceId = static_cast<std::uint32_t>(index);
            target.sbtOffset = source.geometry;
            target.visibilityMask = 255;
            target.flags = OPTIX_INSTANCE_FLAG_NONE;
            target.traversableHandle = geometries.at(source.geometry).handle;
        }

        Buffer instanceData;
        instanceData.upload(instances.data(), instances.size() * sizeof(OptixInstance));
        OptixBuildInput input{};
        input.type = OPTIX_BUILD_INPUT_TYPE_INSTANCES;
        input.instanceArray.instances = instanceData.device();
        input.instanceArray.numInstances = static_cast<std::uint32_t>(instances.size());
        OptixAccelBuildOptions options{};
        options.buildFlags = OPTIX_BUILD_FLAG_PREFER_FAST_TRACE;
        options.operation = OPTIX_BUILD_OPERATION_BUILD;
        OptixAccelBufferSizes sizes{};
        optixCheck(optixAccelComputeMemoryUsage(context, &options, &input, 1, &sizes),
                   "optixAccelComputeMemoryUsage");
        Buffer scratch(sizes.tempSizeInBytes);
        acceleration.resize(sizes.outputSizeInBytes);
        optixCheck(optixAccelBuild(context, nullptr, &options, &input, 1,
                                  scratch.device(), scratch.size(),
                                  acceleration.device(), acceleration.size(),
                                  &sceneHandle, nullptr, 0),
                   "optixAccelBuild");
        cudaCheck(cudaDeviceSynchronize());
    }

    void buildLights(const Scene& scene) {
        std::vector<GpuLight> source;
        for (std::size_t instanceIndex = 0; instanceIndex < scene.instances.size();
             ++instanceIndex) {
            const auto& instance = scene.instances[instanceIndex];
            const auto& geometry = scene.geometries[instance.geometry];
            const auto& material = scene.materials[geometry.material];
            const float3 emission = make_float3(
                material.emissive.x * material.emissiveStrength,
                material.emissive.y * material.emissiveStrength,
                material.emissive.z * material.emissiveStrength);
            if (luminance(emission) <= 0.0f)
                continue;

            for (std::size_t index = 0; index < geometry.indices.size(); index += 3) {
                const float3 a = point(instance.transform,
                    geometry.vertices[geometry.indices[index]].position);
                const float3 b = point(instance.transform,
                    geometry.vertices[geometry.indices[index + 1]].position);
                const float3 c = point(instance.transform,
                    geometry.vertices[geometry.indices[index + 2]].position);
                const float3 areaVector =
                    crossProduct(difference(b, a), difference(c, a));
                const float twiceArea = length(areaVector);
                if (twiceArea <= 1e-8f)
                    continue;
                const float area = 0.5f * twiceArea;
                source.push_back({
                    a, b, c,
                    make_float3(areaVector.x / twiceArea, areaVector.y / twiceArea,
                                areaVector.z / twiceArea),
                    emission, area, area * luminance(emission),
                    0.0f, 0.0f, 0.0f,
                    static_cast<std::uint32_t>(instanceIndex),
                    static_cast<std::uint32_t>(index / 3), 3U
                });
            }
        }

        for (const auto& light : scene.lights) {
            const auto& m = light.transform.values;
            const float3 emission = make_float3(
                light.color.x * light.intensity,
                light.color.y * light.intensity,
                light.color.z * light.intensity);
            const float3 direction = normalized(-m[8], -m[9], -m[10]);
            const float solidAngle = light.type == 1U
                ? 6.2831853f * (1.0f - std::cos(light.outerCone))
                : light.type == 2U ? 12.566371f : 1.0f;
            source.push_back({
                make_float3(m[12], m[13], m[14]), direction, {}, {}, emission,
                0.0f, luminance(emission) * solidAngle, light.range,
                light.innerCone, light.outerCone,
                0xffffffffU, 0xffffffffU, light.type
            });
        }
        if (environmentPower > 0.0f) {
            source.push_back({
                {}, {}, {}, {}, {}, 0.0f, environmentPower,
                0.0f, 0.0f, 0.0f, 0xffffffffU, 0xffffffffU, 4U
            });
        }

        launch.lightWeight = 0.0f;
        for (const auto& light : source)
            launch.lightWeight += light.weight;
        launch.lightCount = static_cast<std::uint32_t>(source.size());
        if (!source.empty()) {
            lights.upload(source.data(), source.size() * sizeof(GpuLight));
            launch.lights = reinterpret_cast<const GpuLight*>(lights.device());
        }
    }

    void configureCaustics(const Scene& scene) {
        const float limit = std::numeric_limits<float>::max();
        float3 minimum = make_float3(limit, limit, limit);
        float3 maximum = make_float3(-limit, -limit, -limit);
        for (const auto& instance : scene.instances) {
            for (const auto& vertex : scene.geometries[instance.geometry].vertices) {
                const float3 value = point(instance.transform, vertex.position);
                minimum = make_float3(std::min(minimum.x, value.x),
                                      std::min(minimum.y, value.y),
                                      std::min(minimum.z, value.z));
                maximum = make_float3(std::max(maximum.x, value.x),
                                      std::max(maximum.y, value.y),
                                      std::max(maximum.z, value.z));
            }
        }
        launch.sceneCenter = make_float3(
            (minimum.x + maximum.x) * 0.5f,
            (minimum.y + maximum.y) * 0.5f,
            (minimum.z + maximum.z) * 0.5f);
        launch.sceneRadius = length(difference(maximum, minimum)) * 0.525f;
        causticRadius = std::max(launch.sceneRadius * 0.018f, 0.01f);
        launch.photonBucketCount = 1U << 17U;
        launch.photonBucketSize = 8U;
        launch.photonEmissions = 1U << 20U;
        photons.resize(static_cast<std::size_t>(launch.photonBucketCount) *
                       launch.photonBucketSize * sizeof(Photon));
        photonBuckets.resize(static_cast<std::size_t>(launch.photonBucketCount) *
                             sizeof(std::uint32_t));
        launch.photons = static_cast<Photon*>(photons.data());
        launch.photonBuckets =
            static_cast<std::uint32_t*>(photonBuckets.data());
    }

    void createPipeline() {
        OptixModuleCompileOptions moduleOptions{};
        moduleOptions.optLevel = OPTIX_COMPILE_OPTIMIZATION_LEVEL_3;
        moduleOptions.debugLevel = OPTIX_COMPILE_DEBUG_LEVEL_NONE;
        pipelineOptions.usesMotionBlur = false;
        pipelineOptions.traversableGraphFlags =
            OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_SINGLE_LEVEL_INSTANCING;
        pipelineOptions.numPayloadValues = 2;
        pipelineOptions.numAttributeValues = 2;
        pipelineOptions.exceptionFlags = OPTIX_EXCEPTION_FLAG_NONE;
        pipelineOptions.pipelineLaunchParamsVariableName = "params";
        pipelineOptions.usesPrimitiveTypeFlags = OPTIX_PRIMITIVE_TYPE_FLAGS_TRIANGLE;

        std::array<char, 4096> log{};
        std::size_t logSize = log.size();
        auto result = optixModuleCreate(
            context, &moduleOptions, &pipelineOptions,
            reinterpret_cast<const char*>(deviceCode), deviceCodeSize,
            log.data(), &logSize, &module);
        if (result != OPTIX_SUCCESS)
            throw std::runtime_error("optixModuleCreate: " +
                                     std::string(log.data(), logSize));

        raygen = createProgram(OPTIX_PROGRAM_GROUP_KIND_RAYGEN, "__raygen__render");
        photonRaygen =
            createProgram(OPTIX_PROGRAM_GROUP_KIND_RAYGEN, "__raygen__photons");
        miss = createProgram(OPTIX_PROGRAM_GROUP_KIND_MISS, "__miss__environment");
        hit = createProgram(OPTIX_PROGRAM_GROUP_KIND_HITGROUP, "__closesthit__surface");

        const std::array groups{raygen, photonRaygen, miss, hit};
        OptixPipelineLinkOptions linkOptions{};
        linkOptions.maxTraceDepth = 1;
        logSize = log.size();
        result = optixPipelineCreate(context, &pipelineOptions, &linkOptions,
                                     groups.data(), groups.size(), log.data(),
                                     &logSize, &pipeline);
        if (result != OPTIX_SUCCESS)
            throw std::runtime_error("optixPipelineCreate: " +
                                     std::string(log.data(), logSize));

        OptixStackSizes stack{};
        for (const auto group : groups)
            optixCheck(optixUtilAccumulateStackSizes(group, &stack, pipeline),
                       "optixUtilAccumulateStackSizes");
        std::uint32_t traversal;
        std::uint32_t state;
        std::uint32_t continuation;
        optixCheck(optixUtilComputeStackSizes(&stack, 1, 0, 0, &traversal, &state,
                                              &continuation),
                   "optixUtilComputeStackSizes");
        optixCheck(optixPipelineSetStackSize(pipeline, traversal, state, continuation, 2),
                   "optixPipelineSetStackSize");
    }

    OptixProgramGroup createProgram(OptixProgramGroupKind kind, const char* entry) {
        OptixProgramGroupDesc description{};
        description.kind = kind;
        if (kind == OPTIX_PROGRAM_GROUP_KIND_RAYGEN) {
            description.raygen.module = module;
            description.raygen.entryFunctionName = entry;
        } else if (kind == OPTIX_PROGRAM_GROUP_KIND_MISS) {
            description.miss.module = module;
            description.miss.entryFunctionName = entry;
        } else {
            description.hitgroup.moduleCH = module;
            description.hitgroup.entryFunctionNameCH = entry;
        }

        OptixProgramGroup group = nullptr;
        OptixProgramGroupOptions options{};
        std::array<char, 2048> log{};
        std::size_t logSize = log.size();
        const auto result = optixProgramGroupCreate(context, &description, 1, &options,
                                                    log.data(), &logSize, &group);
        if (result != OPTIX_SUCCESS)
            throw std::runtime_error("optixProgramGroupCreate: " +
                                     std::string(log.data(), logSize));
        return group;
    }

    void createBindingTable() {
        Record<Empty> raygenRecord{};
        optixCheck(optixSbtRecordPackHeader(raygen, &raygenRecord),
                   "optixSbtRecordPackHeader");
        raygenRecordBuffer.upload(&raygenRecord, sizeof(raygenRecord));

        Record<Empty> photonRecord{};
        optixCheck(optixSbtRecordPackHeader(photonRaygen, &photonRecord),
                   "optixSbtRecordPackHeader");
        photonRaygenRecordBuffer.upload(&photonRecord, sizeof(photonRecord));

        Record<Empty> missRecord{};
        optixCheck(optixSbtRecordPackHeader(miss, &missRecord),
                   "optixSbtRecordPackHeader");
        missRecordBuffer.upload(&missRecord, sizeof(missRecord));

        std::vector<Record<HitGroupData>> hitRecords(geometries.size());
        for (std::size_t index = 0; index < geometries.size(); ++index) {
            auto& record = hitRecords[index];
            optixCheck(optixSbtRecordPackHeader(hit, &record),
                       "optixSbtRecordPackHeader");
            record.data = {
                reinterpret_cast<const GpuVertex*>(geometries[index].vertices.device()),
                reinterpret_cast<const uint3*>(geometries[index].indices.device()),
                geometries[index].material
            };
        }
        hitRecordBuffer.upload(hitRecords.data(), hitRecords.size() * sizeof(hitRecords[0]));

        bindingTable.raygenRecord = raygenRecordBuffer.device();
        bindingTable.missRecordBase = missRecordBuffer.device();
        bindingTable.missRecordStrideInBytes = sizeof(missRecord);
        bindingTable.missRecordCount = 1;
        bindingTable.hitgroupRecordBase = hitRecordBuffer.device();
        bindingTable.hitgroupRecordStrideInBytes = sizeof(hitRecords[0]);
        bindingTable.hitgroupRecordCount =
            static_cast<std::uint32_t>(hitRecords.size());
    }

    void configureCamera(const Camera& camera) {
        const auto& matrix = camera.transform.values;
        const float3 right = normalized(matrix[0], matrix[1], matrix[2]);
        const float3 up = normalized(matrix[4], matrix[5], matrix[6]);
        const float3 direction = normalized(matrix[8], matrix[9], matrix[10]);
        const float3 forward = make_float3(-direction.x, -direction.y, -direction.z);
        const float aspect = camera.aspectRatio > 0.0f
            ? camera.aspectRatio
            : static_cast<float>(width) / static_cast<float>(height);
        const float scale = std::tan(camera.verticalFov * 0.5f);
        launch.eye = make_float3(matrix[12], matrix[13], matrix[14]);
        launch.cameraU = make_float3(right.x * scale * aspect,
                                     right.y * scale * aspect,
                                     right.z * scale * aspect);
        launch.cameraV = make_float3(up.x * scale, up.y * scale, up.z * scale);
        launch.cameraW = forward;
        launch.lensU = right;
        launch.lensV = up;
        launch.aperture = camera.aperture;
        launch.focusDistance = camera.focusDistance;
        launch.width = width;
        launch.height = height;
        launch.accumulation = static_cast<float4*>(accumulation.data());
        launch.albedoGuide = static_cast<float4*>(albedoGuide.data());
        launch.normalGuide = static_cast<float4*>(normalGuide.data());
        launch.materials = reinterpret_cast<const GpuMaterial*>(materials.device());
        launch.textures = reinterpret_cast<const GpuTexture*>(textures.device());
        launch.scene = sceneHandle;
        parameters.resize(sizeof(LaunchParams));
    }

    void reset() {
        cudaCheck(cudaMemset(accumulation.data(), 0, accumulation.size()));
        cudaCheck(cudaMemset(albedoGuide.data(), 0, albedoGuide.size()));
        cudaCheck(cudaMemset(normalGuide.data(), 0, normalGuide.size()));
        sample = 0;
        denoisedReady = false;
    }

    void setCamera(const Camera& camera) {
        configureCamera(camera);
        reset();
    }

    void setProfile(Profile profile) {
        launch.maxDepth = profile == Profile::Final ? 12U : 5U;
        if (profile == Profile::Final)
            buildCaustics();
        else
            launch.photonRadius = 0.0f;
        reset();
    }

    void buildCaustics() {
        launch.photonRadius = causticRadius;
        if (causticsReady)
            return;
        cudaCheck(cudaMemset(photonBuckets.data(), 0, photonBuckets.size()));
        cudaCheck(cudaMemcpy(parameters.data(), &launch, sizeof(launch),
                             cudaMemcpyHostToDevice));
        const auto record = bindingTable.raygenRecord;
        bindingTable.raygenRecord = photonRaygenRecordBuffer.device();
        const auto result = optixLaunch(
            pipeline, nullptr, parameters.device(), sizeof(launch), &bindingTable,
            launch.photonEmissions, 1, 1);
        bindingTable.raygenRecord = record;
        optixCheck(result, "optixLaunch");
        cudaCheck(cudaDeviceSynchronize());
        causticsReady = true;
    }

    void render(void* output) {
        launch.display = nullptr;
        launch.output = static_cast<uchar4*>(output ? output : this->output.data());
        launch.sample = sample;
        cudaCheck(cudaMemcpy(parameters.data(), &launch, sizeof(launch),
                             cudaMemcpyHostToDevice));
        optixCheck(optixLaunch(pipeline, nullptr, parameters.device(), sizeof(launch),
                              &bindingTable, width, height, 1),
                   "optixLaunch");
        cudaCheck(cudaGetLastError());
        ++sample;
    }

    OptixImage2D image(CUdeviceptr data) const {
        return {
            data, width, height, width * static_cast<unsigned int>(sizeof(float4)),
            static_cast<unsigned int>(sizeof(float4)), OPTIX_PIXEL_FORMAT_FLOAT4
        };
    }

    void denoiseImage() {
        if (!denoiser) {
            OptixDenoiserOptions options{};
            options.guideAlbedo = 1;
            options.guideNormal = 1;
            options.denoiseAlpha = OPTIX_DENOISER_ALPHA_MODE_COPY;
            optixCheck(optixDenoiserCreate(context, OPTIX_DENOISER_MODEL_KIND_AOV,
                                           &options, &denoiser),
                       "optixDenoiserCreate");
            OptixDenoiserSizes sizes{};
            optixCheck(optixDenoiserComputeMemoryResources(
                           denoiser, width, height, &sizes),
                       "optixDenoiserComputeMemoryResources");
            denoiserState.resize(sizes.stateSizeInBytes);
            denoiserScratch.resize(std::max(
                sizes.withoutOverlapScratchSizeInBytes,
                sizes.computeAverageColorSizeInBytes));
            denoised.resize(accumulation.size());
            denoiserAverageColor.resize(3 * sizeof(float));
            optixCheck(optixDenoiserSetup(
                           denoiser, nullptr, width, height, denoiserState.device(),
                           denoiserState.size(), denoiserScratch.device(),
                           denoiserScratch.size()),
                       "optixDenoiserSetup");
        }

        const auto input = image(accumulation.device());
        optixCheck(optixDenoiserComputeAverageColor(
                       denoiser, nullptr, &input, denoiserAverageColor.device(),
                       denoiserScratch.device(), denoiserScratch.size()),
                   "optixDenoiserComputeAverageColor");
        OptixDenoiserParams denoiserParams{};
        denoiserParams.hdrAverageColor = denoiserAverageColor.device();
        OptixDenoiserGuideLayer guides{};
        guides.albedo = image(albedoGuide.device());
        guides.normal = image(normalGuide.device());
        OptixDenoiserLayer layer{};
        layer.input = input;
        layer.output = image(denoised.device());
        optixCheck(optixDenoiserInvoke(
                       denoiser, nullptr, &denoiserParams, denoiserState.device(),
                       denoiserState.size(), &guides, &layer, 1, 0, 0,
                       denoiserScratch.device(), denoiserScratch.size()),
                   "optixDenoiserInvoke");

        launch.display = static_cast<const float4*>(denoised.data());
        launch.output = static_cast<uchar4*>(output.data());
        cudaCheck(cudaMemcpy(parameters.data(), &launch, sizeof(launch),
                             cudaMemcpyHostToDevice));
        optixCheck(optixLaunch(pipeline, nullptr, parameters.device(), sizeof(launch),
                              &bindingTable, width, height, 1),
                   "optixLaunch");
        cudaCheck(cudaGetLastError());
        denoisedReady = true;
    }

    std::vector<std::uint8_t> pixels() const {
        std::vector<std::uint8_t> result(output.size());
        cudaCheck(cudaMemcpy(result.data(), output.data(), output.size(),
                             cudaMemcpyDeviceToHost));
        return result;
    }

    void copyOutput(void* destination) const {
        cudaCheck(cudaMemcpy(destination, output.data(), output.size(),
                             cudaMemcpyDeviceToDevice));
    }

    std::vector<float> linearPixels() const {
        const Buffer& source = denoisedReady ? denoised : accumulation;
        std::vector<float> result(source.size() / sizeof(float));
        cudaCheck(cudaMemcpy(result.data(), source.data(), source.size(),
                             cudaMemcpyDeviceToHost));
        return result;
    }

    std::uint32_t width = 0;
    std::uint32_t height = 0;
    std::uint32_t sample = 0;
    OptixDeviceContext context = nullptr;
    OptixModule module = nullptr;
    OptixProgramGroup raygen = nullptr;
    OptixProgramGroup photonRaygen = nullptr;
    OptixProgramGroup miss = nullptr;
    OptixProgramGroup hit = nullptr;
    OptixDenoiser denoiser = nullptr;
    OptixPipeline pipeline = nullptr;
    OptixPipelineCompileOptions pipelineOptions{};
    OptixShaderBindingTable bindingTable{};
    OptixTraversableHandle sceneHandle = 0;
    std::vector<GeometryState> geometries;
    std::vector<ImageState> images;
    Buffer materials;
    Buffer textures;
    Buffer lights;
    Buffer environment;
    Buffer environmentCdf;
    Buffer acceleration;
    Buffer accumulation;
    Buffer albedoGuide;
    Buffer normalGuide;
    Buffer output;
    Buffer denoised;
    Buffer denoiserState;
    Buffer denoiserScratch;
    Buffer denoiserAverageColor;
    Buffer parameters;
    Buffer raygenRecordBuffer;
    Buffer missRecordBuffer;
    Buffer hitRecordBuffer;
    Buffer photonRaygenRecordBuffer;
    Buffer photons;
    Buffer photonBuckets;
    LaunchParams launch{};
    float environmentPower = 0.0f;
    float causticRadius = 0.0f;
    bool denoisedReady = false;
    bool causticsReady = false;
};

Renderer::Renderer(const Scene& scene, std::uint32_t width, std::uint32_t height,
                   Profile profile)
    : impl_(std::make_unique<Impl>()) {
    impl_->initialize(scene, width, height, profile);
}

Renderer::~Renderer() = default;

void Renderer::render(void* output) {
    impl_->render(output);
}

void Renderer::denoise() {
    impl_->denoiseImage();
}

void Renderer::setCamera(const Camera& camera) {
    impl_->setCamera(camera);
}

void Renderer::setProfile(Profile profile) {
    impl_->setProfile(profile);
}

void Renderer::copyOutput(void* output) const {
    impl_->copyOutput(output);
}

std::vector<std::uint8_t> Renderer::pixels() const {
    return impl_->pixels();
}

std::vector<float> Renderer::linearPixels() const {
    return impl_->linearPixels();
}

std::uint32_t Renderer::samples() const {
    return impl_->sample;
}

}
