#include "gpu_shared.hpp"

#include <optix_device.h>

namespace rt {

extern "C" {
__constant__ LaunchParams params;
}

static __forceinline__ __device__ float3 operator+(float3 a, float3 b) {
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

static __forceinline__ __device__ float3 operator+(float3 value, float scalar) {
    return make_float3(value.x + scalar, value.y + scalar, value.z + scalar);
}

static __forceinline__ __device__ float3 operator-(float3 a, float3 b) {
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

static __forceinline__ __device__ float3 operator-(float3 value) {
    return make_float3(-value.x, -value.y, -value.z);
}

static __forceinline__ __device__ float3 operator*(float3 a, float3 b) {
    return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}

static __forceinline__ __device__ float3 operator*(float3 value, float scalar) {
    return make_float3(value.x * scalar, value.y * scalar, value.z * scalar);
}

static __forceinline__ __device__ float3 operator*(float scalar, float3 value) {
    return value * scalar;
}

static __forceinline__ __device__ float3 operator/(float3 value, float scalar) {
    return value * (1.0f / scalar);
}

static __forceinline__ __device__ float3 operator/(float3 a, float3 b) {
    return make_float3(a.x / b.x, a.y / b.y, a.z / b.z);
}

static __forceinline__ __device__ void operator+=(float3& a, float3 b) {
    a = a + b;
}

static __forceinline__ __device__ void operator*=(float3& a, float3 b) {
    a = a * b;
}

static __forceinline__ __device__ void operator/=(float3& value, float scalar) {
    value = value / scalar;
}

static __forceinline__ __device__ float dot(float3 a, float3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static __forceinline__ __device__ float3 cross(float3 a, float3 b) {
    return make_float3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                       a.x * b.y - a.y * b.x);
}

static __forceinline__ __device__ float3 normalize(float3 value) {
    return value / sqrtf(dot(value, value));
}

static __forceinline__ __device__ float3 rgb(float4 value) {
    return make_float3(value.x, value.y, value.z);
}

static __forceinline__ __device__ void packPointer(void* pointer, unsigned int& high,
                                                   unsigned int& low) {
    const auto value = reinterpret_cast<unsigned long long>(pointer);
    high = static_cast<unsigned int>(value >> 32);
    low = static_cast<unsigned int>(value);
}

static __forceinline__ __device__ void* unpackPointer() {
    const auto value = static_cast<unsigned long long>(optixGetPayload_0()) << 32 |
                       optixGetPayload_1();
    return reinterpret_cast<void*>(value);
}

static __forceinline__ __device__ float random(unsigned int& state) {
    state = state * 747796405U + 2891336453U;
    const auto word = ((state >> ((state >> 28U) + 4U)) ^ state) * 277803737U;
    return static_cast<float>((word >> 22U) ^ word) * 0x1p-32f;
}

static __forceinline__ __device__ float3 sky(float3 direction) {
    const float t = 0.5f * (direction.y + 1.0f);
    return make_float3(0.03f, 0.04f, 0.06f) * (1.0f - t) +
           make_float3(0.35f, 0.48f, 0.70f) * t;
}

static __forceinline__ __device__ float3 cosineDirection(float3 normal,
                                                          unsigned int& rng) {
    const float phi = 6.2831853f * random(rng);
    const float radius = sqrtf(random(rng));
    const float z = sqrtf(fmaxf(0.0f, 1.0f - radius * radius));
    const float3 tangent = normalize(fabsf(normal.x) > 0.5f
        ? cross(make_float3(0.0f, 1.0f, 0.0f), normal)
        : cross(make_float3(1.0f, 0.0f, 0.0f), normal));
    const float3 bitangent = cross(normal, tangent);
    return normalize(tangent * (radius * cosf(phi)) +
                     bitangent * (radius * sinf(phi)) + normal * z);
}

static __forceinline__ __device__ float3 aces(float3 color) {
    color = color * (2.51f * color + 0.03f) /
            (color * (2.43f * color + 0.59f) + 0.14f);
    return make_float3(
        powf(fminf(fmaxf(color.x, 0.0f), 1.0f), 1.0f / 2.2f),
        powf(fminf(fmaxf(color.y, 0.0f), 1.0f), 1.0f / 2.2f),
        powf(fminf(fmaxf(color.z, 0.0f), 1.0f), 1.0f / 2.2f));
}

static __forceinline__ __device__ Hit trace(float3 origin, float3 direction) {
    Hit hit{};
    unsigned int high;
    unsigned int low;
    packPointer(&hit, high, low);
    optixTrace(params.scene, origin, direction, 0.001f, 1e16f, 0.0f, 255,
               OPTIX_RAY_FLAG_DISABLE_ANYHIT, 0, 1, 0, high, low);
    return hit;
}

extern "C" __global__ void __raygen__render() {
    const uint3 pixel = optixGetLaunchIndex();
    const unsigned int index = pixel.y * params.width + pixel.x;
    unsigned int rng = index * 9781U + params.sample * 6271U + 0x68bc21ebU;
    const float2 jitter = make_float2(random(rng), random(rng));
    const float2 screen = make_float2(
        (static_cast<float>(pixel.x) + jitter.x) / static_cast<float>(params.width),
        (static_cast<float>(pixel.y) + jitter.y) / static_cast<float>(params.height));

    float3 origin = params.eye;
    float3 direction = normalize(params.cameraW +
        (2.0f * screen.x - 1.0f) * params.cameraU +
        (2.0f * screen.y - 1.0f) * params.cameraV);
    float3 throughput = make_float3(1.0f, 1.0f, 1.0f);
    float3 radiance = make_float3(0.0f, 0.0f, 0.0f);

    for (unsigned int depth = 0; depth < 8; ++depth) {
        const Hit hit = trace(origin, direction);
        if (!hit.found) {
            radiance += throughput * sky(direction);
            break;
        }

        const GpuMaterial material = params.materials[hit.material];
        radiance += throughput * material.emissive * material.emissiveStrength;
        throughput *= rgb(material.baseColor);
        origin = hit.position + hit.normal * 0.001f;
        direction = cosineDirection(hit.normal, rng);

        if (depth > 2) {
            const float survival = fminf(fmaxf(fmaxf(throughput.x, throughput.y),
                                               throughput.z), 0.05f);
            if (random(rng) > survival)
                break;
            throughput /= survival;
        }
    }

    const float4 previous = params.accumulation[index];
    const float weight = 1.0f / static_cast<float>(params.sample + 1U);
    const float3 accumulated = rgb(previous) + (radiance - rgb(previous)) * weight;
    params.accumulation[index] =
        make_float4(accumulated.x, accumulated.y, accumulated.z, 1.0f);
    const float3 mapped = aces(accumulated);
    params.output[index] = make_uchar4(
        static_cast<unsigned char>(mapped.x * 255.0f + 0.5f),
        static_cast<unsigned char>(mapped.y * 255.0f + 0.5f),
        static_cast<unsigned char>(mapped.z * 255.0f + 0.5f), 255);
}

extern "C" __global__ void __miss__environment() {
    static_cast<Hit*>(unpackPointer())->found = false;
}

extern "C" __global__ void __closesthit__surface() {
    auto* hit = static_cast<Hit*>(unpackPointer());
    const auto* data = reinterpret_cast<const HitGroupData*>(optixGetSbtDataPointer());
    const uint3 triangle = data->indices[optixGetPrimitiveIndex()];
    const float2 barycentrics = optixGetTriangleBarycentrics();
    const float b0 = 1.0f - barycentrics.x - barycentrics.y;
    float3 normal = data->vertices[triangle.x].normal * b0 +
                    data->vertices[triangle.y].normal * barycentrics.x +
                    data->vertices[triangle.z].normal * barycentrics.y;
    if (dot(normal, normal) < 1e-12f) {
        const float3 a = data->vertices[triangle.x].position;
        const float3 b = data->vertices[triangle.y].position;
        const float3 c = data->vertices[triangle.z].position;
        normal = cross(b - a, c - a);
    }

    normal = normalize(optixTransformNormalFromObjectToWorldSpace(normal));
    const float3 ray = optixGetWorldRayDirection();
    hit->normal = dot(normal, ray) < 0.0f ? normal : -normal;
    hit->position = optixGetWorldRayOrigin() + optixGetRayTmax() * ray;
    hit->material = data->material;
    hit->found = true;
}

}
