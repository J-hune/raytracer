#include "gpu_shared.hpp"
#include "photon_hash.cuh"

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

static __forceinline__ __device__ void operator*=(float3& value, float scalar) {
    value = value * scalar;
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

static __forceinline__ __device__ float linear(float value) {
    return value <= 0.04045f ? value / 12.92f
                            : powf((value + 0.055f) / 1.055f, 2.4f);
}

static __forceinline__ __device__ int wrapped(int value, unsigned int size,
                                               unsigned int mode) {
    if (mode == 33071U)
        return min(max(value, 0), static_cast<int>(size) - 1);
    if (mode == 33648U) {
        const int period = static_cast<int>(size) * 2;
        int coordinate = (value % period + period) % period;
        return coordinate < static_cast<int>(size)
            ? coordinate : period - coordinate - 1;
    }
    return (value % static_cast<int>(size) + static_cast<int>(size)) %
           static_cast<int>(size);
}

static __forceinline__ __device__ float4 textureTexel(
    const GpuTexture& texture, int x, int y) {
    const uchar4 value = texture.pixels[
        wrapped(y, texture.height, texture.wrapV) * texture.width +
        wrapped(x, texture.width, texture.wrapU)];
    return make_float4(value.x / 255.0f, value.y / 255.0f,
                       value.z / 255.0f, value.w / 255.0f);
}

static __forceinline__ __device__ float4 texture(
    const GpuTextureRef& reference, const Hit& hit, bool srgb) {
    if (reference.texture < 0)
        return make_float4(1.0f, 1.0f, 1.0f, 1.0f);
    const GpuTexture source = params.textures[reference.texture];
    const float2 input = reference.texCoord == 1U ? hit.uv1 : hit.uv;
    const float2 scaled =
        make_float2(input.x * reference.scale.x, input.y * reference.scale.y);
    const float cosine = cosf(reference.rotation);
    const float sine = sinf(reference.rotation);
    const float2 uv = make_float2(
        reference.offset.x + cosine * scaled.x - sine * scaled.y,
        reference.offset.y + sine * scaled.x + cosine * scaled.y);
    const float x = uv.x * source.width - 0.5f;
    const float y = uv.y * source.height - 0.5f;
    const int x0 = static_cast<int>(floorf(x));
    const int y0 = static_cast<int>(floorf(y));
    const float tx = x - floorf(x);
    const float ty = y - floorf(y);
    const float4 a = textureTexel(source, x0, y0);
    const float4 b = textureTexel(source, x0 + 1, y0);
    const float4 c = textureTexel(source, x0, y0 + 1);
    const float4 d = textureTexel(source, x0 + 1, y0 + 1);
    float4 value = make_float4(
        (a.x * (1.0f - tx) + b.x * tx) * (1.0f - ty) +
            (c.x * (1.0f - tx) + d.x * tx) * ty,
        (a.y * (1.0f - tx) + b.y * tx) * (1.0f - ty) +
            (c.y * (1.0f - tx) + d.y * tx) * ty,
        (a.z * (1.0f - tx) + b.z * tx) * (1.0f - ty) +
            (c.z * (1.0f - tx) + d.z * tx) * ty,
        (a.w * (1.0f - tx) + b.w * tx) * (1.0f - ty) +
            (c.w * (1.0f - tx) + d.w * tx) * ty);
    if (srgb) {
        value.x = linear(value.x);
        value.y = linear(value.y);
        value.z = linear(value.z);
    }
    return value;
}

static __forceinline__ __device__ GpuMaterial textured(
    GpuMaterial material, const Hit& hit) {
    const float4 base = texture(material.baseColorTexture, hit, true);
    material.baseColor.x *= base.x;
    material.baseColor.y *= base.y;
    material.baseColor.z *= base.z;
    material.baseColor.w *= base.w;
    const float4 pbr = texture(material.metallicRoughnessTexture, hit, false);
    material.roughness *= pbr.y;
    material.metallic *= pbr.z;
    const float4 emissive = texture(material.emissiveTexture, hit, true);
    material.emissive.x *= emissive.x;
    material.emissive.y *= emissive.y;
    material.emissive.z *= emissive.z;
    material.transmission *= texture(material.transmissionTexture, hit, false).x;
    material.thickness *= texture(material.thicknessTexture, hit, false).y;
    return material;
}

static __forceinline__ __device__ float3 mappedNormal(
    const Hit& hit, const GpuMaterial& material) {
    if (material.normalTexture.texture < 0)
        return hit.normal;
    const float4 sample = texture(material.normalTexture, hit, false);
    const float3 tangent =
        normalize(make_float3(hit.tangent.x, hit.tangent.y, hit.tangent.z));
    const float3 local = normalize(make_float3(
        (sample.x * 2.0f - 1.0f) * material.normalTexture.strength,
        (sample.y * 2.0f - 1.0f) * material.normalTexture.strength,
        sample.z * 2.0f - 1.0f));
    return normalize(tangent * local.x +
        cross(hit.normal, tangent) * (local.y * hit.tangent.w) +
        hit.normal * local.z);
}

static __forceinline__ __device__ float maximum(float3 value) {
    return fmaxf(fmaxf(value.x, value.y), value.z);
}

static __forceinline__ __device__ float saturate(float value) {
    return fminf(fmaxf(value, 0.0f), 1.0f);
}

static __forceinline__ __device__ float3 reflect(float3 direction, float3 normal) {
    return direction - 2.0f * dot(direction, normal) * normal;
}

static __forceinline__ __device__ bool refract(float3 direction, float3 normal,
                                                float eta, float3& result) {
    const float cosine = fminf(dot(-direction, normal), 1.0f);
    const float3 perpendicular = eta * (direction + cosine * normal);
    const float parallelSquared = 1.0f - dot(perpendicular, perpendicular);
    if (parallelSquared < 0.0f)
        return false;
    result = perpendicular - sqrtf(parallelSquared) * normal;
    return true;
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

static __forceinline__ __device__ unsigned int seeded(unsigned int a, unsigned int b) {
    unsigned int value = a * 0x9e3779b9U ^ b * 0x85ebca6bU;
    value ^= value >> 16U;
    value *= 0x7feb352dU;
    value ^= value >> 15U;
    value *= 0x846ca68bU;
    return value ^ (value >> 16U);
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

static __forceinline__ __device__ float2 environmentUv(float3 direction) {
    float u = atan2f(direction.z, direction.x) * 0.15915494f + 0.5f +
              params.environmentRotation * 0.15915494f;
    u -= floorf(u);
    return make_float2(u, acosf(fminf(fmaxf(direction.y, -1.0f), 1.0f)) *
                              0.31830989f);
}

static __forceinline__ __device__ float3 environmentTexel(int x, int y) {
    x = (x % static_cast<int>(params.environmentWidth) +
         static_cast<int>(params.environmentWidth)) %
        static_cast<int>(params.environmentWidth);
    y = min(max(y, 0), static_cast<int>(params.environmentHeight) - 1);
    return rgb(params.environment[y * params.environmentWidth + x]);
}

static __forceinline__ __device__ float3 environment(float3 direction) {
    if (!params.environment)
        return sky(direction);

    const float2 uv = environmentUv(direction);
    const float x = uv.x * params.environmentWidth - 0.5f;
    const float y = uv.y * params.environmentHeight - 0.5f;
    const int x0 = static_cast<int>(floorf(x));
    const int y0 = static_cast<int>(floorf(y));
    const float tx = x - floorf(x);
    const float ty = y - floorf(y);
    const float3 a = environmentTexel(x0, y0) * (1.0f - tx) +
                     environmentTexel(x0 + 1, y0) * tx;
    const float3 b = environmentTexel(x0, y0 + 1) * (1.0f - tx) +
                     environmentTexel(x0 + 1, y0 + 1) * tx;
    return (a * (1.0f - ty) + b * ty) * params.environmentStrength;
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

static __forceinline__ __device__ float roughnessAlpha(float roughness) {
    return fmaxf(roughness * roughness, 1e-4f);
}

// Heitz, "Sampling the GGX Distribution of Visible Normals". Sampling the
// visible lobe rather than the full NDF keeps the microfacet in front of the
// viewer and reduces the weight to the masking ratio below.
static __forceinline__ __device__ float3 ggxNormal(float3 view, float3 normal,
                                                   float alpha, unsigned int& rng) {
    const float3 tangent = normalize(fabsf(normal.x) > 0.5f
        ? cross(make_float3(0.0f, 1.0f, 0.0f), normal)
        : cross(make_float3(1.0f, 0.0f, 0.0f), normal));
    const float3 bitangent = cross(normal, tangent);
    const float3 local = make_float3(dot(view, tangent), dot(view, bitangent),
                                     dot(view, normal));

    const float3 stretched =
        normalize(make_float3(alpha * local.x, alpha * local.y, local.z));
    const float lengthSquared =
        stretched.x * stretched.x + stretched.y * stretched.y;
    const float3 basisX = lengthSquared > 0.0f
        ? make_float3(-stretched.y, stretched.x, 0.0f) / sqrtf(lengthSquared)
        : make_float3(1.0f, 0.0f, 0.0f);
    const float3 basisY = cross(stretched, basisX);

    const float radius = sqrtf(random(rng));
    const float phi = 6.2831853f * random(rng);
    const float x = radius * cosf(phi);
    float y = radius * sinf(phi);
    const float lerp = 0.5f * (1.0f + stretched.z);
    y = (1.0f - lerp) * sqrtf(fmaxf(0.0f, 1.0f - x * x)) + lerp * y;

    const float3 hemisphere = basisX * x + basisY * y +
        stretched * sqrtf(fmaxf(0.0f, 1.0f - x * x - y * y));
    const float3 micro = normalize(make_float3(
        alpha * hemisphere.x, alpha * hemisphere.y, fmaxf(0.0f, hemisphere.z)));
    return normalize(tangent * micro.x + bitangent * micro.y + normal * micro.z);
}

static __forceinline__ __device__ float smithLambda(float cosine, float alpha) {
    const float squared = cosine * cosine;
    const float tangentSquared = fmaxf(1.0f - squared, 0.0f) / fmaxf(squared, 1e-8f);
    return 0.5f * (sqrtf(1.0f + alpha * alpha * tangentSquared) - 1.0f);
}

// G2 / G1, the only term left in the throughput once the visible normals have
// been sampled. It is what puts the shadowing back into rough reflections.
static __forceinline__ __device__ float maskingRatio(float viewCosine,
                                                     float lightCosine, float alpha) {
    const float view = smithLambda(viewCosine, alpha);
    return (1.0f + view) / (1.0f + view + smithLambda(lightCosine, alpha));
}

static __forceinline__ __device__ float fresnel(float cosine, float ior) {
    const float r = (1.0f - ior) / (1.0f + ior);
    const float r2 = r * r;
    return r2 + (1.0f - r2) * powf(fmaxf(1.0f - cosine, 0.0f), 5.0f);
}

static __forceinline__ __device__ float3 schlick(float3 f0, float cosine) {
    const float scale = powf(fmaxf(1.0f - cosine, 0.0f), 5.0f);
    return make_float3(f0.x + (1.0f - f0.x) * scale,
                       f0.y + (1.0f - f0.y) * scale,
                       f0.z + (1.0f - f0.z) * scale);
}

static __forceinline__ __device__ float3 baseReflectance(
    const GpuMaterial& material) {
    const float3 color = rgb(material.baseColor);
    return make_float3(0.04f + (color.x - 0.04f) * material.metallic,
                       0.04f + (color.y - 0.04f) * material.metallic,
                       0.04f + (color.z - 0.04f) * material.metallic);
}

static __forceinline__ __device__ float3 absorption(const GpuMaterial& material,
                                                     float distance) {
    if (!isfinite(material.attenuationDistance) ||
        material.attenuationDistance <= 0.0f)
        return make_float3(1.0f, 1.0f, 1.0f);
    const float scale = distance / material.attenuationDistance;
    return make_float3(
        powf(fmaxf(material.attenuationColor.x, 1e-4f), scale),
        powf(fmaxf(material.attenuationColor.y, 1e-4f), scale),
        powf(fmaxf(material.attenuationColor.z, 1e-4f), scale));
}

static __forceinline__ __device__ float3 lensSample(unsigned int& rng) {
    const float radius = sqrtf(random(rng));
    const float angle = 6.2831853f * random(rng);
    return params.lensU * (radius * cosf(angle) * params.aperture) +
           params.lensV * (radius * sinf(angle) * params.aperture);
}

static __forceinline__ __device__ float3 aces(float3 color) {
    color = color * (2.51f * color + 0.03f) /
            (color * (2.43f * color + 0.59f) + 0.14f);
    return make_float3(
        powf(fminf(fmaxf(color.x, 0.0f), 1.0f), 1.0f / 2.2f),
        powf(fminf(fmaxf(color.y, 0.0f), 1.0f), 1.0f / 2.2f),
        powf(fminf(fmaxf(color.z, 0.0f), 1.0f), 1.0f / 2.2f));
}

static __forceinline__ __device__ Hit trace(float3 origin, float3 direction,
                                            float maximumDistance = 1e16f) {
    Hit hit{};
    unsigned int high;
    unsigned int low;
    packPointer(&hit, high, low);
    optixTrace(params.scene, origin, direction, 0.001f, maximumDistance, 0.0f, 255,
               OPTIX_RAY_FLAG_DISABLE_ANYHIT, 0, 1, 0, high, low);
    return hit;
}

static __forceinline__ __device__ float powerHeuristic(float a, float b) {
    const float a2 = a * a;
    const float b2 = b * b;
    return a2 / fmaxf(a2 + b2, 1e-12f);
}

// Tracking the view-dependent Fresnel keeps the lobe choice proportional to the
// energy each lobe actually carries, so grazing angles do not turn into a rare
// draw weighted by a huge factor.
static __forceinline__ __device__ float specularProbability(
    const GpuMaterial& material, float3 view, float3 normal) {
    const float3 reflectance =
        schlick(baseReflectance(material), fmaxf(dot(view, normal), 0.0f));
    return fminf(fmaxf(maximum(reflectance), 0.05f), 0.95f);
}

struct LightSample {
    float3 direction;
    float3 radiance;
    float distance;
    float pdf;
    unsigned int instance;
    unsigned int primitive;
    bool delta;
    bool environment;
    bool valid;
};

struct EnvironmentSample {
    float3 direction;
    float3 radiance;
    float pdf;
};

static __forceinline__ __device__ const GpuLight* selectLight(
    unsigned int& rng, float& probability) {
    if (params.lightCount == 0U || params.lightWeight <= 0.0f)
        return nullptr;
    float target = random(rng) * params.lightWeight;
    const GpuLight* selected = &params.lights[params.lightCount - 1U];
    for (unsigned int index = 0; index < params.lightCount; ++index) {
        selected = &params.lights[index];
        target -= selected->weight;
        if (target <= 0.0f)
            break;
    }
    probability = selected->weight / params.lightWeight;
    return selected;
}

static __forceinline__ __device__ EnvironmentSample sampleEnvironment(
    unsigned int& rng) {
    const unsigned int count =
        params.environmentWidth * params.environmentHeight;
    const float target = random(rng);
    unsigned int first = 0;
    unsigned int last = count - 1U;
    while (first < last) {
        const unsigned int middle = (first + last) / 2U;
        if (params.environmentCdf[middle] < target)
            first = middle + 1U;
        else
            last = middle;
    }
    const unsigned int index = first;
    const float probability = params.environmentCdf[index] -
        (index ? params.environmentCdf[index - 1U] : 0.0f);
    const float u = (static_cast<float>(index % params.environmentWidth) +
                     random(rng)) / params.environmentWidth;
    const float v = (static_cast<float>(index / params.environmentWidth) +
                     random(rng)) / params.environmentHeight;
    const float theta = 3.14159265f * v;
    const float phi = 6.2831853f * (u - 0.5f) - params.environmentRotation;
    const float sine = sinf(theta);
    const float3 direction =
        make_float3(cosf(phi) * sine, cosf(theta), sinf(phi) * sine);
    const float solidAngle = 19.7392088f * fmaxf(sine, 1e-6f) /
        static_cast<float>(count);
    return {direction, environment(direction), probability / solidAngle};
}

static __forceinline__ __device__ LightSample sampleLight(float3 position,
                                                          unsigned int& rng) {
    LightSample sample{};
    float choice = 0.0f;
    const GpuLight* selected = selectLight(rng, choice);
    if (!selected)
        return sample;

    if (selected->type == 4U) {
        const auto source = sampleEnvironment(rng);
        sample.direction = source.direction;
        sample.radiance = source.radiance;
        sample.distance = 1e16f;
        sample.pdf = choice * source.pdf;
        sample.environment = true;
        sample.valid = true;
        return sample;
    }

    if (selected->type == 3U) {
        const float root = sqrtf(random(rng));
        const float u = 1.0f - root;
        const float v = random(rng) * root;
        const float3 point = selected->a * u + selected->b * v +
                             selected->c * (1.0f - u - v);
        const float3 offset = point - position;
        const float distanceSquared = dot(offset, offset);
        sample.distance = sqrtf(distanceSquared);
        sample.direction = offset / sample.distance;
        const float cosine = fabsf(dot(selected->normal, -sample.direction));
        if (cosine <= 1e-6f)
            return sample;
        sample.pdf = choice * distanceSquared / (selected->area * cosine);
        sample.radiance = selected->emission;
        sample.instance = selected->instance;
        sample.primitive = selected->primitive;
        sample.valid = true;
        return sample;
    }

    sample.delta = true;
    sample.pdf = choice;
    sample.radiance = selected->emission;
    if (selected->type == 0U) {
        sample.direction = -selected->b;
        sample.distance = 1e16f;
    } else {
        const float3 offset = selected->a - position;
        const float distanceSquared = dot(offset, offset);
        sample.distance = sqrtf(distanceSquared);
        if (sample.distance >= selected->range)
            return sample;
        sample.direction = offset / sample.distance;
        sample.radiance /= distanceSquared;
        if (selected->type == 1U) {
            const float cone = dot(selected->b, -sample.direction);
            const float outer = cosf(selected->outerCone);
            const float inner = cosf(selected->innerCone);
            const float falloff = saturate((cone - outer) / fmaxf(inner - outer, 1e-5f));
            sample.radiance = sample.radiance * (falloff * falloff);
            if (falloff <= 0.0f)
                return sample;
        }
    }
    sample.valid = true;
    return sample;
}

static __forceinline__ __device__ bool visible(float3 position, float3 normal,
                                                const LightSample& light) {
    const float limit = light.delta ? light.distance - 0.002f
                                    : light.distance + 0.002f;
    const Hit blocker = trace(position + normal * 0.001f, light.direction, limit);
    if (light.environment)
        return !blocker.found;
    if (light.delta)
        return !blocker.found;
    return blocker.found && blocker.instance == light.instance &&
           blocker.primitive == light.primitive;
}

static __forceinline__ __device__ float environmentPdf(float3 direction) {
    if (!params.environment || params.lightWeight <= 0.0f)
        return 0.0f;
    const float2 uv = environmentUv(direction);
    const unsigned int x = min(static_cast<unsigned int>(
        uv.x * params.environmentWidth), params.environmentWidth - 1U);
    const unsigned int y = min(static_cast<unsigned int>(
        uv.y * params.environmentHeight), params.environmentHeight - 1U);
    const unsigned int index = y * params.environmentWidth + x;
    const float probability = params.environmentCdf[index] -
        (index ? params.environmentCdf[index - 1U] : 0.0f);
    const float theta = 3.14159265f *
        (static_cast<float>(y) + 0.5f) / params.environmentHeight;
    const float solidAngle = 19.7392088f * fmaxf(sinf(theta), 1e-6f) /
        static_cast<float>(params.environmentWidth * params.environmentHeight);
    return params.environmentWeight / params.lightWeight *
           probability / solidAngle;
}

static __forceinline__ __device__ float3 directLighting(
    const Hit& hit, const GpuMaterial& material, float3 view, unsigned int& rng) {
    const float diffuseWeight = (1.0f - material.metallic) *
                                (1.0f - material.transmission);
    if (diffuseWeight <= 0.0f)
        return make_float3(0.0f, 0.0f, 0.0f);

    const LightSample light = sampleLight(hit.position, rng);
    const float cosine = dot(hit.normal, light.direction);
    if (!light.valid || cosine <= 0.0f || !visible(hit.position, hit.normal, light))
        return make_float3(0.0f, 0.0f, 0.0f);

    // The same (1 - F) split the path tracer applies when it picks the diffuse
    // lobe, so next event estimation and BSDF sampling agree and MIS stays sound.
    const float3 reflectance =
        schlick(baseReflectance(material), fmaxf(dot(view, hit.normal), 0.0f));
    const float3 diffuse = make_float3(1.0f - reflectance.x, 1.0f - reflectance.y,
                                       1.0f - reflectance.z);
    const float3 bsdf =
        rgb(material.baseColor) * diffuse * (diffuseWeight / 3.14159265f);
    const float bsdfPdf = (1.0f - material.transmission) *
        (1.0f - specularProbability(material, view, hit.normal)) * cosine /
        3.14159265f;
    const float weight = light.delta ? 1.0f : powerHeuristic(light.pdf, bsdfPdf);
    return bsdf * light.radiance * (cosine * weight / light.pdf);
}

static __forceinline__ __device__ float emissivePdf(float3 origin, const Hit& hit) {
    const float3 offset = hit.position - origin;
    const float distanceSquared = dot(offset, offset);
    const float3 direction = offset / sqrtf(distanceSquared);
    for (unsigned int index = 0; index < params.lightCount; ++index) {
        const GpuLight light = params.lights[index];
        if (light.type != 3U || light.instance != hit.instance ||
            light.primitive != hit.primitive)
            continue;
        const float cosine = fabsf(dot(light.normal, -direction));
        return light.weight / params.lightWeight * distanceSquared /
               fmaxf(light.area * cosine, 1e-8f);
    }
    return 0.0f;
}

struct PhotonEmission {
    float3 origin;
    float3 direction;
    float3 power;
    bool valid;
};

static __forceinline__ __device__ float3 directionAround(
    float3 axis, float cosine, unsigned int& rng) {
    const float sine = sqrtf(fmaxf(0.0f, 1.0f - cosine * cosine));
    const float angle = 6.2831853f * random(rng);
    const float3 tangent = normalize(fabsf(axis.x) > 0.5f
        ? cross(make_float3(0.0f, 1.0f, 0.0f), axis)
        : cross(make_float3(1.0f, 0.0f, 0.0f), axis));
    return normalize(axis * cosine + tangent * (sine * cosf(angle)) +
                     cross(axis, tangent) * (sine * sinf(angle)));
}

static __forceinline__ __device__ float3 diskOffset(
    float3 normal, float radius, unsigned int& rng) {
    const float distance = sqrtf(random(rng)) * radius;
    const float angle = 6.2831853f * random(rng);
    const float3 tangent = normalize(fabsf(normal.x) > 0.5f
        ? cross(make_float3(0.0f, 1.0f, 0.0f), normal)
        : cross(make_float3(1.0f, 0.0f, 0.0f), normal));
    return tangent * (distance * cosf(angle)) +
           cross(normal, tangent) * (distance * sinf(angle));
}

static __forceinline__ __device__ PhotonEmission emitPhoton(
    unsigned int& rng) {
    PhotonEmission photon{};
    float choice = 0.0f;
    const GpuLight* light = selectLight(rng, choice);
    if (!light)
        return photon;
    const float normalization =
        1.0f / (static_cast<float>(params.photonEmissions) * choice);

    if (light->type == 4U) {
        const auto source = sampleEnvironment(rng);
        photon.direction = -source.direction;
        photon.origin = params.sceneCenter -
            photon.direction * params.sceneRadius +
            diskOffset(photon.direction, params.sceneRadius, rng);
        photon.power = source.radiance *
            (3.14159265f * params.sceneRadius * params.sceneRadius *
             normalization / source.pdf);
    } else if (light->type == 3U) {
        const float root = sqrtf(random(rng));
        const float u = 1.0f - root;
        const float v = random(rng) * root;
        photon.origin = light->a * u + light->b * v +
                        light->c * (1.0f - u - v);
        const float3 side = random(rng) < 0.5f ? light->normal : -light->normal;
        photon.direction = cosineDirection(side, rng);
        photon.origin = photon.origin + photon.direction * 0.001f;
        photon.power = light->emission *
            (6.2831853f * light->area * normalization);
    } else if (light->type == 0U) {
        photon.direction = normalize(light->b);
        photon.origin = params.sceneCenter -
            photon.direction * params.sceneRadius +
            diskOffset(photon.direction, params.sceneRadius, rng);
        photon.power = light->emission *
            (3.14159265f * params.sceneRadius * params.sceneRadius *
             normalization);
    } else {
        photon.origin = light->a;
        if (light->type == 1U) {
            const float outer = cosf(light->outerCone);
            const float cosine = 1.0f - random(rng) * (1.0f - outer);
            photon.direction = directionAround(normalize(light->b), cosine, rng);
            const float inner = cosf(light->innerCone);
            const float falloff =
                saturate((cosine - outer) / fmaxf(inner - outer, 1e-5f));
            photon.power = light->emission *
                (6.2831853f * (1.0f - outer) * falloff * falloff *
                 normalization);
        } else {
            const float cosine = 1.0f - 2.0f * random(rng);
            photon.direction =
                directionAround(make_float3(0.0f, 1.0f, 0.0f), cosine, rng);
            photon.power = light->emission * (12.566371f * normalization);
        }
        photon.origin = photon.origin + photon.direction * 0.001f;
    }
    photon.valid = maximum(photon.power) > 0.0f;
    return photon;
}

static __forceinline__ __device__ void storePhoton(
    const Hit& hit, float3 power) {
    const unsigned int slot = atomicAdd(params.photonCount, 1U);
    if (slot >= params.photonCapacity)
        return;
    params.photons[slot] = {hit.position, power, hit.normal};
}

static __forceinline__ __device__ float3 causticLighting(
    const Hit& hit, const GpuMaterial& material) {
    const float diffuseWeight = (1.0f - material.metallic) *
                                (1.0f - material.transmission);
    if (params.photonRadius <= 0.0f || diffuseWeight <= 0.0f)
        return make_float3(0.0f, 0.0f, 0.0f);

    const int3 center = photonCell(hit.position, params.photonRadius);
    const float radiusSquared = params.photonRadius * params.photonRadius;
    float3 flux = make_float3(0.0f, 0.0f, 0.0f);
    for (int z = -1; z <= 1; ++z) {
        for (int y = -1; y <= 1; ++y) {
            for (int x = -1; x <= 1; ++x) {
                const unsigned int bucket = photonBucket(
                    make_int3(center.x + x, center.y + y, center.z + z),
                    params.photonBucketCount);
                const unsigned int end = params.photonCellStart[bucket + 1U];
                for (unsigned int index = params.photonCellStart[bucket];
                     index < end; ++index) {
                    const Photon photon = params.sortedPhotons[index];
                    const float3 offset = photon.position - hit.position;
                    const float distanceSquared = dot(offset, offset);
                    if (distanceSquared >= radiusSquared ||
                        dot(photon.normal, hit.normal) < 0.7f)
                        continue;
                    flux += photon.power *
                        (1.0f - distanceSquared / radiusSquared);
                }
            }
        }
    }
    const float kernel = 2.0f /
        (3.14159265f * radiusSquared * 3.14159265f);
    return rgb(material.baseColor) * (diffuseWeight * kernel) * flux;
}

extern "C" __global__ void __raygen__photons() {
    const unsigned int index = optixGetLaunchIndex().x;
    // Every pass reseeds so that repeated builds explore different photon paths
    // and the caustic keeps converging as camera samples accumulate.
    unsigned int rng = seeded(index, params.photonPass);
    const PhotonEmission emission = emitPhoton(rng);
    if (!emission.valid)
        return;

    float3 origin = emission.origin;
    float3 direction = emission.direction;
    float3 power = emission.power;
    int medium = -1;
    bool specularPath = false;
    for (unsigned int depth = 0; depth < 10U; ++depth) {
        Hit hit = trace(origin, direction);
        if (!hit.found)
            return;
        GpuMaterial material = textured(params.materials[hit.material], hit);
        hit.normal = mappedNormal(hit, material);
        if (medium >= 0)
            power *= absorption(params.materials[medium], hit.distance);

        const float diffuseWeight = (1.0f - material.metallic) *
                                    (1.0f - material.transmission);
        if (specularPath && diffuseWeight > 0.0f) {
            storePhoton(hit, power);
            return;
        }

        if (material.transmission > 0.0f &&
            random(rng) < material.transmission) {
            float ior = material.ior;
            if (material.dispersion > 0.0f) {
                const unsigned int channel =
                    static_cast<unsigned int>(random(rng) * 3.0f) % 3U;
                const float spread =
                    (material.ior - 1.0f) * material.dispersion * 0.5f;
                ior += (static_cast<float>(channel) - 1.0f) * spread;
                power *= channel == 0U ? make_float3(3.0f, 0.0f, 0.0f)
                    : channel == 1U ? make_float3(0.0f, 3.0f, 0.0f)
                                    : make_float3(0.0f, 0.0f, 3.0f);
            }
            const float3 view = -direction;
            const float alpha = roughnessAlpha(material.roughness);
            const float3 microNormal = ggxNormal(view, hit.normal, alpha, rng);
            const float eta = hit.frontFace ? 1.0f / ior : ior;
            const float cosine = fminf(fmaxf(dot(view, microNormal), 0.0f), 1.0f);
            float3 transmitted;
            const bool totalReflection =
                !refract(direction, microNormal, eta, transmitted);
            float reflectance = 1.0f;
            if (!totalReflection) {
                transmitted = normalize(transmitted);
                reflectance = fresnel(
                    eta > 1.0f ? fabsf(dot(transmitted, microNormal)) : cosine, ior);
            }

            const float viewCosine = fabsf(dot(view, hit.normal));
            if (totalReflection || random(rng) < reflectance) {
                direction = reflect(direction, microNormal);
                if (dot(direction, hit.normal) <= 0.0f)
                    return;
                origin = hit.position + hit.normal * 0.001f;
            } else {
                direction = transmitted;
                origin = hit.position - hit.normal * 0.001f;
                medium = hit.frontFace ? static_cast<int>(hit.material) : -1;
            }
            power *= maskingRatio(viewCosine,
                                  fabsf(dot(direction, hit.normal)), alpha);
            specularPath = true;
            continue;
        }
        return;
    }
}

extern "C" __global__ void __raygen__render() {
    const uint3 pixel = optixGetLaunchIndex();
    const unsigned int index = pixel.y * params.width + pixel.x;
    if (params.display) {
        float3 color = rgb(params.display[index]);
        if (params.display2)
            color += rgb(params.display2[index]);
        if (params.display3)
            color += rgb(params.display3[index]);
        if (params.composite)
            params.composite[index] = make_float4(color.x, color.y, color.z, 1.0f);
        const float3 mapped = aces(color * exp2f(params.exposure));
        params.output[index] = make_uchar4(
            static_cast<unsigned char>(mapped.x * 255.0f + 0.5f),
            static_cast<unsigned char>(mapped.y * 255.0f + 0.5f),
            static_cast<unsigned char>(mapped.z * 255.0f + 0.5f), 255);
        return;
    }
    unsigned int rng = index * 9781U + params.sample * 6271U + 0x68bc21ebU;
    const float2 jitter = make_float2(random(rng), random(rng));
    const float2 screen = make_float2(
        (static_cast<float>(pixel.x) + jitter.x) / static_cast<float>(params.width),
        (static_cast<float>(pixel.y) + jitter.y) / static_cast<float>(params.height));

    const float3 pinhole = normalize(params.cameraW +
        (2.0f * screen.x - 1.0f) * params.cameraU +
        (2.0f * screen.y - 1.0f) * params.cameraV);
    const float3 focalPoint = params.eye + pinhole *
        (params.focusDistance / fmaxf(dot(pinhole, params.cameraW), 1e-4f));
    float3 origin = params.eye + lensSample(rng);
    float3 direction = normalize(focalPoint - origin);
    float3 throughput = make_float3(1.0f, 1.0f, 1.0f);
    float3 radiance[3] = {};
    float3 guideAlbedo = make_float3(0.0f, 0.0f, 0.0f);
    float3 guideNormal = make_float3(0.0f, 0.0f, 0.0f);
    unsigned int lobe = 0;
    int medium = -1;
    float lastPdf = 0.0f;
    float3 lastOrigin = origin;
    bool lastDelta = true;
    bool primaryChain = true;
    bool guidePending = true;

    for (unsigned int depth = 0; depth < params.maxDepth; ++depth) {
        Hit hit = trace(origin, direction);
        if (!hit.found) {
            const float pdf = lastDelta ? 0.0f : environmentPdf(direction);
            const float weight = lastDelta ? 1.0f : powerHeuristic(lastPdf, pdf);
            radiance[lobe] += throughput * environment(direction) * weight;
            break;
        }

        const GpuMaterial material = textured(params.materials[hit.material], hit);
        hit.normal = mappedNormal(hit, material);
        const float3 view = -direction;
        if (depth == 0U) {
            // Classify the pixel once, from the surface the camera actually sees.
            // Drawing the lobe per sample instead would scatter a pixel's energy
            // across layers and leave each one far noisier than the beauty pass.
            lobe = material.transmission > 0.5f ? 2U
                 : material.metallic > 0.5f ? 1U : 0U;
            // The denoiser wants camera space, and taking the normal from the
            // primary hit is what gives it the silhouettes to stop blurring at.
            guideNormal = make_float3(dot(hit.normal, params.lensU),
                                      dot(hit.normal, params.lensV),
                                      dot(hit.normal, params.cameraW));
        }
        if (guidePending) {
            const bool clearGlass =
                material.transmission > 0.95f && material.roughness < 0.1f;
            if (!clearGlass) {
                if (material.transmission < 0.05f)
                    guideAlbedo = rgb(material.baseColor);
                guidePending = false;
            }
        }
        if (medium >= 0)
            throughput *= absorption(params.materials[medium], hit.distance);
        const float lightPdf = lastDelta ? 0.0f : emissivePdf(lastOrigin, hit);
        const float emissionWeight = lastDelta ? 1.0f
            : powerHeuristic(lastPdf, lightPdf);
        radiance[lobe] += throughput * material.emissive *
                          (material.emissiveStrength * emissionWeight);
        radiance[lobe] += throughput * directLighting(hit, material, view, rng);
        if (primaryChain)
            radiance[lobe] += throughput * causticLighting(hit, material);

        const bool transmissive = material.transmission > 0.0f &&
                                  random(rng) < material.transmission;
        if (transmissive) {
            lastDelta = true;
            float ior = material.ior;
            if (material.dispersion > 0.0f) {
                const unsigned int channel =
                    static_cast<unsigned int>(random(rng) * 3.0f) % 3U;
                const float spread = (material.ior - 1.0f) * material.dispersion * 0.5f;
                ior += (static_cast<float>(channel) - 1.0f) * spread;
                throughput *= channel == 0U ? make_float3(3.0f, 0.0f, 0.0f)
                    : channel == 1U ? make_float3(0.0f, 3.0f, 0.0f)
                                    : make_float3(0.0f, 0.0f, 3.0f);
            }

            const float alpha = roughnessAlpha(material.roughness);
            const float3 microNormal = ggxNormal(view, hit.normal, alpha, rng);
            const float eta = hit.frontFace ? 1.0f / ior : ior;
            const float cosine = fminf(fmaxf(dot(view, microNormal), 0.0f), 1.0f);
            float3 transmitted;
            const bool totalReflection =
                !refract(direction, microNormal, eta, transmitted);

            // Schlick is only valid on the dense side of the interface, so when
            // the ray is leaving the glass the transmitted angle is the one to
            // feed it.
            float reflectance = 1.0f;
            if (!totalReflection) {
                transmitted = normalize(transmitted);
                reflectance = fresnel(
                    eta > 1.0f ? fabsf(dot(transmitted, microNormal)) : cosine, ior);
            }

            const float viewCosine = fabsf(dot(view, hit.normal));
            if (totalReflection || random(rng) < reflectance) {
                direction = reflect(direction, microNormal);
                if (dot(direction, hit.normal) <= 0.0f)
                    break;
                origin = hit.position + hit.normal * 0.001f;
            } else {
                direction = transmitted;
                origin = hit.position - hit.normal * 0.001f;
                medium = hit.frontFace ? static_cast<int>(hit.material) : -1;
            }
            throughput *= maskingRatio(viewCosine,
                                       fabsf(dot(direction, hit.normal)), alpha);
        } else {
            const float3 color = rgb(material.baseColor);
            const float3 f0 = baseReflectance(material);
            const float viewCosine = fmaxf(dot(view, hit.normal), 0.0f);
            const float probability =
                specularProbability(material, view, hit.normal);
            if (random(rng) < probability) {
                lastDelta = true;
                const float alpha = roughnessAlpha(material.roughness);
                const float3 microNormal = ggxNormal(view, hit.normal, alpha, rng);
                direction = reflect(direction, microNormal);
                const float lightCosine = dot(direction, hit.normal);
                if (lightCosine <= 0.0f)
                    break;
                // Fresnel at the microfacet, not at normal incidence: this is
                // what makes grazing angles approach a mirror instead of a flat
                // four percent.
                throughput *= schlick(f0, fmaxf(dot(view, microNormal), 0.0f)) *
                    (maskingRatio(viewCosine, lightCosine, alpha) / probability);
            } else {
                lastDelta = false;
                primaryChain = false;
                direction = cosineDirection(hit.normal, rng);
                const float3 reflectance = schlick(f0, viewCosine);
                const float3 diffuse =
                    make_float3(1.0f - reflectance.x, 1.0f - reflectance.y,
                                1.0f - reflectance.z);
                throughput *= color * diffuse * ((1.0f - material.metallic) /
                                                 (1.0f - probability));
                lastPdf = (1.0f - material.transmission) *
                          (1.0f - probability) *
                          fmaxf(dot(hit.normal, direction), 0.0f) / 3.14159265f;
            }
            origin = hit.position + hit.normal * 0.001f;
        }
        lastOrigin = hit.position;

        if (depth > 2) {
            const float survival = fminf(fmaxf(maximum(throughput), 0.05f), 0.95f);
            if (random(rng) > survival)
                break;
            throughput /= survival;
        }
    }

    const float4 previous = params.accumulation[index];
    const float weight = 1.0f / static_cast<float>(params.sample + 1U);
    const float3 sample = radiance[0] + radiance[1] + radiance[2];
    const float3 accumulated = rgb(previous) + (sample - rgb(previous)) * weight;
    const float3 albedo = rgb(params.albedoGuide[index]) +
        (guideAlbedo - rgb(params.albedoGuide[index])) * weight;
    const float3 normal = rgb(params.normalGuide[index]) +
        (guideNormal - rgb(params.normalGuide[index])) * weight;
    params.accumulation[index] =
        make_float4(accumulated.x, accumulated.y, accumulated.z, 1.0f);
    const float3 diffuse = rgb(params.diffuse[index]) +
        (radiance[0] - rgb(params.diffuse[index])) * weight;
    const float3 reflection = rgb(params.reflection[index]) +
        (radiance[1] - rgb(params.reflection[index])) * weight;
    const float3 refraction = rgb(params.refraction[index]) +
        (radiance[2] - rgb(params.refraction[index])) * weight;
    params.diffuse[index] = make_float4(diffuse.x, diffuse.y, diffuse.z, 1.0f);
    params.reflection[index] =
        make_float4(reflection.x, reflection.y, reflection.z, 1.0f);
    params.refraction[index] =
        make_float4(refraction.x, refraction.y, refraction.z, 1.0f);
    params.albedoGuide[index] = make_float4(albedo.x, albedo.y, albedo.z, 1.0f);
    params.normalGuide[index] = make_float4(normal.x, normal.y, normal.z, 1.0f);
    const float3 mapped = aces(accumulated * exp2f(params.exposure));
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
    const GpuVertex v0 = data->vertices[triangle.x];
    const GpuVertex v1 = data->vertices[triangle.y];
    const GpuVertex v2 = data->vertices[triangle.z];
    float3 normal = v0.normal * b0 + v1.normal * barycentrics.x +
                    v2.normal * barycentrics.y;
    if (dot(normal, normal) < 1e-12f) {
        normal = cross(v1.position - v0.position, v2.position - v0.position);
    }

    const float3 objectNormal = normalize(normal);
    normal = normalize(optixTransformNormalFromObjectToWorldSpace(normal));
    const float3 ray = optixGetWorldRayDirection();
    hit->frontFace = dot(normal, ray) < 0.0f;
    hit->normal = hit->frontFace ? normal : -normal;
    float4 tangent = make_float4(
        v0.tangent.x * b0 + v1.tangent.x * barycentrics.x +
            v2.tangent.x * barycentrics.y,
        v0.tangent.y * b0 + v1.tangent.y * barycentrics.x +
            v2.tangent.y * barycentrics.y,
        v0.tangent.z * b0 + v1.tangent.z * barycentrics.x +
            v2.tangent.z * barycentrics.y,
        v0.tangent.w * b0 + v1.tangent.w * barycentrics.x +
            v2.tangent.w * barycentrics.y);
    float3 tangentDirection = make_float3(tangent.x, tangent.y, tangent.z);
    if (dot(tangentDirection, tangentDirection) < 1e-12f) {
        const float2 duv1 = make_float2(v1.uv.x - v0.uv.x, v1.uv.y - v0.uv.y);
        const float2 duv2 = make_float2(v2.uv.x - v0.uv.x, v2.uv.y - v0.uv.y);
        const float determinant = duv1.x * duv2.y - duv1.y * duv2.x;
        tangentDirection = fabsf(determinant) > 1e-8f
            ? ((v1.position - v0.position) * duv2.y -
               (v2.position - v0.position) * duv1.y) / determinant
            : cross(fabsf(objectNormal.x) > 0.5f
                        ? make_float3(0.0f, 1.0f, 0.0f)
                        : make_float3(1.0f, 0.0f, 0.0f),
                    objectNormal);
        tangent.w = 1.0f;
    }
    tangentDirection =
        optixTransformVectorFromObjectToWorldSpace(tangentDirection);
    tangentDirection = normalize(
        tangentDirection - hit->normal * dot(tangentDirection, hit->normal));
    hit->tangent =
        make_float4(tangentDirection.x, tangentDirection.y, tangentDirection.z,
                    tangent.w < 0.0f ? -1.0f : 1.0f);
    hit->uv = make_float2(v0.uv.x * b0 + v1.uv.x * barycentrics.x +
                              v2.uv.x * barycentrics.y,
                          v0.uv.y * b0 + v1.uv.y * barycentrics.x +
                              v2.uv.y * barycentrics.y);
    hit->uv1 = make_float2(v0.uv1.x * b0 + v1.uv1.x * barycentrics.x +
                               v2.uv1.x * barycentrics.y,
                           v0.uv1.y * b0 + v1.uv1.y * barycentrics.x +
                               v2.uv1.y * barycentrics.y);
    hit->position = optixGetWorldRayOrigin() + optixGetRayTmax() * ray;
    hit->distance = optixGetRayTmax();
    hit->material = data->material;
    hit->instance = optixGetInstanceId();
    hit->primitive = optixGetPrimitiveIndex();
    hit->found = true;
}

}
