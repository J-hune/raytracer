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

static __forceinline__ __device__ float3 ggxNormal(float3 normal, float roughness,
                                                   unsigned int& rng) {
    const float alpha = fmaxf(roughness * roughness, 0.001f);
    const float phi = 6.2831853f * random(rng);
    const float sample = random(rng);
    const float cosine = sqrtf((1.0f - sample) /
                               (1.0f + (alpha * alpha - 1.0f) * sample));
    const float sine = sqrtf(fmaxf(0.0f, 1.0f - cosine * cosine));
    const float3 tangent = normalize(fabsf(normal.x) > 0.5f
        ? cross(make_float3(0.0f, 1.0f, 0.0f), normal)
        : cross(make_float3(1.0f, 0.0f, 0.0f), normal));
    return normalize(tangent * (sine * cosf(phi)) +
                     cross(normal, tangent) * (sine * sinf(phi)) + normal * cosine);
}

static __forceinline__ __device__ float fresnel(float cosine, float ior) {
    const float r = (1.0f - ior) / (1.0f + ior);
    const float r2 = r * r;
    return r2 + (1.0f - r2) * powf(1.0f - cosine, 5.0f);
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

static __forceinline__ __device__ float specularProbability(
    const GpuMaterial& material) {
    const float3 color = rgb(material.baseColor);
    return fminf(fmaxf(maximum(make_float3(
        0.04f + (color.x - 0.04f) * material.metallic,
        0.04f + (color.y - 0.04f) * material.metallic,
        0.04f + (color.z - 0.04f) * material.metallic)), 0.05f), 0.95f);
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

static __forceinline__ __device__ LightSample sampleLight(float3 position,
                                                          unsigned int& rng) {
    LightSample sample{};
    if (params.lightCount == 0U || params.lightWeight <= 0.0f)
        return sample;

    float target = random(rng) * params.lightWeight;
    const GpuLight* selected = &params.lights[params.lightCount - 1U];
    for (unsigned int index = 0; index < params.lightCount; ++index) {
        selected = &params.lights[index];
        target -= selected->weight;
        if (target <= 0.0f)
            break;
    }
    const float choice = selected->weight / params.lightWeight;

    if (selected->type == 4U) {
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
        sample.direction =
            make_float3(cosf(phi) * sine, cosf(theta), sinf(phi) * sine);
        const float solidAngle = 19.7392088f * fmaxf(sine, 1e-6f) /
            static_cast<float>(count);
        sample.radiance = environment(sample.direction);
        sample.distance = 1e16f;
        sample.pdf = choice * probability / solidAngle;
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
    const Hit& hit, const GpuMaterial& material, unsigned int& rng) {
    const float diffuseWeight = (1.0f - material.metallic) *
                                (1.0f - material.transmission);
    if (diffuseWeight <= 0.0f)
        return make_float3(0.0f, 0.0f, 0.0f);

    const LightSample light = sampleLight(hit.position, rng);
    const float cosine = dot(hit.normal, light.direction);
    if (!light.valid || cosine <= 0.0f || !visible(hit.position, hit.normal, light))
        return make_float3(0.0f, 0.0f, 0.0f);

    const float3 bsdf = rgb(material.baseColor) * (diffuseWeight / 3.14159265f);
    const float bsdfPdf = (1.0f - material.transmission) *
        (1.0f - specularProbability(material)) * cosine / 3.14159265f;
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

extern "C" __global__ void __raygen__render() {
    const uint3 pixel = optixGetLaunchIndex();
    const unsigned int index = pixel.y * params.width + pixel.x;
    if (params.display) {
        const float3 mapped =
            aces(rgb(params.display[index]) * exp2f(params.exposure));
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
    float3 radiance = make_float3(0.0f, 0.0f, 0.0f);
    int medium = -1;
    float lastPdf = 0.0f;
    float3 lastOrigin = origin;
    bool lastDelta = true;

    for (unsigned int depth = 0; depth < params.maxDepth; ++depth) {
        Hit hit = trace(origin, direction);
        if (!hit.found) {
            const float pdf = lastDelta ? 0.0f : environmentPdf(direction);
            const float weight = lastDelta ? 1.0f : powerHeuristic(lastPdf, pdf);
            radiance += throughput * environment(direction) * weight;
            break;
        }

        const GpuMaterial material = textured(params.materials[hit.material], hit);
        hit.normal = mappedNormal(hit, material);
        if (medium >= 0)
            throughput *= absorption(params.materials[medium], hit.distance);
        const float lightPdf = lastDelta ? 0.0f : emissivePdf(lastOrigin, hit);
        const float emissionWeight = lastDelta ? 1.0f
            : powerHeuristic(lastPdf, lightPdf);
        radiance += throughput * material.emissive *
                    (material.emissiveStrength * emissionWeight);
        radiance += throughput * directLighting(hit, material, rng);

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

            float3 microNormal = ggxNormal(hit.normal, material.roughness, rng);
            if (dot(-direction, microNormal) < 0.0f)
                microNormal = -microNormal;
            const float eta = hit.frontFace ? 1.0f / ior : ior;
            const float cosine = fminf(dot(-direction, microNormal), 1.0f);
            float3 transmitted;
            const bool totalReflection = !refract(direction, microNormal, eta, transmitted);
            if (totalReflection || random(rng) < fresnel(cosine, ior)) {
                direction = reflect(direction, microNormal);
                origin = hit.position + hit.normal * 0.001f;
            } else {
                direction = normalize(transmitted);
                origin = hit.position - hit.normal * 0.001f;
                medium = hit.frontFace ? static_cast<int>(hit.material) : -1;
            }
        } else {
            const float3 color = rgb(material.baseColor);
            const float3 f0 = make_float3(
                0.04f + (color.x - 0.04f) * material.metallic,
                0.04f + (color.y - 0.04f) * material.metallic,
                0.04f + (color.z - 0.04f) * material.metallic);
            const float probability = fminf(fmaxf(maximum(f0), 0.05f), 0.95f);
            if (random(rng) < probability) {
                lastDelta = true;
                const float3 microNormal = ggxNormal(hit.normal, material.roughness, rng);
                direction = reflect(direction, microNormal);
                if (dot(direction, hit.normal) <= 0.0f)
                    direction = reflect(direction, hit.normal);
                throughput *= f0 / probability;
            } else {
                lastDelta = false;
                direction = cosineDirection(hit.normal, rng);
                throughput *= color * ((1.0f - material.metallic) /
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
    const float3 accumulated = rgb(previous) + (radiance - rgb(previous)) * weight;
    params.accumulation[index] =
        make_float4(accumulated.x, accumulated.y, accumulated.z, 1.0f);
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
