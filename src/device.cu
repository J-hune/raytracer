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

    for (unsigned int depth = 0; depth < 8; ++depth) {
        const Hit hit = trace(origin, direction);
        if (!hit.found) {
            const float pdf = lastDelta ? 0.0f : environmentPdf(direction);
            const float weight = lastDelta ? 1.0f : powerHeuristic(lastPdf, pdf);
            radiance += throughput * environment(direction) * weight;
            break;
        }

        const GpuMaterial material = params.materials[hit.material];
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
    hit->frontFace = dot(normal, ray) < 0.0f;
    hit->normal = hit->frontFace ? normal : -normal;
    hit->position = optixGetWorldRayOrigin() + optixGetRayTmax() * ray;
    hit->distance = optixGetRayTmax();
    hit->material = data->material;
    hit->instance = optixGetInstanceId();
    hit->primitive = optixGetPrimitiveIndex();
    hit->found = true;
}

}
