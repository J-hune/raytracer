#pragma once

#include <cuda_runtime.h>
#include <optix.h>

namespace rt {

struct GpuVertex {
    float3 position;
    float3 normal;
    float4 tangent;
    float2 uv;
};

struct GpuMaterial {
    float4 baseColor;
    float3 emissive;
    float3 attenuationColor;
    float metallic;
    float roughness;
    float transmission;
    float ior;
    float thickness;
    float attenuationDistance;
    float emissiveStrength;
    float dispersion;
};

struct Hit {
    float3 position;
    float3 normal;
    float2 uv;
    float distance;
    unsigned int material;
    unsigned int instance;
    unsigned int primitive;
    bool frontFace;
    bool found;
};

struct HitGroupData {
    const GpuVertex* vertices;
    const uint3* indices;
    unsigned int material;
};

struct GpuLight {
    float3 a;
    float3 b;
    float3 c;
    float3 normal;
    float3 emission;
    float area;
    float weight;
    float range;
    float innerCone;
    float outerCone;
    unsigned int instance;
    unsigned int primitive;
    unsigned int type;
};

struct LaunchParams {
    uchar4* output;
    float4* accumulation;
    const GpuMaterial* materials;
    const GpuLight* lights;
    OptixTraversableHandle scene;
    float3 eye;
    float3 cameraU;
    float3 cameraV;
    float3 cameraW;
    float3 lensU;
    float3 lensV;
    float aperture;
    float focusDistance;
    float lightWeight;
    unsigned int lightCount;
    unsigned int width;
    unsigned int height;
    unsigned int sample;
};

}
