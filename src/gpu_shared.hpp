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
    float metallic;
    float roughness;
    float transmission;
    float ior;
    float emissiveStrength;
};

struct Hit {
    float3 position;
    float3 normal;
    float2 uv;
    unsigned int material;
    bool found;
};

struct HitGroupData {
    const GpuVertex* vertices;
    const uint3* indices;
    unsigned int material;
};

struct LaunchParams {
    uchar4* output;
    float4* accumulation;
    const GpuMaterial* materials;
    OptixTraversableHandle scene;
    float3 eye;
    float3 cameraU;
    float3 cameraV;
    float3 cameraW;
    unsigned int width;
    unsigned int height;
    unsigned int sample;
};

}
