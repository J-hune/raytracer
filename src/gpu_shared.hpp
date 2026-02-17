#pragma once

#include <cuda_runtime.h>
#include <optix.h>

namespace rt {

struct GpuVertex {
    float3 position;
    float3 normal;
    float4 tangent;
    float2 uv;
    float2 uv1;
};

struct GpuTextureRef {
    int texture;
    unsigned int texCoord;
    float2 offset;
    float2 scale;
    float rotation;
    float strength;
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
    GpuTextureRef baseColorTexture;
    GpuTextureRef metallicRoughnessTexture;
    GpuTextureRef normalTexture;
    GpuTextureRef emissiveTexture;
    GpuTextureRef transmissionTexture;
    GpuTextureRef thicknessTexture;
};

struct GpuTexture {
    const uchar4* pixels;
    unsigned int width;
    unsigned int height;
    unsigned int wrapU;
    unsigned int wrapV;
};

struct Hit {
    float3 position;
    float3 normal;
    float4 tangent;
    float2 uv;
    float2 uv1;
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

struct Photon {
    float3 position;
    float3 power;
    float3 normal;
};

struct LaunchParams {
    uchar4* output;
    float4* accumulation;
    float4* diffuse;
    float4* reflection;
    float4* refraction;
    float4* albedoGuide;
    float4* normalGuide;
    const float4* display;
    const float4* display2;
    const float4* display3;
    float4* composite;
    const GpuMaterial* materials;
    const GpuTexture* textures;
    const GpuLight* lights;
    const float4* environment;
    const float* environmentCdf;
    Photon* photons;
    unsigned int* photonCount;
    const Photon* sortedPhotons;
    const unsigned int* photonCellStart;
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
    float environmentWeight;
    float environmentRotation;
    float environmentStrength;
    float exposure;
    float3 sceneCenter;
    float sceneRadius;
    float photonRadius;
    unsigned int lightCount;
    unsigned int environmentWidth;
    unsigned int environmentHeight;
    unsigned int width;
    unsigned int height;
    unsigned int sample;
    unsigned int maxDepth;
    unsigned int photonBucketCount;
    unsigned int photonCapacity;
    unsigned int photonEmissions;
    unsigned int photonPass;
};

}
