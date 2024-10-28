#include "Lighting.h"
#include "Material.h"
#include <algorithm>

Vec3 Lighting::computePhongComponents(const Vec3 &lightDir, const Vec3 &viewDir, const Vec3 &normal, const Material &material, const Light &light) {
    const Vec3 H = (lightDir + viewDir).normalize();
    const float NdotL = std::max(Vec3::dot(normal, lightDir), 0.0f);
    const float NdotH = std::max(Vec3::dot(normal, H), 0.0f);

    const Vec3 ambient = Vec3::compProduct(material.ambient_material, light.material);
    const Vec3 diffuse = Vec3::compProduct(material.diffuse_material, light.material) * NdotL;
    const Vec3 specular = Vec3::compProduct(material.specular_material, light.material) * static_cast<float>(std::pow(NdotH, material.shininess));

    return ambient + diffuse + specular;
}

Vec3 Lighting::samplePointOnQuad(const Light &light, std::mt19937 &rng) {
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    const float u = dist(rng);
    const float v = dist(rng);

    return light.quad.vertices[0].position * (1 - u) * (1 - v) +
           light.quad.vertices[1].position * u * (1 - v) +
           light.quad.vertices[2].position * u * v +
           light.quad.vertices[3].position * (1 - u) * v;
}

// Compute the Fresnel effect using Schlick's approximation
float Lighting::computeFresnelEffect(const Vec3 &I, const Vec3 &N, const float ior) {
    float cosI = std::clamp(-1.0f, 1.0f, Vec3::dot(I, N));
    float etaI = 1.0f, etaT = ior;
    if (cosI > 0) { std::swap(etaI, etaT); }

    // Compute the sine of the angle using Snell's Law
    const float sinT = etaI / etaT * sqrtf(std::max(0.f, 1 - cosI * cosI));

    // Total internal reflection
    if (sinT >= 1) return 1.0f;

    const float cosT = sqrtf(std::max(0.f, 1 - sinT * sinT));
    cosI = fabsf(cosI);
    const float Rs = ((etaT * cosI) - (etaI * cosT)) / ((etaT * cosI) + (etaI * cosT));
    const float Rp = ((etaI * cosI) - (etaT * cosT)) / ((etaI * cosI) + (etaT * cosT));
    return (Rs * Rs + Rp * Rp) / 2;
}

// Compute the reflected direction using the normal
Vec3 Lighting::computeReflectedDirection(const Vec3 &direction, const Vec3 &normal) {
    return direction - 2.0f * Vec3::dot(direction, normal) * normal;
}

// Compute the refracted direction using Snell's Law
Vec3 Lighting::computeRefractedDirection(const Vec3 &incident, const Vec3 &normal, const float &ior) {
    // Snell's Law: etaI * sin(thetaI) = etaT * sin(thetaT)
    float cosI = std::clamp(Vec3::dot(incident, normal), -1.0f, 1.0f);
    float etaI = 1.0f, etaT = ior;
    Vec3 n = normal;

    // If the ray is entering the material, invert the normal and swap the refractive indices
    if (cosI < 0) { cosI = -cosI; }
    else { std::swap(etaI, etaT); n = -n; }
    const float etaRatio = etaI / etaT;
    const float k = 1 - etaRatio * etaRatio * (1 - cosI * cosI);

    // If k is negative, there is total internal reflection, return a zero vector
    return k < 0 ? Vec3(0, 0, 0) : etaRatio * incident + (etaRatio * cosI - sqrtf(k)) * n;
}