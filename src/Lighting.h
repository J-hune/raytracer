#ifndef LIGHTING_H
#define LIGHTING_H

#include <random>
#include <bits/stl_algo.h>

#include "Vec3.h"
#include "Material.h"
#include "Light.h"

class Lighting {
public:
    static Vec3 computePhongComponents(const Vec3 &lightDir, const Vec3 &viewDir, const Vec3 &normal, const Material &material, const Light &light);
    static Vec3 samplePointOnQuad(const Light &light, std::mt19937 &rng);
    static float computeFresnelEffect(const Vec3 &I, const Vec3 &N, float ior);
    static Vec3 computeReflectedDirection(const Vec3 &direction, const Vec3 &normal);
    static Vec3 computeRefractedDirection(const Vec3 &incident, const Vec3 &normal, const float &ior);
};

#endif