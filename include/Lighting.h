#ifndef LIGHTING_H
#define LIGHTING_H

#include "Vec3.h"
#include "Light.h"
#include <random>

/**
 * Class providing lighting computations for a scene.
 */
class Lighting {
public:
    /**
     * Computes the Phong reflection components for a given light and material.
     * @param lightDir The direction of the light.
     * @param viewDir The direction of the viewer.
     * @param normal The normal at the point of intersection.
     * @param material The material properties of the surface.
     * @param light The light source.
     * @return The computed Phong reflection components as a Vec3.
     */
    static Vec3 computePhongComponents(const Vec3 &lightDir, const Vec3 &viewDir, const Vec3 &normal, const Material &material, const Light &light);

    /**
     * Samples a point on a quad light source.
     * @param light The quad light source.
     * @param rng A random number generator.
     * @return The sampled point on the quad as a Vec3.
     */
    static Vec3 samplePointOnQuad(const Light &light, std::mt19937 &rng);

    /**
     * Computes the Fresnel effect for a given incident vector and normal.
     * @param I The incident vector.
     * @param N The normal vector.
     * @param ior The index of refraction.
     * @return The computed Fresnel effect as a float.
     */
    static float computeFresnelEffect(const Vec3 &I, const Vec3 &N, float ior);

    /**
     * Computes the reflected direction for a given direction and normal.
     * @param direction The incident direction.
     * @param normal The normal vector.
     * @return The computed reflected direction as a Vec3.
     */
    static Vec3 computeReflectedDirection(const Vec3 &direction, const Vec3 &normal);

    /**
     * Computes the refracted direction for a given incident vector, normal, and index of refraction.
     * @param incident The incident vector.
     * @param normal The normal vector.
     * @param ior The index of refraction.
     * @return The computed refracted direction as a Vec3.
     */
    static Vec3 computeRefractedDirection(const Vec3 &incident, const Vec3 &normal, const float &ior);
};

#endif