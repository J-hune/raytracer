#ifndef LIGHT_H
#define LIGHT_H

#include "Mesh.h"

/**
 * Enum representing the type of light.
 */
enum LightType {
    LightType_Spherical, ///< Spherical light type.
    LightType_Quad ///< Quad light type.
};

/**
 * Structure representing a light source.
 */
struct Light {
    Vec3 position;          ///< Position of the light.
    LightType type;         ///< Type of the light.
    Vec3 material;          ///< Material properties of the light.
    float radius{};         ///< Radius of the light (for spherical lights).
    float powerCorrection;  ///< Power correction factor for the light.
    bool isInCamSpace{};    ///< Flag indicating if the light is in camera space.
    mutable Mesh quad;      ///< Mesh representing the quad light.

    Light() : type(), powerCorrection(1.0) {}
    Light(
        const Vec3 &pos, const LightType type, const Vec3 &material,
        const float radius, const float powerCorrection, const bool isInCamSpace)
        : position(pos), type(type), material(material), radius(radius),
          powerCorrection(powerCorrection), isInCamSpace(isInCamSpace) {
    }
};

#endif //LIGHT_H
