#ifndef LIGHT_H
#define LIGHT_H

#include "Mesh.h"

enum LightType {
    LightType_Spherical,
    LightType_Quad
};

struct Light {
    Vec3 position;
    LightType type;
    Vec3 material;
    float radius{};
    float powerCorrection;
    bool isInCamSpace{};
    mutable Mesh quad;

    Light() : type(), powerCorrection(1.0) {
    }

    Light(
        const Vec3 &pos, const LightType type, const Vec3 &material,
        const float radius, const float powerCorrection, const bool isInCamSpace)
        : position(pos), type(type), material(material), radius(radius),
          powerCorrection(powerCorrection), isInCamSpace(isInCamSpace) {
    }
};


#endif //LIGHT_H
