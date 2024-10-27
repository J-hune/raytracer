#ifndef LIGHT_H
#define LIGHT_H

#include "Mesh.h"

enum LightType {
    LightType_Spherical,
    LightType_Quad
};

struct Light {
    Vec3 material;
    bool isInCamSpace{};
    LightType type;

    Vec3 pos;
    float radius{};

    Mesh quad;

    float powerCorrection;

    Light() : type(), powerCorrection(1.0) {
    }
};


#endif //LIGHT_H