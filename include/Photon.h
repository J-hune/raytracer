#ifndef PHOTON_H
#define PHOTON_H

#include "Material.h"
#include "Vec3.h"

struct Photon {
    Vec3 position;
    Vec3 direction;
    Vec3 color;
    MaterialType materialType;

    [[nodiscard]] Vec3 getPosition() const { return position; }
};

#endif //PHOTON_H
