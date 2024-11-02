#ifndef PHOTON_H
#define PHOTON_H

#include "Material.h"
#include "Vec3.h"

struct Photon {
    Vec3 position;
    Vec3 direction;
    Vec3 color;
    MaterialType materialType;
};

#endif //PHOTON_H
