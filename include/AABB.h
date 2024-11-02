#ifndef AABB_H
#define AABB_H

#include "Vec3.h"

struct AABB {
    Vec3 min;
    Vec3 max;

    AABB() : min(0.0f), max(0.0f) {}
    AABB(const Vec3& _min, const Vec3& _max) : min(_min), max(_max) {}
    AABB(const AABB& aabb) = default;
};

#endif //AABB_H
