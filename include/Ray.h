#ifndef RAY_H
#define RAY_H

#include "Line.h"
#include "AABB.h"
#include <cfloat>

struct RayIntersection {
    bool intersectionExists{};
    float t{};
    Vec3 intersection;
    Vec3 normal;
};

class Ray : public Line {
public:
    float tMax;

    Ray() : Line(), tMax(FLT_MAX) {}
    Ray(const Vec3 &o, const Vec3 &d, const float tMax = FLT_MAX) : Line(o, d), tMax(tMax) {}

    [[nodiscard]] bool intersectAABB(AABB aabb) const {
        Vec3 invDirection = Vec3(1.f) / direction();
        const float t1 = (aabb.min[0] - origin()[0]) * invDirection[0];
        const float t2 = (aabb.max[0] - origin()[0]) * invDirection[0];
        const float t3 = (aabb.min[1] - origin()[1]) * invDirection[1];
        const float t4 = (aabb.max[1] - origin()[1]) * invDirection[1];
        const float t5 = (aabb.min[2] - origin()[2]) * invDirection[2];
        const float t6 = (aabb.max[2] - origin()[2]) * invDirection[2];

        const float tMin = std::max(std::max(std::min(t1, t2), std::min(t3, t4)), std::min(t5, t6));
        const float tMax = std::min(std::min(std::max(t1, t2), std::max(t3, t4)), std::max(t5, t6));

        return tMax >= tMin && tMax >= 0 && tMin <= tMax;
    }
};
#endif
