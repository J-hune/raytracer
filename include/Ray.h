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

    [[nodiscard]] bool intersectAABB(const AABB &aabb) const {
        Vec3 invDir = Vec3(1.f) / direction();
        Vec3 orig = origin();

        float tMin = -FLT_MAX;
        float tMax = FLT_MAX;

        for (int i = 0; i < 3; ++i) {
            float t1 = (aabb.getMin()[i] - orig[i]) * invDir[i];
            float t2 = (aabb.getMax()[i] - orig[i]) * invDir[i];

            tMin = std::max(tMin, std::min(t1, t2));
            tMax = std::min(tMax, std::max(t1, t2));
        }

        return tMax >= tMin && tMax >= 0;
    }
};
#endif
