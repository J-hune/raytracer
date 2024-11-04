#ifndef RAY_H
#define RAY_H

#include "Line.h"
#include "AABB.h"
#include <cfloat>

// -------------------------------------------
// RayIntersection Structure
// -------------------------------------------

/**
 * Structure representing the result of a ray intersection.
 */
struct RayIntersection {
    bool intersectionExists{};   ///< Indicates if an intersection exists.
    float t{};                   ///< The distance from the ray origin to the intersection point.
    Vec3 intersection;           ///< The intersection point.
    Vec3 normal;                 ///< The normal at the intersection point.
};

// -------------------------------------------
// Ray Class
// -------------------------------------------

/**
 * Class representing a ray, inheriting from Line.
 */
class Ray : public Line {
public:
    float tMax;     ///< The maximum distance the ray can travel.

    Ray() : Line(), tMax(FLT_MAX) {}

    /**
     * Constructor for Ray with origin, direction, and maximum distance.
     * @param o The origin of the ray.
     * @param d The direction of the ray.
     * @param tMax The maximum distance the ray can travel.
     */
    Ray(const Vec3 &o, const Vec3 &d, const float tMax = FLT_MAX) : Line(o, d), tMax(tMax) {}

    /**
     * Checks if the ray intersects with an axis-aligned bounding box (AABB).
     * @param aabb The axis-aligned bounding box to check for intersection.
     * @return True if the ray intersects with the AABB, false otherwise.
     */
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

#endif // RAY_H