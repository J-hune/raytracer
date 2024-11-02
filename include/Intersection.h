#ifndef INTERSECTION_H
#define INTERSECTION_H

#include <vector>
#include "Ray.h"
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"
#include "Triangle.h"
#include <cfloat>

struct RaySceneIntersection {
    bool intersectionExists;
    unsigned int typeOfIntersectedObject{};
    unsigned int objectIndex{};
    float t;
    RayTriangleIntersection rayMeshIntersection;
    RaySphereIntersection raySphereIntersection;
    RaySquareIntersection raySquareIntersection;

    RaySceneIntersection() : intersectionExists(false), t(FLT_MAX) {}
};

class Intersection {
public:
    static RaySceneIntersection computeIntersection(const Ray &ray, const std::vector<Sphere> &spheres, const std::vector<Square> &squares, const std::vector<Mesh> &meshes, float z_near);
    static std::tuple<Vec3, Vec3, Material> parseIntersection(const RaySceneIntersection &intersection, const std::vector<Sphere> &spheres, const std::vector<Square> &squares, const std::vector<Mesh> &meshes);
};

#endif // INTERSECTION_H