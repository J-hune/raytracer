#include "../include/Intersection.h"
#include <cfloat>

RaySceneIntersection Intersection::computeIntersection(const Ray &ray, const std::vector<Sphere> &spheres, const std::vector<Square> &squares, const std::vector<Mesh> &meshes, float z_near) {
    RaySceneIntersection result;
    result.t = FLT_MAX;

    // For each object in the scene, compute the intersection with the ray
    for (unsigned int i = 0; i < spheres.size(); i++) {
        // We first check if the ray intersects an AABB around the sphere
        if (!spheres[i].intersectAABB(ray)) {
            continue;
        }

        const RayIntersection intersection = spheres[i].intersect(ray);

        // If the intersection exists and is closer than the previous one, keep it
        if (intersection.intersectionExists && intersection.t <= result.t) {
            result.intersectionExists = true;
            result.raySphereIntersection = static_cast<RaySphereIntersection>(intersection);
            result.t = intersection.t;
            result.objectIndex = i;
            result.typeOfIntersectedObject = 0;
        }
    }

    for (unsigned int i = 0; i < squares.size(); i++) {
        // We first check if the ray intersects an AABB around the square
        if (!squares[i].intersectAABB(ray)) {
            continue;
        }

        const RayIntersection intersection = squares[i].intersect(ray);

        // If the intersection exists and is closer than the previous one, keep it
        if (intersection.intersectionExists && intersection.t <= result.t && intersection.t > z_near) {
            result.intersectionExists = true;
            result.raySquareIntersection = static_cast<RaySquareIntersection>(intersection);
            result.t = intersection.t;
            result.objectIndex = i;
            result.typeOfIntersectedObject = 1;
        }
    }

    for (unsigned int i = 0; i < meshes.size(); i++) {
        // We first check if the ray intersects an AABB around the mesh
        if (!meshes[i].intersectAABB(ray)) {
            continue;
        }

        const RayIntersection intersection = meshes[i].intersect(ray);

        // If the intersection exists and is closer than the previous one, keep it
        if (intersection.intersectionExists && intersection.t <= result.t && intersection.t > z_near) {
            result.intersectionExists = true;
            result.rayMeshIntersection = static_cast<RayTriangleIntersection>(intersection);
            result.t = intersection.t;
            result.objectIndex = i;
            result.typeOfIntersectedObject = 2;
        }
    }

    return result;
}

std::tuple<Vec3, Vec3, Material> Intersection::parseIntersection(const RaySceneIntersection &intersection,
    const std::vector<Sphere> &spheres, const std::vector<Square> &squares, const std::vector<Mesh> &meshes) {
    Vec3 intersectionPoint, normal;
    Material material;

    // Determine the intersected object (sphere, square, mesh)
    if (intersection.typeOfIntersectedObject == 0) {
        // Sphere
        intersectionPoint = intersection.raySphereIntersection.intersection;
        normal = intersection.raySphereIntersection.normal;
        material = spheres[intersection.objectIndex].material;
    } else if (intersection.typeOfIntersectedObject == 1) {
        // Square
        intersectionPoint = intersection.raySquareIntersection.intersection;
        normal = intersection.raySquareIntersection.normal;
        material = squares[intersection.objectIndex].material;
    } else {
        // Mesh
        intersectionPoint = intersection.rayMeshIntersection.intersection;
        normal = intersection.rayMeshIntersection.normal;
        material = meshes[intersection.objectIndex].material;
    }

    return {intersectionPoint, normal, material};
}
