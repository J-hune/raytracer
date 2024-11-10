#include "../include/Intersection.h"
#include <cfloat>

#include "../include/MeshKDTree.h"
#include "../include/Settings.h"

RaySceneIntersection Intersection::computeIntersection(const Ray &ray, const std::vector<Sphere> &spheres,
const std::vector<Square> &squares, const std::vector<Mesh> &meshes, const MeshKDTree &kd_tree, const float z_near) {
    const Settings &settings = Settings::getInstance();
    RaySceneIntersection result;
    result.t = FLT_MAX;

    // For each object in the scene, compute the intersection with the ray
    for (const auto &sphere : spheres) {
        // We first check if the ray intersects an AABB around the sphere
        if (!sphere.intersectAABB(ray)) {
            continue;
        }

        const RayIntersection intersection = sphere.intersect(ray);

        // If the intersection exists and is closer than the previous one, keep it
        if (intersection.intersectionExists && intersection.t <= result.t) {
            result = RaySceneIntersection(intersection);
            result.material = sphere.material;
            result.textureColor = sphere.sampleTexture(intersection.u, intersection.v);
        }
    }

    for (const auto &square : squares) {
        // We first check if the ray intersects an AABB around the square
        if (!square.intersectAABB(ray)) {
            continue;
        }

        const RayIntersection intersection = square.intersect(ray);

        // If the intersection exists and is closer than the previous one, keep it
        if (intersection.intersectionExists && intersection.t <= result.t && intersection.t > z_near) {
            result = RaySceneIntersection(intersection);
            result.material = square.material;
            result.textureColor = square.sampleTexture(intersection.u, intersection.v);
        }
    }

    // Compute the intersection with the KD tree
    if (settings.useKDTree) {
        RayTriangleIntersection intersection;
        if (!kd_tree.isEmpty() && kd_tree.intersect(ray, intersection)) {
            // If the intersection exists and is closer than the previous one, keep it
            if (intersection.intersectionExists && intersection.t <= result.t && intersection.t > z_near) {
                result = RaySceneIntersection(intersection);
                result.material = meshes[intersection.tIndex].material;
                result.textureColor = meshes[intersection.tIndex].sampleTexture(intersection.u, intersection.v);
            }
        }
    }

    // I'm keeping the old mesh intersection implementation (without KD-Tree) for rendering time comparisons.
    // And also because I'm a little sentimental about it. :')
    else {
        for (const auto &mesh : meshes) {
            // We first check if the ray intersects an AABB around the mesh
            if (!mesh.intersectAABB(ray)) {
                continue;
            }

            const RayIntersection intersection = mesh.intersect(ray);

            // If the intersection exists and is closer than the previous one, keep it
            if (intersection.intersectionExists && intersection.t <= result.t && intersection.t > z_near) {
                result = RaySceneIntersection(intersection);
                result.material = mesh.material;
            }
        }
    }

    return result;
}