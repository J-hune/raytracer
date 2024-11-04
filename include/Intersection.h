#ifndef INTERSECTION_H
#define INTERSECTION_H

#include <vector>
#include "Ray.h"
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"
#include "Triangle.h"
#include <cfloat>

#include "MeshKDTree.h"

/**
 * Structure representing the intersection of a ray with a scene.
 */
struct RaySceneIntersection {
    bool intersectionExists;                        ///< Indicates if an intersection exists.
    unsigned int typeOfIntersectedObject{};         ///< Type of the intersected object.
    unsigned int objectIndex{};                     ///< Index of the intersected object.
    float t;                                        ///< Distance to the intersection point.
    RayTriangleIntersection rayMeshIntersection;    ///< Intersection details with a mesh.
    RaySphereIntersection raySphereIntersection;    ///< Intersection details with a sphere.
    RaySquareIntersection raySquareIntersection;    ///< Intersection details with a square.

    RaySceneIntersection() : intersectionExists(false), t(FLT_MAX) {}
};

/**
 * Class providing methods to compute and parse intersections in a scene.
 */
class Intersection {
public:
    /**
     * Computes the intersection of a ray with the objects in the scene.
     * @param ray The ray to test for intersections.
     * @param spheres The list of spheres in the scene.
     * @param squares The list of squares in the scene.
     * @param meshes The list of meshes in the scene.
     * @param kd_tree The KD-tree for mesh acceleration.
     * @param z_near The near clipping plane distance.
     * @return A RaySceneIntersection object containing the intersection details.
     */
    static RaySceneIntersection computeIntersection(
        const Ray &ray, const std::vector<Sphere> &spheres, const std::vector<Square> &squares, const std::vector<Mesh> &meshes,
        const MeshKDTree &kd_tree, float z_near
    );

    /**
     * Parses the intersection details to extract position, normal, and material.
     * @param intersection The RaySceneIntersection object containing the intersection details.
     * @param spheres The list of spheres in the scene.
     * @param squares The list of squares in the scene.
     * @param meshes The list of meshes in the scene.
     * @return A tuple containing the position, normal, and material of the intersection.
     */
    static std::tuple<Vec3, Vec3, Material> parseIntersection(
        const RaySceneIntersection &intersection,
        const std::vector<Sphere> &spheres,
        const std::vector<Square> &squares,
        const std::vector<Mesh> &meshes
    );
};

#endif // INTERSECTION_H