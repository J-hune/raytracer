#ifndef INTERSECTION_H
#define INTERSECTION_H

#include <vector>
#include "Mesh.h"
#include "Ray.h"
#include "Sphere.h"
#include "Square.h"

#include "MeshKDTree.h"

/**
 * Structure representing the intersection of a ray with a scene.
 */
struct RaySceneIntersection : RayIntersection{
    Vec3 textureColor;          ///< Texture color at the intersection point.
    Material material;          ///< Material at the intersection point.
    RaySceneIntersection() : RayIntersection() {}
    explicit RaySceneIntersection(const RayIntersection &intersection) : RayIntersection(intersection) {}
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
     * Tests whether anything blocks the segment between the ray origin and maxDistance.
     *
     * A shadow ray only needs to know whether the light is occluded, never by what. Compared to
     * computeIntersection this skips the texture lookup and the material copy performed on every
     * candidate hit, and returns as soon as one blocker is found instead of looking for the nearest one.
     *
     * @param ray The shadow ray to test.
     * @param spheres The list of spheres in the scene.
     * @param squares The list of squares in the scene.
     * @param kd_tree The KD-tree for mesh acceleration.
     * @param maxDistance Distance to the light sample; a hit beyond it does not block.
     * @param epsilon Near cutoff, guards against the surface shadowing itself.
     * @return True if something blocks the light.
     */
    static bool isOccluded(
        const Ray &ray, const std::vector<Sphere> &spheres, const std::vector<Square> &squares,
        const MeshKDTree &kd_tree, float maxDistance, float epsilon
    );
};

#endif // INTERSECTION_H