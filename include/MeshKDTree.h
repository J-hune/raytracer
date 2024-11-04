#ifndef MESHKD_TREE_H
#define MESHKD_TREE_H

#include <vector>
#include <memory>

#include "Mesh.h"

// -------------------------------------------
// MeshKDNode Structure
// -------------------------------------------

/**
 * Structure representing a node in the KD tree.
 */
struct MeshKDNode {
    AABB aabb;                          ///< Axis-aligned bounding box of the node.
    std::vector<Triangle> triangles;    ///< Triangles contained in the node.
    std::unique_ptr<MeshKDNode> left;   ///< Pointer to the left child node.
    std::unique_ptr<MeshKDNode> right;  ///< Pointer to the right child node.

    /**
     * Constructor for MeshKDNode.
     * @param aabb Axis-aligned bounding box of the node.
     * @param triangles Triangles contained in the node.
     */
    explicit MeshKDNode(const AABB &aabb, const std::vector<Triangle> &triangles) : aabb(aabb), triangles(triangles) {}
};

// -------------------------------------------
// MeshKDTree Class
// -------------------------------------------

/**
 * Class representing the KD tree.
 */
class MeshKDTree {
public:
    MeshKDTree() = default;
    explicit MeshKDTree(const std::vector<Mesh>& elements) {
        root = std::unique_ptr<MeshKDNode>(buildBalancedTree(elements));
    }

    /**
     * Checks if the KD tree is empty.
     * @return True if the tree is empty, false otherwise.
     */
    [[nodiscard]] bool isEmpty() const { return !root; }

    /**
     * Converts the KD tree to a vector of triangles.
     * @return Vector of triangles contained in the tree.
     */
    [[nodiscard]] std::vector<Triangle> toVector() const;

    /**
     * Checks if a ray intersects with any triangle in the KD tree.
     * @param ray The ray to check for intersection.
     * @param intersection The intersection information.
     * @return True if the ray intersects with any triangle, false otherwise.
     */
    [[nodiscard]] bool intersect(const Ray &ray, RayTriangleIntersection &intersection) const;

    /**
     * Draws the KD tree.
     */
    void draw();

private:
    std::vector<AABB> AABBtoDraw;       ///< Vector of axis-aligned bounding boxes to draw.
    std::unique_ptr<MeshKDNode> root;   ///< Root node of the KD tree.

    /**
     * Recursively retrieves elements from the KD tree.
     * @param node The current node.
     * @param elements Vector to store the retrieved elements.
     */
    void getElementsRecursive(const MeshKDNode *node, std::vector<Triangle> &elements) const;

    /**
     * Builds a balanced KD tree from a vector of meshes.
     * @param elements Vector of meshes to build the tree from.
     * @return Pointer to the root node of the built KD tree.
     */
    MeshKDNode *buildBalancedTree(const std::vector<Mesh> &elements);

    /**
     * Recursively builds a balanced KD tree.
     * @param vector Vector of triangles to build the tree from.
     * @param aabb Axis-aligned bounding box of the current node.
     * @param depth Current depth of the tree.
     * @param maxDepth Maximum depth of the tree.
     * @return Unique pointer to the built node.
     */
    std::unique_ptr<MeshKDNode> buildRecursiveBalancedTree(const std::vector<Triangle> &vector, const AABB &aabb, int depth, int maxDepth);

    /**
     * Recursively traverses the KD tree to check for ray intersections.
     * @param node The current node.
     * @param ray The ray to check for intersection.
     * @param intersection The intersection information.
     * @return True if the ray intersects with any triangle, false otherwise.
     */
    bool traverseTree(const MeshKDNode *node, const Ray &ray, RayTriangleIntersection &intersection) const;

    /**
     * Recursively retrieves axis-aligned bounding boxes to draw.
     * @param node The current node.
     * @return Vector of axis-aligned bounding boxes.
     */
    std::vector<AABB> drawRecursive(const MeshKDNode *node) const;
};

#endif // MESHKD_TREE_H