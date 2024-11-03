#ifndef MESHKD_TREE_H
#define MESHKD_TREE_H

#include <vector>
#include <memory>

#include "Mesh.h"

// Represents a node in the KD tree
struct MeshKDNode {
    AABB aabb;
    std::vector<Triangle> triangles;
    std::unique_ptr<MeshKDNode> left;
    std::unique_ptr<MeshKDNode> right;

    explicit MeshKDNode(const AABB &aabb, const std::vector<Triangle> &triangles) : aabb(aabb), triangles(triangles) {}
};

// Represents the KD tree
class MeshKDTree {
public:
    MeshKDTree() = default;
    explicit MeshKDTree(const std::vector<Mesh>& elements) {
        root = std::unique_ptr<MeshKDNode>(buildBalancedTree(elements));
    }

    [[nodiscard]] bool isEmpty() const { return !root; }
    [[nodiscard]] std::vector<Triangle> toVector() const;
    [[nodiscard]] bool intersect(const Ray &ray, RayTriangleIntersection &intersection) const;
    void draw();


private:
    std::vector<AABB> AABBtoDraw;
    std::unique_ptr<MeshKDNode> root;
    void getElementsRecursive(const MeshKDNode *node, std::vector<Triangle> &elements) const;
    MeshKDNode *buildBalancedTree(const std::vector<Mesh> &elements);
    std::unique_ptr<MeshKDNode> buildRecursiveBalancedTree(const std::vector<Triangle> &vector, const AABB &aabb, int depth, int maxDepth);
    bool traverseTree(const MeshKDNode *node, const Ray &ray, RayTriangleIntersection &intersection) const;
    std::vector<AABB> drawRecursive(const MeshKDNode *node) const;
};

#endif // MESHKD_TREE_H