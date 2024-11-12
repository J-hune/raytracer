#include "../include/MeshKDTree.h"
#include "../include/Triangle.h"
#include "../include/Settings.h"

std::vector<Triangle> MeshKDTree::toVector() const {
    std::vector<Triangle> elements;
    getElementsRecursive(root.get(), elements);
    return elements;
}

void MeshKDTree::getElementsRecursive(const MeshKDNode *node, std::vector<Triangle> &elements) const {
    if (!node) return;

    // Add triangles from the current node
    elements.insert(elements.end(), node->triangles.begin(), node->triangles.end());

    // Recursively gather triangles from child nodes
    getElementsRecursive(node->left.get(), elements);
    getElementsRecursive(node->right.get(), elements);
}

std::unique_ptr<MeshKDNode> MeshKDTree::buildBalancedTree(const std::vector<Mesh> &elements) {
    const Settings &settings = Settings::getInstance();
    if (elements.empty()) return nullptr;

    // Collect all triangles from the meshes
    std::vector<Triangle> sceneTriangles;
    for (const auto &mesh : elements) {
        std::vector<Triangle> meshTriangles = mesh.getTriangles();
        sceneTriangles.insert(sceneTriangles.end(), meshTriangles.begin(), meshTriangles.end());
    }

    // Calculate the axis-aligned bounding box (AABB) of the scene
    const AABB sceneAABB = Triangle::getAABB(sceneTriangles);
    auto node = buildRecursiveBalancedTree(sceneTriangles, sceneAABB, 0, settings.maxKdTreeDepth);
    return node;
}

std::unique_ptr<MeshKDNode> MeshKDTree::buildRecursiveBalancedTree(const std::vector<Triangle> &vector, const AABB &aabb, const int depth, int maxDepth) {
    // Stop recursion if maximum depth is reached or no triangles left
    if (depth >= maxDepth || vector.empty()) {
        return nullptr;
    }

    const int axis = depth % 3; // Determine the axis to split on
    std::vector<Triangle> leftTriangles, rightTriangles;

    for (const auto &triangle : vector) {
        const Vec3 &v0 = triangle.getVertex(0);
        const Vec3 &v1 = triangle.getVertex(1);
        const Vec3 &v2 = triangle.getVertex(2);
        const Vec3 min = Vec3::min(v0, Vec3::min(v1, v2));

        if (min[axis] < aabb.center()[axis]) {
            leftTriangles.push_back(triangle);
        } else {
            rightTriangles.push_back(triangle);
        }
    }

    // Calculate AABBs for left and right triangles
    AABB leftAABB, rightAABB;
    for (const auto &triangle : leftTriangles) {
        leftAABB = AABB::merge(leftAABB, triangle.getAABB());
    }
    for (const auto &triangle : rightTriangles) {
        rightAABB = AABB::merge(rightAABB, triangle.getAABB());
    }
    auto node = std::make_unique<MeshKDNode>(AABB::merge(leftAABB, rightAABB), vector);

    // Create child nodes
    node->left = buildRecursiveBalancedTree(leftTriangles, leftAABB, depth + 1, maxDepth);
    node->right = buildRecursiveBalancedTree(rightTriangles, rightAABB, depth + 1, maxDepth);

    return node;
}

bool MeshKDTree::intersect(const Ray &ray, RayTriangleIntersection &intersection) const {
    intersection.t = std::numeric_limits<float>::max();
    return traverseTree(root.get(), ray, intersection);
}

bool MeshKDTree::traverseTree(const MeshKDNode *node, const Ray &ray, RayTriangleIntersection &intersection) const {
    if (!node) return false;

    // Check if the ray intersects the AABB of this node
    if (!ray.intersectAABB(node->aabb)) {
        return false;
    }

    // If it's a leaf node, check for intersections with its triangles
    if (node->left == nullptr && node->right == nullptr) {
        bool intersectionFound = false;
        for (const auto &triangle : node->triangles) {
            const RayTriangleIntersection currentIntersection = triangle.intersect(ray);
            if (currentIntersection.intersectionExists && currentIntersection.t < intersection.t) {
                intersection = currentIntersection;
                intersectionFound = true;
            }
        }
        return intersectionFound;
    }

    // Check both child nodes for intersections
    const bool leftIntersection = traverseTree(node->left.get(), ray, intersection);
    const bool rightIntersection = traverseTree(node->right.get(), ray, intersection);

    return leftIntersection || rightIntersection;
}

void MeshKDTree::draw() {
    // If the AABB vector is empty, fill it with the AABBs to draw (it's expensive)
    if (AABBtoDraw.empty()) AABBtoDraw = drawRecursive(root.get());

    // Draw each AABB
    for (const auto &aabb : AABBtoDraw) {
        aabb.draw();
    }
}

std::vector<AABB> MeshKDTree::drawRecursive(const MeshKDNode *node) const {
    if (!node) return {};

    std::vector<AABB> aabbs;
    aabbs.push_back(node->aabb);

    // Recursively gather AABBs from left and right child nodes
    if (node->left) {
        std::vector<AABB> leftAABBs = drawRecursive(node->left.get());
        aabbs.insert(aabbs.end(), leftAABBs.begin(), leftAABBs.end());
    }

    if (node->right) {
        std::vector<AABB> rightAABBs = drawRecursive(node->right.get());
        aabbs.insert(aabbs.end(), rightAABBs.begin(), rightAABBs.end());
    }

    return aabbs;
}