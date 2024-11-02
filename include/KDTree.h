#ifndef KDTREE_H
#define KDTREE_H

#include "Vec3.h"
#include <vector>
#include <memory>

// Represents a node in the KD tree
template<typename T>
struct KDNode {
    T data;
    std::unique_ptr<KDNode> left;
    std::unique_ptr<KDNode> right;

    explicit KDNode(T data) : data(data) {}
};

// Represents the KD tree
template<typename T>
class KDTree {
public:
    KDTree() = default;
    explicit KDTree(std::vector<T>& elements);
    [[nodiscard]] std::vector<T> toVector() const;
    std::vector<T> findNearestNeighbors(const Vec3 &point, float maxDistance) const;

private:
    std::unique_ptr<KDNode<T>> root;
    KDNode<T> *buildBalancedTree(std::vector<T> &elements, int start, int end, int depth);
    void getElementsRecursive(const KDNode<T> *node, std::vector<T> &elements) const;
    void findNearestNeighborsRecursive(const KDNode<T> *node, const Vec3 &point, float maxDistSq, int depth, std::vector<std::pair<float, T>> &nearestElements) const;
};


/*********************************************************************************************************************/
/************************ Implementations (in cpp, generics must be defined in the same file) ************************/
/*********************************************************************************************************************/

template<typename T>
KDTree<T>::KDTree(std::vector<T> &elements) {
    root = std::unique_ptr<KDNode<T> >(buildBalancedTree(elements, 0, elements.size(), 0));
}

template<typename T>
KDNode<T> *KDTree<T>::buildBalancedTree(std::vector<T> &elements, int start, int end, int depth) {
    if (start >= end) return nullptr;

    int axis = depth % 3;
    const int mid = (start + end) / 2;
    std::nth_element(
        elements.begin() + start, elements.begin() + mid, elements.begin() + end,
        [axis](const T &a, const T &b) {
            return a.position[axis] < b.position[axis];
        }
    );

    auto node = std::make_unique<KDNode<T>>(elements[mid]);
    node->left = std::unique_ptr<KDNode<T>>(buildBalancedTree(elements, start, mid, depth + 1));
    node->right = std::unique_ptr<KDNode<T>>(buildBalancedTree(elements, mid + 1, end, depth + 1));

    return node.release();
}

template<typename T>
std::vector<T> KDTree<T>::toVector() const {
    std::vector<T> elements;
    getElementsRecursive(root.get(), elements);
    return elements;
}

template<typename T>
void KDTree<T>::getElementsRecursive(const KDNode<T>* node, std::vector<T>& elements) const {
    if (!node) return;
    elements.push_back(node->data);
    getElementsRecursive(node->left.get(), elements);
    getElementsRecursive(node->right.get(), elements);
}

template<typename T>
std::vector<T> KDTree<T>::findNearestNeighbors(const Vec3 &point, const float maxDistance) const {
    std::vector<std::pair<float, T>> nearestElements;
    findNearestNeighborsRecursive(root.get(), point, maxDistance * maxDistance, 0, nearestElements);

    std::vector<T> result;
    for (const auto& pair : nearestElements) {
        result.push_back(pair.second);
    }
    return result;
}

template<typename T>
void KDTree<T>::findNearestNeighborsRecursive(const KDNode<T> *node, const Vec3 &point, const float maxDistSq, const int depth,
                                              std::vector<std::pair<float, T> > &nearestElements)
const {
    if (!node) return;

    // Calculate the squared distance
    float distSq = (node->data.position - point).squareLength();

    // Add element if within the distance
    if (distSq < maxDistSq) {
        nearestElements.emplace_back(distSq, node->data);
    }

    // Determine the axis and the difference
    const int axis = depth % 3;
    const float diff = point[axis] - node->data.position[axis];

    // Determine which subtree to explore first
    const KDNode<T>* nearNode = diff < 0 ? node->left.get() : node->right.get();
    const KDNode<T>* farNode = diff < 0 ? node->right.get() : node->left.get();

    // Search in the nearer subtree
    findNearestNeighborsRecursive(nearNode, point, maxDistSq, depth + 1, nearestElements);

    // If the distance to the axis is less than the maximum distance, explore the other subtree
    if (diff * diff < maxDistSq) {
        findNearestNeighborsRecursive(farNode, point, maxDistSq, depth + 1, nearestElements);
    }
}

#endif // KDTREE_H