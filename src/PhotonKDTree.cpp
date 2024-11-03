#include "../include/PhotonKDTree.h"

std::vector<Photon> PhotonKDTree::toVector() const {
    std::vector<Photon> elements;
    getElementsRecursive(root.get(), elements);
    return elements;
}

void PhotonKDTree::getElementsRecursive(const PhotonKDNode *node, std::vector<Photon> &elements) {
    if (!node) return;

    // Add photons from the current node
    elements.push_back(node->data);

    // Recursively gather photons from child nodes
    getElementsRecursive(node->left.get(), elements);
    getElementsRecursive(node->right.get(), elements);
}

PhotonKDNode *PhotonKDTree::buildPhotonBalancedTree(std::vector<Photon> &elements, const int start, const int end, const int depth) {
    if (start >= end) return nullptr;

    int axis = depth % 3; // Determine axis to split on
    const int mid = (start + end) / 2;

    // Partition elements around the median along the selected axis
    std::nth_element(
        elements.begin() + start, elements.begin() + mid, elements.begin() + end,
        [axis](const Photon &a, const Photon &b) {
            return a.position[axis] < b.position[axis];
        }
    );

    // Create a node for the median photon
    auto node = std::make_unique<PhotonKDNode>(elements[mid]);

    // Recursively build left and right subtrees
    node->left = std::unique_ptr<PhotonKDNode>(buildPhotonBalancedTree(elements, start, mid, depth + 1));
    node->right = std::unique_ptr<PhotonKDNode>(buildPhotonBalancedTree(elements, mid + 1, end, depth + 1));

    return node.release();
}

std::vector<Photon> PhotonKDTree::findNearestNeighbors(const Vec3 &point, const float maxDistance) const {
    std::vector<std::pair<float, Photon>> nearestElements;
    findNearestNeighborsRecursive(root.get(), point, maxDistance * maxDistance, 0, nearestElements);

    std::vector<Photon> result;
    for (const auto& pair : nearestElements) {
        result.push_back(pair.second); // Collect only the photon data
    }
    return result;
}

void PhotonKDTree::findNearestNeighborsRecursive(const PhotonKDNode *node, const Vec3 &point, const float maxDistSq, const int depth, std::vector<std::pair<float, Photon> > &nearestElements)
const {
    if (!node) return;

    // Compute squared distance to the current node's photon
    float distSq = (node->data.position - point).squareLength();

    // If within the maximum distance, add the photon to the nearest elements list
    if (distSq < maxDistSq) {
        nearestElements.emplace_back(distSq, node->data);
    }

    // Determine the current splitting axis and the difference along this axis
    const int axis = depth % 3;
    const float diff = point[axis] - node->data.position[axis];

    // Choose which subtree to search first
    const PhotonKDNode* nearNode = diff < 0 ? node->left.get() : node->right.get();
    const PhotonKDNode* farNode = diff < 0 ? node->right.get() : node->left.get();

    // Search in the nearer subtree
    findNearestNeighborsRecursive(nearNode, point, maxDistSq, depth + 1, nearestElements);

    // If the squared distance along the axis is within maxDistSq, search the farther subtree
    if (diff * diff < maxDistSq) {
        findNearestNeighborsRecursive(farNode, point, maxDistSq, depth + 1, nearestElements);
    }
}