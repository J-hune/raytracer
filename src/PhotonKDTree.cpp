#include "../include/PhotonKDTree.h"

#include <cmath>
#include <queue>

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

std::vector<Photon> PhotonKDTree::findNearestNeighbors(const Vec3 &point, const Vec3 &normal, const float maxDistance, const int maxCount) const {
    // Priority queue to store the nearest photons, with the farthest on top
    std::priority_queue<std::pair<float, Photon>, std::vector<std::pair<float, Photon>>, PhotonDistanceComparator> nearestElements;
    findNearestNeighborsRecursive(root.get(), point, normal, maxDistance * maxDistance, 0, maxCount, nearestElements);

    // Extract photons from the priority queue
    std::vector<Photon> result;
    while (!nearestElements.empty()) {
        result.push_back(nearestElements.top().second); // Collect only the photon data
        nearestElements.pop();
    }

    // Warning: the result vector is in reverse order (farthest to nearest)
    return result;
}

void PhotonKDTree::findNearestNeighborsRecursive(
    const PhotonKDNode *node, const Vec3 &point, const Vec3 &normal, const float maxDistSq, const int depth, const int maxCount,
    std::priority_queue<std::pair<float, Photon>, std::vector<std::pair<float, Photon>>, PhotonDistanceComparator> &nearestElements) {

    if (!node) return;

    // Calculate the vector from the point to the photon
    const Vec3 toPhoton = node->data.position - point;

    // Calculate the squared distance to the photon center
    float distSqToCenter = toPhoton.squareLength();

    // Calculate the distance to the plane defined by the normal
    const float distToPlane = Vec3::dot(toPhoton, normal);

    // Check if the photon is within the disk region (both distance to center and distance to the plane)
    if (distSqToCenter < maxDistSq && distToPlane < 1e-5f) {
        // Add photon to the results if it's within the allowed distance
        if (static_cast<int>(nearestElements.size()) < maxCount) {
            nearestElements.emplace(distSqToCenter, node->data);
        } else if (!nearestElements.empty() && distSqToCenter < nearestElements.top().first) {
            // Replace the farthest photon if this one is closer
            nearestElements.pop();
            nearestElements.emplace(distSqToCenter, node->data);
        }
    }

    // Determine the current splitting axis and the difference along this axis
    const int axis = depth % 3;
    const float diff = point[axis] - node->data.position[axis];

    // Choose which subtree to search first
    const PhotonKDNode* nearNode = diff < 0 ? node->left.get() : node->right.get();
    const PhotonKDNode* farNode = diff < 0 ? node->right.get() : node->left.get();

    // Search in the nearer subtree
    findNearestNeighborsRecursive(nearNode, point, normal, maxDistSq, depth + 1, maxCount, nearestElements);

    // If the squared distance along the axis is within maxDistSq, search the farther subtree
    if (diff * diff < maxDistSq && (static_cast<int>(nearestElements.size()) < maxCount || diff * diff < nearestElements.top().first)) {
        findNearestNeighborsRecursive(farNode, point, normal, maxDistSq, depth + 1, maxCount, nearestElements);
    }
}