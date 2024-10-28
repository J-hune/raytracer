#include "PhotonKDTree.h"
#include <algorithm>

KDNode::KDNode(Photon photon) : photon(std::move(photon)), left(nullptr), right(nullptr) {}

PhotonKDTree::PhotonKDTree(std::vector<Photon>& photons) {
    root = std::unique_ptr<KDNode>(buildBalancedTree(photons, 0, static_cast<int>(photons.size()), 0));
}

std::vector<Photon> PhotonKDTree::toVector() const {
    std::vector<Photon> photons;
    getPhotonsRecursive(root.get(), photons);
    return photons;
}

std::vector<Photon> PhotonKDTree::findNearestNeighbors(const Vec3& point, const float maxDistance) const {
    std::vector<std::pair<float, Photon>> nearestPhotons;
    findNearestNeighborsRecursive(root.get(), point, maxDistance * maxDistance, 0, nearestPhotons);
    std::vector<Photon> result;
    for (const auto&[_, snd] : nearestPhotons) {
        result.push_back(snd);
    }
    return result;
}

KDNode* PhotonKDTree::buildBalancedTree(std::vector<Photon>& photons, const int start, const int end, const int depth) {
    if (start >= end) return nullptr;
    int axis = depth % 3; // Choose the axis based on the depth
    const int mid = (start + end) / 2;
    std::nth_element(
        photons.begin() + start, photons.begin() + mid, photons.begin() + end,
        [axis](const Photon &a, const Photon &b) {
            return a.position[axis] < b.position[axis];
        });

    auto node = std::make_unique<KDNode>(photons[mid]);
    node->left = std::unique_ptr<KDNode>(buildBalancedTree(photons, start, mid, depth + 1));
    node->right = std::unique_ptr<KDNode>(buildBalancedTree(photons, mid + 1, end, depth + 1));
    return node.release();
}

void PhotonKDTree::getPhotonsRecursive(const KDNode* node, std::vector<Photon>& photons) const {
    if (!node) return;
    photons.push_back(node->photon);
    getPhotonsRecursive(node->left.get(), photons);
    getPhotonsRecursive(node->right.get(), photons);
}

void PhotonKDTree::findNearestNeighborsRecursive(
    const KDNode* node, const Vec3& point, const float maxDistSq, const int depth,
    std::vector<std::pair<float, Photon>>& nearestPhotons)
const {
    if (!node) return;

    // Calculate the squared distance
    float distSq = node->photon.position.distanceSquared(point);

    // Add the photon if within the distance
    if (distSq < maxDistSq) {
        nearestPhotons.emplace_back(distSq, node->photon);
    }

    // Determine the axis and the difference
    const int axis = depth % 3;
    const float diff = point[axis] - node->photon.position[axis];

    // Determine which subtree to explore first
    const KDNode* nearChild = diff < 0 ? node->left.get() : node->right.get();
    const KDNode* farChild = diff < 0 ? node->right.get() : node->left.get();

    // Search in the nearer subtree
    findNearestNeighborsRecursive(nearChild, point, maxDistSq, depth + 1, nearestPhotons);

    // If the distance to the axis is less than the maximum distance, explore the other subtree
    if (diff * diff < maxDistSq) {
        findNearestNeighborsRecursive(farChild, point, maxDistSq, depth + 1, nearestPhotons);
    }
}