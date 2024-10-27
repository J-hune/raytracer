#ifndef PHOTONKDTREE_H
#define PHOTONKDTREE_H

#include <utility>
#include <vector>
#include <algorithm>
#include <memory>

#include "Photon.h"
#include "Vec3.h"

// Class representing a node in the KD tree
struct KDNode {
    Photon photon;
    std::unique_ptr<KDNode> left;
    std::unique_ptr<KDNode> right;

    explicit KDNode(Photon photon) : photon(std::move(photon)), left(nullptr), right(nullptr) {}
};

// Class representing the KD tree for photons
class PhotonKDTree {
public:
    PhotonKDTree() = default;
    explicit PhotonKDTree(std::vector<Photon>& photons) {
        root = std::unique_ptr<KDNode>(buildBalancedTree(photons, 0, photons.size(), 0));
    }

    // Convert the KD tree to a vector of photons
    std::vector<Photon> toVector() const {
        std::vector<Photon> photons;
        getPhotonsRecursive(root.get(), photons);
        return photons;
    }

    // Find the nearest neighbors to a given point within a maximum distance
    std::vector<Photon> findNearestNeighbors(const Vec3& point, const float maxDistance) const {
        std::vector<std::pair<float, Photon>> nearestPhotons;

        // Call the recursive function
        findNearestNeighborsRecursive(root.get(), point, maxDistance * maxDistance, 0, nearestPhotons);

        // Sort the photons by squared distance and extract only the photons
        std::ranges::sort(
            nearestPhotons,
            [](const std::pair<float, Photon>& a, const std::pair<float, Photon>& b) {
                return a.first < b.first;
            }
        );

        // Keep only the photons
        std::vector<Photon> result;
        for (const auto &[fst, snd] : nearestPhotons) {
            result.push_back(snd);
        }
        return result;
    }

private:
    std::unique_ptr<KDNode> root;

    // Build a balanced KD tree from a vector of photons
    KDNode* buildBalancedTree(std::vector<Photon>& photons, const int start, const int end, const int depth) {
        if (start >= end) return nullptr;

        int axis = depth % 3; // Choose the axis based on the depth
        const int median = (start + end) / 2;

        // Use std::nth_element to place the median at its correct position
        std::nth_element(
            photons.begin() + start, photons.begin() + median, photons.begin() + end,
            [axis](const Photon& a, const Photon& b) {
                return a.position[axis] < b.position[axis];
        });

        auto node = std::make_unique<KDNode>(photons[median]);
        node->left = std::unique_ptr<KDNode>(buildBalancedTree(photons, start, median, depth + 1));
        node->right = std::unique_ptr<KDNode>(buildBalancedTree(photons, median + 1, end, depth + 1));

        return node.release();
    }

    // Recursively get photons from the KD tree
    void getPhotonsRecursive(const KDNode* node, std::vector<Photon>& photons) const {
        if (!node) return;
        photons.push_back(node->photon);
        getPhotonsRecursive(node->left.get(), photons);
        getPhotonsRecursive(node->right.get(), photons);
    }

    // Recursively find the nearest neighbors
    void findNearestNeighborsRecursive(
        const KDNode* node,
        const Vec3& point,
        const float maxDistSq,
        const int depth,
        std::vector<std::pair<float, Photon>>& nearestPhotons) const
    {
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
        const KDNode* nearer = (diff < 0) ? node->left.get() : node->right.get();
        const KDNode* further = (diff < 0) ? node->right.get() : node->left.get();

        // Search in the nearer subtree
        findNearestNeighborsRecursive(nearer, point, maxDistSq, depth + 1, nearestPhotons);

        // If the distance to the axis is less than the maximum distance, explore the other subtree
        if (diff * diff < maxDistSq) {
            findNearestNeighborsRecursive(further, point, maxDistSq, depth + 1, nearestPhotons);
        }
    }
};

#endif // PHOTONKDTREE_H