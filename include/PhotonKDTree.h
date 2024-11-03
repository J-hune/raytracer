#ifndef PHOTONKD_TREE_H
#define PHOTONKD_TREE_H

#include "Vec3.h"
#include "Photon.h"

#include <vector>
#include <memory>

// Represents a node in the KD tree
struct PhotonKDNode {
    Photon data;
    std::unique_ptr<PhotonKDNode> left;
    std::unique_ptr<PhotonKDNode> right;

    explicit PhotonKDNode(const Photon &data) : data(data) {}
};

// Represents the KD tree
class PhotonKDTree {
public:
    PhotonKDTree() = default;
    explicit PhotonKDTree(std::vector<Photon>& elements) {
        root = std::unique_ptr<PhotonKDNode>(buildPhotonBalancedTree(elements, 0, static_cast<int>(elements.size()), 0));
    }

    [[nodiscard]] bool isEmpty() const { return !root; }
    [[nodiscard]] std::vector<Photon> toVector() const;
    [[nodiscard]] std::vector<Photon> findNearestNeighbors(const Vec3 &point, float maxDistance) const;

private:
    std::unique_ptr<PhotonKDNode> root;
    static void getElementsRecursive(const PhotonKDNode *node, std::vector<Photon> &elements);
    PhotonKDNode *buildPhotonBalancedTree(std::vector<Photon> &elements, int start, int end, int depth);
    void findNearestNeighborsRecursive(const PhotonKDNode *node, const Vec3 &point, float maxDistSq, int depth, std::vector<std::pair<float, Photon>> &nearestElements) const;
};

#endif // PHOTONKD_TREE_H