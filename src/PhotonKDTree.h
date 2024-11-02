#ifndef PHOTONKDTREE_H
#define PHOTONKDTREE_H

#include "Photon.h"
#include "Vec3.h"
#include <vector>
#include <memory>
#include <queue>

// Represents a node in the KD tree for photons
struct KDNode {
    Photon photon;
    std::unique_ptr<KDNode> left;
    std::unique_ptr<KDNode> right;

    explicit KDNode(Photon photon);
};

// Represents the KD tree for photons
class PhotonKDTree {
public:
    PhotonKDTree() = default;
    explicit PhotonKDTree(std::vector<Photon>& photons);

    [[nodiscard]] std::vector<Photon> toVector() const;
    [[nodiscard]] std::vector<Photon> findNearestNeighbors(const Vec3 &point, float maxDistance) const;

private:
    std::unique_ptr<KDNode> root;

    KDNode* buildBalancedTree(std::vector<Photon>& photons, int start, int end, int depth);
    void getPhotonsRecursive(const KDNode* node, std::vector<Photon>& photons) const;
    void findNearestNeighborsRecursive(const KDNode *node, const Vec3 &point, float maxDistSq, int depth,
                                       std::vector<std::pair<float, Photon>> &nearestPhotons) const;
};

#endif // PHOTONKDTREE_H