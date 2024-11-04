#ifndef PHOTONKD_TREE_H
#define PHOTONKD_TREE_H

#include "Vec3.h"
#include "Photon.h"

#include <vector>
#include <memory>

// -------------------------------------------
// PhotonKDNode Structure
// -------------------------------------------

/**
 * Structure representing a node in the KD tree.
 */
struct PhotonKDNode {
    Photon data; ///< Photon data contained in the node.
    std::unique_ptr<PhotonKDNode> left; ///< Pointer to the left child node.
    std::unique_ptr<PhotonKDNode> right; ///< Pointer to the right child node.

    /**
     * Constructor for PhotonKDNode.
     * @param data Photon data to be stored in the node.
     */
    explicit PhotonKDNode(const Photon &data) : data(data) {}
};

// -------------------------------------------
// PhotonKDTree Class
// -------------------------------------------

/**
 * Class representing the KD tree for photons.
 */
class PhotonKDTree {
public:
    PhotonKDTree() = default;
    explicit PhotonKDTree(std::vector<Photon>& elements) {
        root = std::unique_ptr<PhotonKDNode>(buildPhotonBalancedTree(elements, 0, static_cast<int>(elements.size()), 0));
    }

    /**
     * Checks if the KD tree is empty.
     * @return True if the tree is empty, false otherwise.
     */
    [[nodiscard]] bool isEmpty() const { return !root; }

    /**
     * Converts the KD tree to a vector of photons.
     * @return Vector of photons contained in the tree.
     */
    [[nodiscard]] std::vector<Photon> toVector() const;

    /**
     * Finds the nearest neighbors of a given point within a specified maximum distance.
     * @param point The point to find the nearest neighbors for.
     * @param maxDistance The maximum distance to search for neighbors.
     * @return Vector of nearest photons within the specified distance.
     */
    [[nodiscard]] std::vector<Photon> findNearestNeighbors(const Vec3 &point, float maxDistance) const;

private:
    std::unique_ptr<PhotonKDNode> root; ///< Root node of the KD tree.

    /**
     * Recursively retrieves elements from the KD tree.
     * @param node The current node.
     * @param elements Vector to store the retrieved elements.
     */
    static void getElementsRecursive(const PhotonKDNode *node, std::vector<Photon> &elements);

    /**
     * Builds a balanced KD tree from a vector of photons.
     * @param elements Vector of photons to build the tree from.
     * @param start Start index of the current segment.
     * @param end End index of the current segment.
     * @param depth Current depth of the tree.
     * @return Pointer to the root node of the built KD tree.
     */
    PhotonKDNode *buildPhotonBalancedTree(std::vector<Photon> &elements, int start, int end, int depth);

    /**
     * Recursively finds the nearest neighbors of a given point within a specified maximum distance.
     * @param node The current node.
     * @param point The point to find the nearest neighbors for.
     * @param maxDistSq The maximum distance squared to search for neighbors.
     * @param depth Current depth of the tree.
     * @param nearestElements Vector to store the nearest neighbors found.
     */
    void findNearestNeighborsRecursive(const PhotonKDNode *node, const Vec3 &point, float maxDistSq, int depth, std::vector<std::pair<float, Photon>> &nearestElements) const;
};

#endif // PHOTONKD_TREE_H