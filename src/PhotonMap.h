#ifndef PHOTONMAP_H
#define PHOTONMAP_H

#include <vector>
#include <random>


#include "PhotonKDTree.h"
#include "Light.h"
#include "Sphere.h"
#include "Square.h"
#include "Mesh.h"

class PhotonMap {
public:
    PhotonKDTree mirrorPhotonTree;
    PhotonKDTree glassPhotonTree;

    void emitPhotons(const std::vector<Light> &lights, const std::vector<Sphere> &spheres, const std::vector<Square> &squares, const std::vector<Mesh> &meshes, int photons);
    [[nodiscard]] Vec3 computeCaustics(const Vec3 &position, const Material &material) const;
    void debugDrawPhotons() const;

private:
    static Vec3 randomDirection(std::mt19937 &rng);
    static Photon createInitialPhoton(const std::vector<Light> &lights, std::mt19937 &rng);
    static void processGlassPhoton(Photon &photon, const Vec3 &normal, const Material &material, std::vector<Photon> &photons);
    static void processMirrorPhoton(Photon &photon, const Vec3 &normal, const Material &material);
};

#endif //PHOTONMAP_H
