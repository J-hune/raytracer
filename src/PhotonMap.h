#ifndef PHOTONMAP_H
#define PHOTONMAP_H

#include <vector>
#include <random>

#include "PhotonKDTree.h"
#include "Photon.h"
#include "Light.h"
#include "Sphere.h"
#include "Square.h"
#include "Mesh.h"
#include "Intersection.h"
#include "Lighting.h"

class PhotonMap {
public:
    PhotonKDTree photonGlassTree;
    PhotonKDTree photonMirrorTree;
    PhotonKDTree photonDiffuseTree;

    void emitPhotons(const std::vector<Light> &lights, const std::vector<Sphere> &spheres, const std::vector<Square> &squares, const std::vector<Mesh> &meshes, int photons);
    Vec3 renderCaustics(const Vec3 &position, const Vec3 &normal, const Material &material) const;
    void debugDrawPhotons() const;
    std::vector<Photon> findNearestGlassPhotons(const Vec3 &point, float maxDistance) const;
    std::vector<Photon> findNearestMirrorPhotons(const Vec3 &point, float maxDistance) const;
    std::vector<Photon> findNearestDiffusePhotons(const Vec3 &point, float maxDistance) const;

private:
    static Vec3 randomDirection(std::mt19937 &rng);
};

#endif //PHOTONMAP_H
