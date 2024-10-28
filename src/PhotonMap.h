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
    PhotonKDTree photonGlassTree;
    PhotonKDTree photonMirrorTree;
    PhotonKDTree photonDiffuseTree;

    void emitPhotons(const std::vector<Light> &lights, const std::vector<Sphere> &spheres, const std::vector<Square> &squares, const std::vector<Mesh> &meshes, int photons);
    [[nodiscard]] Vec3 renderCaustics(const Vec3 &position, const Vec3 &normal, const Material &material) const;
    [[nodiscard]] std::vector<Photon> findNearestGlassPhotons(const Vec3 &point, float maxDistance) const;
    [[nodiscard]] std::vector<Photon> findNearestMirrorPhotons(const Vec3 &point, float maxDistance) const;
    [[nodiscard]] std::vector<Photon> findNearestDiffusePhotons(const Vec3 &point, float maxDistance) const;
    void debugDrawPhotons() const;

private:
    static Vec3 randomDirection(std::mt19937 &rng);
    static Vec3 calculateCausticsContribution(std::vector<Photon> &photons, const Vec3 &position, const Vec3 &normal, const Material &material, float maxRadius, float sigma) ;
};

#endif //PHOTONMAP_H
