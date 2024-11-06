#ifndef PHOTONMAP_H
#define PHOTONMAP_H

#include <random>
#include <vector>

#include "Light.h"
#include "Mesh.h"
#include "MeshKDTree.h"
#include "Photon.h"
#include "PhotonKDTree.h"
#include "Sphere.h"
#include "Square.h"

// -------------------------------------------
// PhotonMap Class
// -------------------------------------------

/**
 * Class representing a photon map.
 */
class PhotonMap {
public:
    PhotonKDTree globalPhotonTree;      ///< KD tree for global photons.
    PhotonKDTree causticsPhotonTree;    ///< KD tree for caustics photons.

    /**
     * Emits photons from the given lights and objects.
     * @param lights Vector of lights to emit photons from.
     * @param spheres Vector of spheres in the scene.
     * @param squares Vector of squares in the scene.
     * @param meshes Vector of meshes in the scene.
     * @param kdTree KD tree of the meshes.
     * @param photons Number of photons to emit.
     */
    void emitPhotons(
        const std::vector<Light> &lights, const std::vector<Sphere> &spheres, const std::vector<Square> &squares,
        const std::vector<Mesh> &meshes, const MeshKDTree &kdTree, int photons
    );

    /**
     * Computes caustics at the given position based on the material.
     * @param position The position to compute caustics at.
     * @param material
     * @return The computed caustics as a Vec3.
     */
    [[nodiscard]] Vec3 computeCaustics(const Vec3 &position, const Material &material) const;

    /**
     * Draws the photons for debugging purposes.
     */
    void debugDrawPhotons();

    /**
     * Draws the light paths for debugging purposes.
     */
    void debugDrawLightPaths() const;

private:
    std::vector<Photon> initialPhotons;     ///< Vector of initial photons.
    std::vector<Photon> globalPhotons;      ///< Vector of global photons.

    /**
     * Emits photons for a specific thread.
     * @param photonCount Number of photons to emit.
     * @param lights Vector of lights to emit photons from.
     * @param spheres Vector of spheres in the scene.
     * @param squares Vector of squares in the scene.
     * @param meshes Vector of meshes in the scene.
     * @param kdTree KD tree of the meshes.
     * @param photonMutex Mutex for synchronizing photon emission.
     * @param globalPhotons Vector to store emitted mirror photons.
     * @param photonsToEmit Vector to store photons to be emitted.
     */
    void emitPhotonsForThread(
        int photonCount, const std::vector<Light> &lights, const std::vector<Sphere> &spheres, const std::vector<Square> &squares,
        const std::vector<Mesh> &meshes, const MeshKDTree &kdTree, std::mutex &photonMutex, std::vector<Photon> &globalPhotons, std::vector<Photon> &photonsToEmit
    );

    /**
     * Generates a random direction for photon emission.
     * @param rng Random number generator.
     * @param normal Normal at the point of emission. Default is Vec3(0.0f).
     * @return A random direction as a Vec3.
     */
    static Vec3 randomDirection(std::mt19937 &rng, const Vec3 &normal = Vec3(0.0f));

    /**
     * Generates a random direction towards objects.
     * @param rng Random number generator.
     * @param spheres Vector of spheres in the scene.
     * @param meshes Vector of meshes in the scene.
     * @return A random direction as a Vec3.
     */
    static Vec3 randomDirectionTowardsObjects(std::mt19937 &rng, const std::vector<Sphere> &spheres,const std::vector<Mesh> &meshes);

    /**
     * Creates an initial photon from the given lights and objects.
     * @param lights Vector of lights to emit photons from.
     * @param rng Random number generator.
     * @return The created initial photon.
     */
    static Photon createInitialPhoton(const std::vector<Light> &lights, std::mt19937 &rng);

    /**
     * Computes the probabilities of reflection, diffuse reflection, and specular reflection.
     * @param material The material at the point of interaction.
     * @param Pr Probability of reflection.
     * @param Pd Probability of diffuse reflection.
     * @param Ps Probability of specular reflection.
     */
    static void computeReflectionProbabilities(const Material &material, float &Pr, float &Pd, float &Ps);

    /**
     * Processes a glass photon.
     * @param photon The photon to process.
     * @param normal The normal at the point of interaction.
     * @param material The material at the point of interaction.
     * @param photons Vector to store the processed photons.
     */
    static void processGlassPhoton(Photon &photon, const Vec3 &normal, const Material &material, std::vector<Photon> &photons);

    /**
     * Processes a mirror photon.
     * @param photon The photon to process.
     * @param normal The normal at the point of interaction.
     * @param material The material at the point of interaction.
     */
    static void processMirrorPhoton(Photon &photon, const Vec3 &normal, const Material &material);
};

#endif //PHOTONMAP_H