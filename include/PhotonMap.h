#ifndef PHOTONMAP_H
#define PHOTONMAP_H

#include <random>
#include <vector>

#include "Light.h"
#include "Mesh.h"
#include "MeshKDTree.h"
#include "Photon.h"
#include "PhotonKDTree.h"
#include "Settings.h"
#include "Sphere.h"
#include "Square.h"

// -------------------------------------------
// PhotonMap Class
// -------------------------------------------

struct RaySceneIntersection;
/**
 * Class representing a photon map.
 */
class PhotonMap {
public:
    PhotonKDTree globalPhotonTree; ///< KD tree for global photons.
    PhotonKDTree causticsPhotonTree; ///< KD tree for caustics photons.

    /**
     * Emits photons from the given lights and objects.
     * @param lights Vector of lights to emit photons from.
     * @param spheres Vector of spheres in the scene.
     * @param squares Vector of squares in the scene.
     * @param meshes Vector of meshes in the scene.
     * @param kdTree KD tree of the meshes.
     * @param settings The settings of the application.
     */
    void emitPhotons(const std::vector<Light> &lights, const std::vector<Sphere> &spheres,
                     const std::vector<Square> &squares, const std::vector<Mesh> &meshes, const MeshKDTree &kdTree,
                     const Settings &settings
    );

    /**
     * Computes the caustics effect at a given position.
     * @param position The position to compute the caustics effect.
     * @param material The material at the position.
     * @return The color of the caustics effect.
     */
    [[nodiscard]] Vec3 computeCaustics(const Vec3 &position, const Material &material) const;

    /**
     * Draws the photons for debugging purposes.
     * @param type The type of photons to draw. 1: Caustics, 2: Global + Caustics, 3: Initial + Global + Caustics.
     */
    void debugDrawPhotons(int type);

    /**
     * Draws the light paths for debugging purposes.
     */
    void debugDrawLightPaths() const;

private:
    std::vector<Photon> initialGlobalPhotons; ///< Vector of initial global photons.
    std::vector<Photon> initialCausticsPhotons; ///< Vector of initial caustics photons.
    std::vector<Photon> globalPhotons; ///< Vector of global photons.
    std::vector<Photon> causticsPhotons; ///< Vector of caustics photons.

    /**
     * Emits photons from the given lights and objects.
     * @param lights Vector of lights to emit photons from.
     * @param spheres Vector of spheres in the scene.
     * @param squares Vector of squares in the scene.
     * @param meshes Vector of meshes in the scene.
     * @param kdTree KD tree of the meshes.
     * @param photonCount The number of photons to emit.
     * @param photonType The type of photons to emit.
     * @return The emitted photons as a KD tree.
     */
    PhotonKDTree emitPhotonsWithType(const std::vector<Light> &lights, const std::vector<Sphere> &spheres,
        const std::vector<Square> &squares, const std::vector<Mesh> &meshes,
        const MeshKDTree &kdTree, int photonCount, int photonType
    );

    /**
     * Processes the path of a photon.
     * @param photon The photon to process.
     * @param spheres Vector of spheres in the scene.
     * @param squares Vector of squares in the scene.
     * @param meshes Vector of meshes in the scene.
     * @param kdTree KD tree of the meshes.
     * @param dist Uniform real distribution for random numbers.
     * @param rng Random number generator.
     * @param localPhotonsToEmit
     * @param localPhotons Vector to store the processed photons.
     * @param photonType The type of photons to process.
     */
    void processPhotonPath(Photon &photon, const std::vector<Sphere> &spheres, const std::vector<Square> &squares,
                           const std::vector<Mesh> &meshes, const MeshKDTree &kdTree,
                           std::uniform_real_distribution<float> &dist, std::mt19937 &rng,
                           std::vector<Photon> &localPhotonsToEmit, std::vector<Photon> &localPhotons, int photonType);

    /**
     * Processes the interaction of a photon.
     * @param photon The photon to process.
     * @param intersection The intersection details.
     * @param dist Uniform real distribution for random numbers.
     * @param rng Random number generator.
     * @param localPhotonsToEmit
     * @param localPhotons Vector to store the processed photons.
     * @param photonType The type of photons to process.
     * @param absorbed Flag to indicate if the photon was absorbed.
     */
    void processPhotonInteraction(Photon &photon, const RaySceneIntersection &intersection,
                                  std::uniform_real_distribution<float> &dist, std::mt19937 &rng,
                                  std::vector<Photon> &localPhotonsToEmit, std::vector<Photon> &localPhotons,
                                  int photonType, bool &absorbed);

    /**
     * Processes a diffuse photon.
     * @param photon The photon to process.
     * @param intersection The intersection details.
     * @param dist Uniform real distribution for random numbers.
     * @param rng Random number generator.
     * @param localPhotons Vector to store the processed photons.
     * @param photonType The type of photons to process.
     * @param absorbed Flag to indicate if the photon was absorbed.
     */
    static void processDiffusePhoton(Photon &photon, const RaySceneIntersection &intersection,
        std::uniform_real_distribution<float> &dist, std::mt19937 &rng,
        std::vector<Photon> &localPhotons, int photonType, bool &absorbed
    );

    /**
     * Normalizes the colors of the photons.
     * @param photons Vector of photons to normalize.
     * @param emittedPhotons The number of emitted photons.
     */
    static void normalizePhotonColors(std::vector<Photon> &photons, int emittedPhotons);

    /**
     * Emits photons for a specific thread.
     * @param photonCount The number of photons to emit.
     * @param photonType The type of photons to emit.
     * @param lights Vector of lights to emit photons from.
     * @param spheres Vector of spheres in the scene.
     * @param squares Vector of squares in the scene.
     * @param meshes Vector of meshes in the scene.
     * @param kdTree KD tree of the meshes.
     * @param photonMutex Mutex for synchronizing photon emission.
     * @param photons Vector to store emitted photons.
     * @param photonsToEmit Vector to store photons to be emitted.
     */
    void emitPhotonsForThread(int photonCount, int photonType, const std::vector<Light> &lights,
                              const std::vector<Sphere> &spheres, const std::vector<Square> &squares,
                              const std::vector<Mesh> &meshes, const MeshKDTree &kdTree, std::mutex &photonMutex,
                              std::vector<Photon> &photons, std::vector<Photon> &photonsToEmit);

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
    static Vec3 randomDirectionTowardsObjects(std::mt19937 &rng, const std::vector<Sphere> &spheres,
                                              const std::vector<Mesh> &meshes);

    /**
     * Creates an initial photon from the given lights and objects.
     * @param lights Vector of lights to emit photons from.
     * @param rng Random number generator.
     * @return The created initial photon.
     */
    static Photon createInitialPhoton(const std::vector<Light> &lights, std::mt19937 &rng);

    /**
     * Creates an initial photon towards objects from the given lights and
     * objects.
     * @param lights Vector of lights to emit photons from.
     * @param rng Random number generator.
     * @param spheres Vector of spheres in the scene.
     * @param meshes Vector of meshes in the scene.
     * @return The created initial photon.
     */
    static Photon createInitialPhotonTowardsObjects(const std::vector<Light> &lights, std::mt19937 &rng,
                                                    const std::vector<Sphere> &spheres,
                                                    const std::vector<Mesh> &meshes);

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
    static void processGlassPhoton(Photon &photon, const Vec3 &normal, const Material &material,
                                   std::vector<Photon> &photons);

    /**
     * Processes a mirror photon.
     * @param photon The photon to process.
     * @param normal The normal at the point of interaction.
     * @param material The material at the point of interaction.
     */
    static void processMirrorPhoton(Photon &photon, const Vec3 &normal, const Material &material);
};

#endif // PHOTONMAP_H
