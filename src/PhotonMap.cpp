#include "../include/PhotonMap.h"
#include "../include/Intersection.h"
#include "../include/Lighting.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <thread>
#include <vector>
#include <GL/gl.h>

/******************************************************************************************************************/
/******************************************** PHOTON EMISSION INTERFACE *******************************************/
/******************************************************************************************************************/

void PhotonMap::emitPhotons(
    const std::vector<Light>& lights, const std::vector<Sphere>& spheres, const std::vector<Square>& squares,
    const std::vector<Mesh>& meshes, const MeshKDTree& kdTree, const Settings& settings)
{
    if (lights.empty()) {
        throw std::invalid_argument("No lights in the scene.");
    }

    /* Emit Global Photons */
    globalPhotonTree = emitPhotonsWithType(lights, spheres, squares, meshes, kdTree, settings.globalPhotons, 0);

    /* Emit Caustics Photons */
    causticsPhotonTree = emitPhotonsWithType(lights, spheres, squares, meshes, kdTree, settings.causticsPhotons, 1);
}

PhotonKDTree PhotonMap::emitPhotonsWithType(
    const std::vector<Light>& lights, const std::vector<Sphere>& spheres, const std::vector<Square>& squares,
    const std::vector<Mesh>& meshes, const MeshKDTree& kdTree, const int photonCount, int photonType)
{
    std::vector<Photon> photons;
    std::mutex photonMutex;

    const int numThreads = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
    const int photonsPerThread = std::max(1, photonCount / numThreads);

    // One counter slot per thread. A single shared counter written by every thread would be a data race,
    // and the total is needed below to normalise the photon energy.
    std::vector<int> emittedPerThread(static_cast<size_t>(numThreads), 0);

    std::vector<std::thread> threads;
    threads.reserve(static_cast<size_t>(numThreads));
    for (int i = 0; i < numThreads; ++i) {
        threads.emplace_back([this, i, &emittedPerThread, photonsPerThread, photonType,
                              &lights, &spheres, &squares, &meshes, &kdTree, &photonMutex, &photons] {
            emittedPerThread[static_cast<size_t>(i)] = emitPhotonsForThread(
                photonsPerThread, photonType, lights, spheres, squares, meshes, kdTree, photonMutex, photons);
        });
    }

    for (auto& thread : threads) {
        thread.join();
    }

    // A stored photon must carry the light power divided by the number of photons actually emitted.
    // Without this normalisation the illumination scales with the photon count, so emitting more
    // photons brightens the image instead of reducing its noise.
    const int totalEmitted = std::accumulate(emittedPerThread.begin(), emittedPerThread.end(), 0);
    if (totalEmitted > 0) {
        const auto emittedCount = static_cast<float>(totalEmitted);
        for (Photon& photon : photons) {
            photon.color /= emittedCount;
        }
    }

    return PhotonKDTree(photons);
}

/******************************************************************************************************************/
/******************************************* PER-THREAD PHOTON EMISSION *******************************************/
/******************************************************************************************************************/

int PhotonMap::emitPhotonsForThread(
    const int photonCount, const int photonType,
    const std::vector<Light>& lights, const std::vector<Sphere>& spheres, const std::vector<Square>& squares, const std::vector<Mesh>& meshes, const MeshKDTree& kdTree,
    std::mutex& photonMutex, std::vector<Photon>& photons)
{
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution dist(0.0f, 1.0f);

    std::vector<Photon> localPhotons;
    std::vector<Photon> localPhotonsToEmit;
    std::vector<Photon> localInitialPhotons;
    localPhotonsToEmit.reserve(64);

    // Guard against a pathological chain of specular splits. Energy decay already bounds the depth,
    // this only keeps a degenerate scene from growing the queue without limit.
    constexpr size_t maxQueuedSecondaries = 4096;

    int emittedPhotons = 0;
    // Loop on the number of photons *emitted*, not on the number stored. The energy of each photon is
    // normalised by the emitted count, so that is the quantity that has to be controlled here.
    while (emittedPhotons < photonCount) {
        Photon photon = photonType == 0 ? createInitialPhoton(lights, rng) : createInitialPhotonTowardsObjects(lights, rng, spheres, meshes);
        localInitialPhotons.emplace_back(photon);
        emittedPhotons++;
        processPhotonPath(photon, spheres, squares, meshes, kdTree, dist, rng, localPhotonsToEmit, localPhotons, photonType);

        // Trace the secondary photons produced by specular splitting, i.e. the Fresnel-reflected half
        // of a glass interaction. This queue used to be filled and never read, so the reflected part of
        // every glass hit was silently discarded and glass contributed refraction only.
        while (!localPhotonsToEmit.empty()) {
            if (localPhotonsToEmit.size() > maxQueuedSecondaries) {
                localPhotonsToEmit.clear();
                break;
            }
            Photon secondary = localPhotonsToEmit.back();
            localPhotonsToEmit.pop_back();
            processPhotonPath(secondary, spheres, squares, meshes, kdTree, dist, rng, localPhotonsToEmit, localPhotons, photonType);
        }
    }

    std::lock_guard lock(photonMutex);
    photons.insert(photons.end(), localPhotons.begin(), localPhotons.end());
    if (photonType == 0) {
        initialGlobalPhotons.insert(initialGlobalPhotons.end(), localInitialPhotons.begin(), localInitialPhotons.end());
    } else {
        initialCausticsPhotons.insert(initialCausticsPhotons.end(), localInitialPhotons.begin(), localInitialPhotons.end());
    }

    return emittedPhotons;
}

void PhotonMap::processPhotonPath(
    Photon& photon, const std::vector<Sphere>& spheres, const std::vector<Square>& squares, const std::vector<Mesh>& meshes, const MeshKDTree& kdTree,
    std::uniform_real_distribution<float>& dist, std::mt19937& rng, std::vector<Photon>& localPhotonsToEmit, std::vector<Photon>& localPhotons, const int photonType)
{
    bool absorbed = false;
    int bounces = 0;
    constexpr int maxBounces = 500;

    while (!absorbed && bounces < maxBounces) {
        if (photon.color.length() < 0.001f) break; // Stop if the photon has low energy

        Ray ray(photon.position, photon.direction);
        RaySceneIntersection intersection = Intersection::computeIntersection(ray, spheres, squares, meshes, kdTree, 0.0f);
        if (!intersection.intersectionExists) break;

        photon.position = intersection.intersection;
        processPhotonInteraction(photon, intersection, dist, rng, localPhotonsToEmit, localPhotons, photonType, absorbed, bounces);
        bounces++;
    }
}

void PhotonMap::processPhotonInteraction(
    Photon& photon, const RaySceneIntersection& intersection, std::uniform_real_distribution<float>& dist, std::mt19937& rng,
    std::vector<Photon>& localPhotonsToEmit, std::vector<Photon>& localPhotons, const int photonType, bool& absorbed, const int bounces)
{
    switch (intersection.material.type) {
        case Material_Glass:
            processGlassPhoton(photon, intersection.normal, intersection.material, localPhotonsToEmit);
            break;
        case Material_Mirror:
            processMirrorPhoton(photon, intersection.normal, intersection.material);
            break;
        case Material_Diffuse_Blinn_Phong:
            processDiffusePhoton(photon, intersection, dist, rng, localPhotons, photonType, absorbed, bounces);
            break;
    }
}


/******************************************************************************************************************/
/******************************************** PHOTON PROCESSES PER MATERIAL ***************************************/
/******************************************************************************************************************/

void PhotonMap::processGlassPhoton(Photon& photon, const Vec3& normal, const Material& material, std::vector<Photon>& photons) {
    const float fresnelEffect = Lighting::computeFresnelEffect(photon.direction, normal, material.index_medium);

    // Two effects: one part goes into reflection, the other into refraction
    // The percentage of reflection is calculated by the Fresnel effect
    const Vec3 refractedDirection = Lighting::computeRefractedDirection(photon.direction, normal, material.index_medium).normalize();
    const Vec3 reflectedDirection = Lighting::computeReflectedDirection(photon.direction, normal).normalize();

    // The color of the reflected photon takes into account the color of the object and the percentage of reflection
    Photon reflectedPhoton = photon;
    reflectedPhoton.direction = reflectedDirection;
    reflectedPhoton.color = fresnelEffect * Vec3::compProduct(photon.color, material.specular_material);
    reflectedPhoton.hasSpecularBounce = true;
    photons.emplace_back(reflectedPhoton);

    // Add object color to the photon color (with transparency)
    photon.direction = refractedDirection;
    photon.color = (1 - fresnelEffect) * (material.transparency * photon.color + (1 - material.transparency) * material.specular_material);
    photon.hasSpecularBounce = true;
}

void PhotonMap::processMirrorPhoton(Photon& photon, const Vec3& normal, const Material& material) {
    // Change the direction of the photon to make it bounce, keep the color at 90% to simulate reflection
    photon.direction = Lighting::computeReflectedDirection(photon.direction, normal).normalize();
    photon.color = Vec3::compProduct(photon.color, material.specular_material);
    photon.hasSpecularBounce = true;
}

void PhotonMap::processDiffusePhoton(
    Photon& photon, const RaySceneIntersection& intersection, std::uniform_real_distribution<float>& dist, std::mt19937& rng,
    std::vector<Photon>& localPhotons, const int photonType, bool& absorbed, const int bounces)
{
    float Pr, Pd, Ps;
    computeReflectionProbabilities(intersection.material, Pr, Pd, Ps);

    // Deposit the power that actually reaches the surface, before any Russian roulette compensation.
    // The two storing branches used to deposit different quantities: the diffuse one after dividing by
    // Pd, the absorbing one without. Two identical photons therefore stored energies differing by a
    // factor 1/Pd depending on a coin flip. The compensation belongs to the continuing path, not to the
    // deposit.
    //
    // The global map takes every diffuse deposit after the first bounce. The caustics map only takes
    // photons that reached this surface through a mirror or a transparent object, that is an L(S|D)*S+D
    // path. Without that condition purely diffuse paths ended up in the caustics map, where they add a
    // low-frequency haze that a small gather radius cannot resolve.
    if (bounces > 0 && (photonType == 0 || photon.hasSpecularBounce)) {
        localPhotons.emplace_back(photon);
    }

    // A caustics photon stops at the first diffuse surface: what follows is diffuse interreflection,
    // which belongs to the global map.
    if (photonType == 1) {
        absorbed = true;
        return;
    }

    const float xi = dist(rng);

    if (xi < Pd) { // Diffuse reflection ξ ∈ [0, Pd]
        photon.direction = randomDirection(rng, intersection.normal);
        photon.color = Vec3::compProduct(photon.color, intersection.material.diffuse_material) / Pd;
    } else if (xi < Pd + Ps) { // Specular reflection ξ ∈ ]Pd, Pd + Ps]
        // Compensate by the probability of the event actually drawn (Ps), not by the total reflection
        // probability Pr = Pd + Ps, which is larger and therefore drained energy from every specular path.
        photon.direction = Lighting::computeReflectedDirection(photon.direction, intersection.normal).normalize();
        photon.color = Vec3::compProduct(photon.color, intersection.material.specular_material) / Ps;
    } else { // Absorption ξ ∈ ]Pd + Ps, 1]
        absorbed = true;
    }
}


/******************************************************************************************************************/
/******************************************** RENDERING PHOTON PATHS **********************************************/
/******************************************************************************************************************/

namespace {
    /**
     * @brief Cone-filtered density estimate over a set of gathered photons.
     *
     * Two corrections with respect to a naive estimate:
     *
     * - The normalisation uses the distance to the farthest photon actually gathered, not the configured
     *   maximum radius. In a sparse region the search returns far fewer photons than requested, and
     *   dividing by the full disc area would darken that region for no physical reason. This is what
     *   produced blotches and dark corners.
     * - It does not divide by the number of photons gathered. Dividing the sum of powers by the disc area
     *   is already what turns a flux into a density; dividing a second time turns it into a mean per
     *   photon, which is dimensionally wrong and makes the luminance depend on a quality setting.
     *
     * @param photons The photons gathered around the shading point.
     * @param point The point being shaded.
     * @param albedo The diffuse albedo of the surface at that point.
     * @param k Cone filter constant.
     * @param maxRadius Upper bound on the gather radius, from the settings.
     * @return The estimated indirect radiance.
     */
    Vec3 coneFilteredEstimate(const std::vector<Photon> &photons, const Vec3 &point, const Vec3 &albedo,
                              const float k, const float maxRadius) {
        if (photons.empty()) return {0.0f, 0.0f, 0.0f};

        // Effective radius: the distance to the farthest photon retained by the search.
        float radiusSq = 0.0f;
        for (const Photon &photon: photons) {
            radiusSq = std::max(radiusSq, (photon.position - point).squareLength());
        }
        const float radius = std::min(std::sqrt(radiusSq), maxRadius);
        if (radius <= 1e-6f) return {0.0f, 0.0f, 0.0f};

        Vec3 sum(0.0f, 0.0f, 0.0f);
        for (const Photon &photon: photons) {
            const float dp = (photon.position - point).length();
            const float wpc = 1.0f - dp / (k * radius);
            if (wpc <= 0.0f) continue;
            sum += wpc * Vec3::compProduct(photon.color, albedo);
        }

        sum /= (1.0f - 2.0f / (3.0f * k)) * M_PIf * radius * radius;
        return sum;
    }
}

Vec3 PhotonMap::computeIndirectIllumination(const RaySceneIntersection &intersection, const Settings &settings) const {
    Vec3 illumination(0.0f);

    // Add global photons if enabled
    if (settings.indirectIllumination) {
        constexpr float k = 1;
        const std::vector<Photon> nearbyGlobalPhotons = globalPhotonTree.findNearestNeighbors(
                intersection.intersection, intersection.normal, settings.maxIndirectDistance, settings.photonCountForIndirectColorEstimation
        );

        illumination += coneFilteredEstimate(nearbyGlobalPhotons, intersection.intersection,
                                             intersection.material.diffuse_material, k, settings.maxIndirectDistance);
    }

    // Add caustics photons if enabled
    if (settings.caustics) {
        constexpr float k = 1.5;
        const std::vector<Photon> nearbyCausticsPhotons = causticsPhotonTree.findNearestNeighbors(
            intersection.intersection, intersection.normal, settings.maxCausticsDistance, settings.photonCountForCausticsColorEstimation
        );

        illumination += coneFilteredEstimate(nearbyCausticsPhotons, intersection.intersection,
                                             intersection.material.diffuse_material, k, settings.maxCausticsDistance);
    }

    return illumination;
}


/******************************************************************************************************************/
/********************************************** UTILITY FUNCTIONS *************************************************/
/******************************************************************************************************************/

Photon PhotonMap::createInitialPhoton(const std::vector<Light>& lights, std::mt19937& rng) {
    Photon photon;
    std::uniform_int_distribution lightDist(0, static_cast<int>(lights.size()) - 1);
    const Light& light = lights[lightDist(rng)];

    photon.position = light.position;
    photon.color = light.material;
    photon.direction = randomDirection(rng);
    return photon;
}

Photon PhotonMap::createInitialPhotonTowardsObjects(const std::vector<Light>& lights, std::mt19937& rng,
    const std::vector<Sphere>& spheres, const std::vector<Mesh>& meshes) {
    Photon photon;
    std::uniform_int_distribution lightDist(0, static_cast<int>(lights.size()) - 1);
    const Light& light = lights[lightDist(rng)];

    photon.position = light.position;
    photon.color = light.material;
    photon.direction = randomDirectionTowardsObjects(rng, spheres, meshes);
    return photon;
}

void PhotonMap::computeReflectionProbabilities(const Material& material, float &Pr, float &Pd, float &Ps) {
    // Compute the probabilities of reflection, diffuse reflection, and specular reflection
    const float maxReflection = std::max({ // Pr = max(d + s)
        material.diffuse_material[0] + material.specular_material[0],
        material.diffuse_material[1] + material.specular_material[1],
        material.diffuse_material[2] + material.specular_material[2]
    });
    Pr = maxReflection;
    const float totalReflection = // dr + dg + db + sr + sg + sb
        material.diffuse_material[0] + material.diffuse_material[1] + material.diffuse_material[2] +
        material.specular_material[0] + material.specular_material[1] + material.specular_material[2];

    Pd = (material.diffuse_material[0] + material.diffuse_material[1] + material.diffuse_material[2]) / totalReflection * Pr;
    Ps = Pr - Pd;
}

Vec3 PhotonMap::randomDirection(std::mt19937 &rng, const Vec3 &normal) {
    std::uniform_real_distribution dist(0.0f, 1.0f);
    const float theta = 2.0f * M_PIf * dist(rng); // Azimuthal angle
    const float phi = std::acos(2.0f * dist(rng) - 1.0f); // Polar angle for full sphere
    const float x = std::sin(phi) * std::cos(theta);
    const float y = std::cos(phi);
    const float z = std::sin(phi) * std::sin(theta);

    Vec3 direction(x, y, z);

    // Ensure the direction is within the hemisphere defined by the normal
    if (normal != Vec3(0.0f) && Vec3::dot(direction, normal) < 0.0f) {
        direction = -direction;
    }

    return direction;
}

Vec3 PhotonMap::randomDirectionTowardsObjects(std::mt19937 &rng, const std::vector<Sphere>& spheres, const std::vector<Mesh>& meshes) {
    std::uniform_real_distribution dist(0.0f, 1.0f);

    // Select a random object type
    const int totalObjects = static_cast<int>(spheres.size() + meshes.size());
    std::uniform_int_distribution objectDist(0, totalObjects - 1);
    const int objectIndex = objectDist(rng);

    Vec3 target;

    if (objectIndex < static_cast<int>(spheres.size())) {
        // Select a random point on the sphere
        const Sphere& sphere = spheres[objectIndex];
        const float theta = dist(rng) * 2.0f * M_PIf;
        const float phi = std::acos(2.0f * dist(rng) - 1.0f);
        const float x = sphere.m_radius * std::sin(phi) * std::cos(theta);
        const float y = sphere.m_radius * std::sin(phi) * std::sin(theta);
        const float z = sphere.m_radius * std::cos(phi);
        target = sphere.m_center + Vec3(x, y, z);
    } else {
        // Select a random point on the mesh
        const Mesh& mesh = meshes[objectIndex - spheres.size()];
        target = mesh.getRandomPointOnSurface(rng).position;
    }

    return target;
}


/******************************************************************************************************************/
/******************************************* DEBUG AND VISUALIZATION **********************************************/
/******************************************************************************************************************/

void PhotonMap::debugDrawPhotons(const int type) {
    glDisable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);

    if (globalPhotons.empty()) globalPhotons = globalPhotonTree.toVector();
    if (causticsPhotons.empty()) causticsPhotons = causticsPhotonTree.toVector();

    /* --------------------------------- Draw Photons --------------------------------- */
    glBegin(GL_POINTS);

    if (type > 2) {
        for (const Photon &photon: globalPhotons) {
            glColor3f(photon.color[0], photon.color[1], photon.color[2]);
            glVertex3f(photon.position[0], photon.position[1], photon.position[2]);
        }
    }

    if (type == 1 || type == 2 || type > 4) {
        for (const Photon &photon: causticsPhotons) {
            glColor3f(photon.color[0], photon.color[1], photon.color[2]);
            glVertex3f(photon.position[0], photon.position[1], photon.position[2]);
        }
    }
    glEnd();

    /* --------------------------------- Draw Photon Paths --------------------------------- */
    if ((type == 4 || type == 6) && !initialGlobalPhotons.empty()) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glColor4f(initialGlobalPhotons[0].color[0], initialGlobalPhotons[0].color[1], initialGlobalPhotons[0].color[2], 0.1f);
        glBegin(GL_LINES);
        for (int i = 0; i < std::min(10000, static_cast<int>(initialGlobalPhotons.size())); i++) {
            const Photon &photon = initialGlobalPhotons[i];
            glVertex3f(photon.position[0], photon.position[1], photon.position[2]);
            glVertex3f(photon.direction[0] + photon.position[0], photon.direction[1] + photon.position[1], photon.direction[2] + photon.position[2]);
        }
        glEnd();

        glDisable(GL_BLEND);
    }

    if ((type == 2 || type == 6) && !initialCausticsPhotons.empty()) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glColor4f(initialCausticsPhotons[0].color[0], initialCausticsPhotons[0].color[1], initialCausticsPhotons[0].color[2], 0.1f);
        glBegin(GL_LINES);
        for (int i = 0; i < std::min(10000, static_cast<int>(initialCausticsPhotons.size())); i++) {
            const Photon &photon = initialCausticsPhotons[i];
            glVertex3f(photon.position[0], photon.position[1], photon.position[2]);
            glVertex3f(photon.direction[0], photon.direction[1], photon.direction[2]);
        }
        glEnd();

        glDisable(GL_BLEND);
    }

    glEnable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
}