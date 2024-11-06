#include "../include/PhotonMap.h"
#include "../include/Intersection.h"
#include "../include/Lighting.h"
#include <cmath>
#include <thread>
#include <GL/gl.h>

/******************************************************************************************************************/
/******************************************** PHOTON EMISSION INTERFACE *******************************************/
/******************************************************************************************************************/

void PhotonMap::emitPhotons(
    const std::vector<Light>& lights, const std::vector<Sphere>& spheres, const std::vector<Square>& squares,
    const std::vector<Mesh>& meshes, const MeshKDTree& kdTree, const Settings& settings)
{
    if (settings.globalPhotons <= 0) {
        throw std::invalid_argument("Number of photons must be greater than zero.");
    }
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
    std::vector<Photon> photonsToEmit;
    std::mutex photonMutex;

    const int numThreads = static_cast<int>(std::thread::hardware_concurrency());
    int photonsPerThread = photonCount / numThreads;

    auto emitTask = [this, &lights, &spheres, &squares, &meshes, &kdTree, &photonMutex, &photons, &photonsToEmit, photonsPerThread, photonType]() {
        emitPhotonsForThread(
            photonsPerThread, photonType, lights, spheres, squares, meshes, kdTree,
            photonMutex, photons, photonsToEmit);
    };

    std::vector<std::thread> threads(numThreads);
    for (auto& thread : threads) {
        thread = std::thread(emitTask);
    }

    for (auto& thread : threads) {
        thread.join();
    }

    return PhotonKDTree(photons);
}

/******************************************************************************************************************/
/******************************************* PER-THREAD PHOTON EMISSION *******************************************/
/******************************************************************************************************************/

void PhotonMap::emitPhotonsForThread(
    const int photonCount, const int photonType,
    const std::vector<Light>& lights, const std::vector<Sphere>& spheres, const std::vector<Square>& squares, const std::vector<Mesh>& meshes, const MeshKDTree& kdTree,
    std::mutex& photonMutex, std::vector<Photon>& photons, std::vector<Photon>& photonsToEmit)
{
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution dist(0.0f, 1.0f);

    std::vector<Photon> localPhotons;
    std::vector<Photon> localPhotonsToEmit;
    std::vector<Photon> localInitialPhotons;
    localPhotonsToEmit.reserve(photonCount * 2);

    int emittedPhotons = 0;
    while (localPhotons.size() < static_cast<size_t>(photonCount)) {
        Photon photon = photonType == 0 ? createInitialPhoton(lights, rng) : createInitialPhotonTowardsObjects(lights, rng, spheres, meshes);
        localPhotonsToEmit.emplace_back(photon);
        localInitialPhotons.emplace_back(photon);
        emittedPhotons++;
        processPhotonPath(photon, spheres, squares, meshes, kdTree, dist, rng, localPhotonsToEmit, localPhotons, photonType);
    }

    std::cout << localPhotons.size() << " photons emitted " << emittedPhotons << std::endl;
    normalizePhotonColors(localPhotons, emittedPhotons);

    std::lock_guard lock(photonMutex);
    photons.insert(photons.end(), localPhotons.begin(), localPhotons.end());
    photonsToEmit.insert(photonsToEmit.end(), localPhotonsToEmit.begin(), localPhotonsToEmit.end());
    if (photonType == 0) {
        initialGlobalPhotons.insert(initialGlobalPhotons.end(), localInitialPhotons.begin(), localInitialPhotons.end());
    } else {
        initialCausticsPhotons.insert(initialCausticsPhotons.end(), localInitialPhotons.begin(), localInitialPhotons.end());
    }
}

void PhotonMap::processPhotonPath(
    Photon& photon, const std::vector<Sphere>& spheres, const std::vector<Square>& squares, const std::vector<Mesh>& meshes, const MeshKDTree& kdTree,
    std::uniform_real_distribution<float>& dist, std::mt19937& rng, std::vector<Photon>& localPhotonsToEmit, std::vector<Photon>& localPhotons, const int photonType)
{
    bool absorbed = false;
    int bounces = 0;
    constexpr int maxBounces = 500;

    while (!absorbed && bounces < maxBounces) {
        Ray ray(photon.position, photon.direction);
        RaySceneIntersection intersection = Intersection::computeIntersection(ray, spheres, squares, meshes, kdTree, 0.0f);
        if (!intersection.intersectionExists) break;

        photon.position = intersection.intersection;
        processPhotonInteraction(photon, intersection, dist, rng, localPhotonsToEmit, localPhotons, photonType, absorbed);
        bounces++;
    }
}

void PhotonMap::processPhotonInteraction(
    Photon& photon, const RaySceneIntersection& intersection, std::uniform_real_distribution<float>& dist, std::mt19937& rng,
    std::vector<Photon>& localPhotonsToEmit, std::vector<Photon>& localPhotons, const int photonType, bool& absorbed)
{
    switch (intersection.material.type) {
        case Material_Glass:
            processGlassPhoton(photon, intersection.normal, intersection.material, localPhotonsToEmit);
            break;
        case Material_Mirror:
            processMirrorPhoton(photon, intersection.normal, intersection.material);
            break;
        case Material_Diffuse_Blinn_Phong:
            processDiffusePhoton(photon, intersection, dist, rng, localPhotons, photonType, absorbed);
            break;
    }
}

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
    reflectedPhoton.flag = 1;
    photons.emplace_back(reflectedPhoton);

    // Add object color to the photon color (with transparency)
    photon.direction = refractedDirection;
    photon.flag = 1;
    photon.color = (1 - fresnelEffect) * (material.transparency * photon.color + (1 - material.transparency) * material.specular_material);
}

void PhotonMap::processMirrorPhoton(Photon& photon, const Vec3& normal, const Material& material) {
    // Change the direction of the photon to make it bounce, keep the color at 90% to simulate reflection
    photon.direction = Lighting::computeReflectedDirection(photon.direction, normal).normalize();
    photon.color = 0.9f * Vec3::compProduct(photon.color, material.specular_material);
    photon.flag = 1;
}

void PhotonMap::processDiffusePhoton(
    Photon& photon, const RaySceneIntersection& intersection, std::uniform_real_distribution<float>& dist, std::mt19937& rng,
    std::vector<Photon>& localPhotons, const int photonType, bool& absorbed)
{
    float Pr, Pd, Ps;
    computeReflectionProbabilities(intersection.material, Pr, Pd, Ps);
    const float xi = dist(rng);

    if (xi < Pd) {  // Diffuse reflection ξ ∈[0, Pd]
        photon.direction = randomDirection(rng, intersection.normal);
        photon.color = Vec3::compProduct(photon.color, intersection.material.diffuse_material) / Pd;
        if (photonType == 1) absorbed = true;
    } else if (xi < Pd + Ps) {  // Specular reflection ξ ∈]Pd, Ps + Pd]
        photon.direction = Lighting::computeReflectedDirection(photon.direction, intersection.normal).normalize();
        photon.color = Vec3::compProduct(photon.color, intersection.material.specular_material) / Pr;
        photon.flag = 1;
    } else {  // Absorption ξ ∈]Ps + Pd, 1]
        absorbed = true;
    }

    if (photonType == 0 || (photonType == 1 && photon.flag == 1)) {
        localPhotons.emplace_back(photon);
    }
}

/******************************************************************************************************************/
/******************************************** PHOTON PROCESSES PER MATERIAL ***************************************/
/******************************************************************************************************************/

void PhotonMap::normalizePhotonColors(std::vector<Photon> &photons, const int emittedPhotons) {
    for (Photon& photon : photons) {
        photon.debugColor = photon.color;
        photon.color /= static_cast<float>(emittedPhotons);
    }
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
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
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
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    // Select a random object type
    const int totalObjects = static_cast<int>(spheres.size() + meshes.size());
    std::uniform_int_distribution<int> objectDist(0, totalObjects - 1);
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

    // Generate a direction vector from the origin to the target point
    Vec3 direction = target - Vec3(0.0f, 0.0f, 0.0f);
    return direction;
}

/******************************************************************************************************************/
/******************************************** RENDERING PHOTON PATHS **********************************************/
/******************************************************************************************************************/
Vec3 PhotonMap::computeCaustics(const Vec3 &position, const Material &material) const {
    // Find nearby photons in the global photon map
    constexpr float globalDistance = 0.5f;
    const std::vector<Photon> nearbyGlobalPhotons = globalPhotonTree.findNearestNeighbors(position, globalDistance);

    Vec3 illumination(0.0f);
    if (nearbyGlobalPhotons.empty()) return illumination;

    // Accumulate the color of the nearby photons
    for (const Photon &photon : nearbyGlobalPhotons) {
        illumination += Vec3::compProduct(photon.color, material.diffuse_material);
    }

    return illumination;
}


/******************************************************************************************************************/
/******************************************* DEBUG AND VISUALIZATION **********************************************/
/******************************************************************************************************************/

void PhotonMap::debugDrawPhotons(const int type) {
    glDisable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);

    if (globalPhotons.empty()) globalPhotons = globalPhotonTree.toVector();
    if (causticsPhotons.empty()) causticsPhotons = causticsPhotonTree.toVector();

    if (type > 2) {
        for (const Photon &photon: globalPhotons) {
            glColor3f(photon.debugColor[0], photon.debugColor[1], photon.debugColor[2]);
            glBegin(GL_POINTS);
            glVertex3f(photon.position[0], photon.position[1], photon.position[2]);
            glEnd();
        }
    }

    if (type < 3 || type > 4) {
        for (const Photon &photon: causticsPhotons) {
            glColor3f(photon.debugColor[0], photon.debugColor[1], photon.debugColor[2]);
            glBegin(GL_POINTS);
            glVertex3f(photon.position[0], photon.position[1], photon.position[2]);
            glEnd();
        }
    }

    if (type == 4 || type == 6) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        for (const Photon &photon : initialGlobalPhotons) {
            glColor4f(initialGlobalPhotons[0].color[0], initialGlobalPhotons[0].color[1], initialGlobalPhotons[0].color[2], 0.2f);
            glBegin(GL_LINES);
            glVertex3f(photon.position[0], photon.position[1], photon.position[2]);
            glVertex3f(photon.direction[0] + photon.position[0], photon.direction[1] + photon.position[1], photon.direction[2] + photon.position[2]);
            glEnd();
        }

        glDisable(GL_BLEND);
    }

    if (type == 2 || type == 6) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        for (const Photon &photon : initialCausticsPhotons) {
            glColor4f(initialCausticsPhotons[0].color[0], initialCausticsPhotons[0].color[1], initialCausticsPhotons[0].color[2], 0.2f);
            glBegin(GL_LINES);
            glVertex3f(photon.position[0], photon.position[1], photon.position[2]);
            glVertex3f(photon.direction[0], photon.direction[1], photon.direction[2]);
            glEnd();
        }

        glDisable(GL_BLEND);
    }


    glEnable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
}
