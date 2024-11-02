#include "PhotonMap.h"
#include "Intersection.h"
#include "Lighting.h"
#include <cmath>
#include <utility>
#include <GL/gl.h>

/******************************************************************************************************************/
/******************************************** PHOTON EMISSION FUNCTIONS *******************************************/
/******************************************************************************************************************/

void PhotonMap::emitPhotons(
    const std::vector<Light> &lights, const std::vector<Sphere> &spheres,
    const std::vector<Square> &squares, const std::vector<Mesh> &meshes, const int photons)
{
    std::mt19937 rng(std::random_device{}());
    if (lights.empty()) {
        std::cerr << "Error: No lights in the scene, cannot emit photons." << std::endl;
        return;
    }

    std::vector<Photon> mirrorPhotons;
    std::vector<Photon> glassPhotons;

    // Create a vector of photons
    std::vector<Photon> photonsToEmit;
    photonsToEmit.reserve(photons * 2); // In the worst case, we double the number of photons (reflection + refraction)
    for (int i = 0; i < photons; ++i) {
        Photon photon = createInitialPhoton(lights, rng);
        photonsToEmit.emplace_back(photon);
    }

    for (auto &photon : photonsToEmit) {
        bool absorbed = false;

        for (int j = 0; j < 10 && !absorbed; j++) {
            Ray ray(photon.position, photon.direction);
            RaySceneIntersection intersection = Intersection::computeIntersection(ray, spheres, squares, meshes, 0.0f);

            if (!intersection.intersectionExists) break;
            auto [intersectionPoint, normal, material] = Intersection::parseIntersection(intersection, spheres, squares, meshes);

            // Calculate the distance traveled by the photon and apply attenuation based on it
            const float distanceTraveled = (intersectionPoint - photon.position).length();
            photon.position = intersectionPoint;
            float attenuationFactor = 1.0f;

            if (material.type == Material_Diffuse_Blinn_Phong || material.type == Material_Mirror) { // Distance traveled in air
                constexpr float absorptionCoefficientAir = 0.02f; // Absorption coefficient for air in m^-1
                attenuationFactor = std::exp(-distanceTraveled * absorptionCoefficientAir);
            } else if (material.type == Material_Glass) { // Distance traveled in glass
                constexpr float absorptionCoefficientGlass = 0.1f; // Absorption coefficient for glass in m^-1
                attenuationFactor = std::exp(-distanceTraveled * absorptionCoefficientGlass);
            }

            photon.color *= attenuationFactor;

            // We do not process the colors of the photons when they are absorbed, it is handled during rendering
            switch (material.type) {
                case Material_Diffuse_Blinn_Phong:
                    if (photon.materialType == Material_Glass) {
                        glassPhotons.emplace_back(photon);
                        absorbed = true;
                    } else if (photon.materialType == Material_Mirror) {
                        mirrorPhotons.emplace_back(photon);
                        absorbed = true;
                    }
                break;
                case Material_Glass:
                    processGlassPhoton(photon, normal, material, photonsToEmit);
                    break;
                case Material_Mirror:
                    processMirrorPhoton(photon, normal, material);
                    break;
            }
        }
    }

    mirrorPhotonTree = PhotonKDTree(mirrorPhotons);
    glassPhotonTree = PhotonKDTree(glassPhotons);
}

Photon PhotonMap::createInitialPhoton(const std::vector<Light>& lights, std::mt19937& rng) {
    Photon photon;
    std::uniform_int_distribution<int> lightDist(0, static_cast<int>(lights.size()) - 1);
    const Light& light = lights[lightDist(rng)];

    photon.position = light.position;
    photon.color = light.material;
    photon.direction = randomDirection(rng);
    photon.materialType = Material_Diffuse_Blinn_Phong;
    return photon;
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
    reflectedPhoton.materialType = material.type;
    photons.emplace_back(reflectedPhoton);

    // Add object color to the photon color (with transparency)
    photon.direction = refractedDirection;
    photon.color = (1 - fresnelEffect) * (material.transparency * photon.color + (1 - material.transparency) * material.specular_material);
    photon.materialType = material.type;
}

void PhotonMap::processMirrorPhoton(Photon& photon, const Vec3& normal, const Material& material) {
    // Change the direction of the photon to make it bounce, keep the color at 90% to simulate reflection
    photon.direction = Lighting::computeReflectedDirection(photon.direction, normal).normalize();
    photon.color = 0.9f * Vec3::compProduct(photon.color, material.specular_material);
    photon.materialType = material.type;
}


/******************************************************************************************************************/
/******************************************* PHOTON RENDERING FUNCTIONS *******************************************/
/******************************************************************************************************************/

Vec3 PhotonMap::computeCaustics(const Vec3 &position, const Material &material) const {
    constexpr float maxMirrorDistance = 1.5f;
    constexpr float maxGlassDistance = 0.2f;

    std::vector<Photon> nearbyPhotons = mirrorPhotonTree.findNearestNeighbors(position, maxMirrorDistance);
    std::vector<Photon> nearbyGlassPhotons = glassPhotonTree.findNearestNeighbors(position, maxGlassDistance);

    Vec3 illumination(0.0f);
    if (nearbyPhotons.empty() && nearbyGlassPhotons.empty()) return illumination;

    for (const Photon &photon : nearbyPhotons) {
        illumination += Vec3::compProduct(photon.color, material.diffuse_material);
    }

    for (const Photon &photon : nearbyGlassPhotons) {
        illumination += Vec3::compProduct(photon.color, material.diffuse_material);
    }

    return illumination;
}

Vec3 PhotonMap::randomDirection(std::mt19937 &rng) {
    std::uniform_real_distribution dist(0.0f, 1.0f);
    const float theta = 2.0f * M_PIf * dist(rng);
    const float phi = std::acos(1.0f - 2.0f * dist(rng)); // Uniform distribution over the hemisphere

    float x = std::sin(phi) * std::cos(theta);
    float y = -std::cos(phi); // Negative to limit the direction to the upper hemisphere
    float z = std::sin(phi) * std::sin(theta);

    return {x, y, z};
}


/******************************************************************************************************************/
/********************************************** DEBUG DRAW FUNCTIONS **********************************************/
/******************************************************************************************************************/

void PhotonMap::debugDrawPhotons() const {
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

    // Debug: Draw the photons
    std::vector<Photon> photons = mirrorPhotonTree.toVector();
    for (const auto &[position, direction, color, materialType] : photons) {
        glColor3f(0.0f, 0.4f, 0.4f);
        glBegin(GL_POINTS);
        glVertex3f(position[0], position[1], position[2]);
        glEnd();
    }

    std::vector<Photon> glassPhotons = glassPhotonTree.toVector();
    for (const auto &[position, direction, color, materialType] : glassPhotons) {
        glColor3f(0.0f, 0.3f, 0.8f);
        glBegin(GL_POINTS);
        glVertex3f(position[0], position[1], position[2]);
        glEnd();
    }

    glDisable(GL_COLOR_MATERIAL);
}
