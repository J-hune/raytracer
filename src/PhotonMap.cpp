#include "PhotonMap.h"
#include "Intersection.h"
#include "Lighting.h"
#include <utility>
#include <cmath>
#include <GL/gl.h>

void PhotonMap::emitPhotons(const std::vector<Light> &lights, const std::vector<Sphere> &spheres, const std::vector<Square> &squares, const std::vector<Mesh> &meshes, int photons) {
    std::mt19937 rng(std::random_device{}());
    std::vector<Photon> photonsGlass, photonsMirror, photonsDiffuse;

    for (int i = 0; i < photons; ++i) {
        if (lights.empty()) {
            std::cerr << "No lights in the scene." << std::endl;
            return;
        }

        Photon photon;
        photon.position = lights[0].position; // Photon starts at the light position
        photon.color = lights[0].material; // Photon color is the light color
        photon.direction = randomDirection(rng); // Random direction for the photon

        for (int j = 0; j < 5; j++) {
            Ray ray(photon.position, photon.direction);
            RaySceneIntersection intersection = Intersection::computeIntersection(ray, spheres, squares, meshes, 0.0f);

            if (intersection.intersectionExists) {
                auto [intersectionPoint, normal, material] = Intersection::parseIntersection(intersection, spheres, squares, meshes);

                // Glass material: refraction with transparency
                if (material.type == Material_Glass) {
                    const Vec3 refractedDirection = Lighting::computeRefractedDirection(photon.direction, normal, material.index_medium);
                    photon.position = intersectionPoint;
                    photon.direction = refractedDirection;

                    // Add object color to the photon color (with transparency)
                    photon.color = material.transparency * photon.color + (1 - material.transparency) * material.diffuse_material;
                    photonsGlass.emplace_back(photon);
                }

                // Mirror material: reflection
                else if (material.type == Material_Mirror) {
                    photon.position = intersectionPoint + normal * 1e-4f; // small bias to avoid self-intersection
                    photon.direction = Lighting::computeReflectedDirection(photon.direction, normal).normalize();
                    photon.color = Vec3::compProduct(photon.color, material.diffuse_material);
                    photonsMirror.emplace_back(photon);
                }

                // Diffuse material
                else {
                    photon.color = Vec3::compProduct(photon.color, material.diffuse_material);
                    photonsDiffuse.emplace_back(photon);
                    break; // Absorption for non-dielectric materials
                }
            } else {
                break; // No intersection, stop the photon
            }
        }
    }

    photonGlassTree = PhotonKDTree(photonsGlass);
    photonMirrorTree = PhotonKDTree(photonsMirror);
    photonDiffuseTree = PhotonKDTree(photonsDiffuse);
}

/**
 * This function, renderCaustics, calculates the caustics color at a specific position on a surface based on nearby photons,
 * simulating realistic lighting effects such as concentrated light patterns (caustics) generated by materials like glass.
 *
 * The caustics effect is calculated by:
 * 1. Searching for nearby photons within a specified radius to determine the relevant photons contributing to the caustic effect.
 * 2. Applying a Gaussian weight to each photon, favoring those closer to the target position.
 * 3. Normalizing the accumulated photon contributions based on their density, using logarithmic compression to smooth the intensity of densely packed photons.
 * 4. Clamping the final color values to avoid overly bright or unnatural artifacts.
 *
 * This method is physically consistent as it considers photon density, light angle, and distance, leading to a more realistic caustics effect.
 **/
Vec3 PhotonMap::renderCaustics(const Vec3 &position, const Vec3 &normal, const Material &material) const {
    constexpr float maxMirrorRadius = 1.5f; // Search radius for nearby mirror photons
    constexpr float sigmaMirror = 0.2f; // Gaussian weight for controlling photon influence by distance

    constexpr float maxGlassRadius = 1.0f; // Search radius for finding relevant photons
    constexpr float sigmaGlass = 0.12f;

    // Calculate the caustics color based on mirror and glass photon contributions
    Vec3 causticsColor(0.0f, 0.0f, 0.0f); //
    std::vector<Photon> nearbyMirrorPhotons = findNearestMirrorPhotons(position, maxMirrorRadius);
    std::vector<Photon> nearbyGlassPhotons = findNearestGlassPhotons(position, maxGlassRadius);

    causticsColor += calculateCausticsContribution(nearbyMirrorPhotons, position, normal, material, maxMirrorRadius, sigmaMirror);
    causticsColor += calculateCausticsContribution(nearbyGlassPhotons, position, normal, material, maxGlassRadius, sigmaGlass);
    return causticsColor;
}

Vec3 PhotonMap::calculateCausticsContribution(
    std::vector<Photon> &photons, const Vec3 &position, const Vec3 &normal,
    const Material &material, const float maxRadius, const float sigma)
{
    Vec3 causticsColor(0.0f, 0.0f, 0.0f); // Initialize the caustics color
    constexpr float minWeightThreshold = std::numeric_limits<float>::min();
    // Minimum threshold to prevent near-zero division artifacts

    // Only calculate if there are any photons found
    if (!photons.empty()) {
        float totalWeight = 0.0f; // Accumulates the total weight for normalization

        // Process each photon in the nearby list
        for (const auto &[photonPosition, direction, color]: photons) {
            const float cosTheta = Vec3::dot(normal, direction);
            // Dot product between surface normal and photon direction

            // Only consider photons that hit the surface at a positive angle
            if (cosTheta > 0) {
                const float distance = (photonPosition - position).length(); // Distance from photon to target position
                const float weight = std::exp(-distance * distance / (2.0f * sigma * sigma));
                // Gaussian weight based on distance

                // Accumulate the weighted contribution of this photon to the caustics color
                causticsColor += Vec3::compProduct(color, material.diffuse_material) * cosTheta * weight;
                totalWeight += weight; // Accumulate total weight for normalization
            }
        }

        // Normalize the color intensity by the density of photons in the area, applying logarithmic smoothing
        if (totalWeight > minWeightThreshold) {
            float densityFactor = std::sqrt(static_cast<float>(photons.size()) / (M_PIf * maxRadius * maxRadius));
            densityFactor = std::log(1.0f + densityFactor); // Logarithmic compression of density for smoother intensity

            // Apply normalized density factor to the accumulated color
            causticsColor *= densityFactor / totalWeight;
        } else {
            causticsColor = Vec3(0.0f); // If weight is too low, set caustics color to zero (no visible effect)
        }
    }

    // Clamp the color to avoid overly bright artifacts in the final result
    causticsColor = causticsColor.min(Vec3(1.0f, 1.0f, 1.0f));

    return causticsColor;
}

std::vector<Photon> PhotonMap::findNearestGlassPhotons(const Vec3 &point, const float maxDistance) const {
    return photonGlassTree.findNearestNeighbors(point, maxDistance);
}

std::vector<Photon> PhotonMap::findNearestMirrorPhotons(const Vec3 &point, const float maxDistance) const {
    return photonMirrorTree.findNearestNeighbors(point, maxDistance);
}

std::vector<Photon> PhotonMap::findNearestDiffusePhotons(const Vec3 &point, const float maxDistance) const {
    return photonDiffuseTree.findNearestNeighbors(point, maxDistance);
}

Vec3 PhotonMap::randomDirection(std::mt19937 &rng) {
    std::uniform_real_distribution dist(0.0f, 1.0f);
    const float theta = 2.0f * M_PIf * dist(rng);
    constexpr float coneAngle = M_PIf / 6.f * 2.5f; // Cone angle to limit the spread of the photons
    const float phi = coneAngle * dist(rng); // Limit phi between 0 and coneAngle

    float x = std::sin(phi) * std::cos(theta);
    float y = -std::cos(phi); // Negative to limit the direction to the upper hemisphere
    float z = std::sin(phi) * std::sin(theta);

    return {x, y, z};
}

void PhotonMap::debugDrawPhotons() const {
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

    // Debug : Draw the photons
    std::vector<Photon> glassPhotons = photonGlassTree.toVector();
    for (const auto &[position, direction, color]: glassPhotons) {
        glColor3f(0.0f, 0.3f, 0.8f); // Dark blue for Glass Photons
        glBegin(GL_LINES);
        glVertex3f(position[0], position[1], position[2]);
        glVertex3f(position[0] + direction[0], position[1] + direction[1], position[2] + direction[2]);
        glEnd();
    }

    std::vector<Photon> MirrorPhotons = photonMirrorTree.toVector();
    for (const auto &[position, direction, color]: MirrorPhotons) {
        glColor3f(0.3f, 0.3f, 0.3f); // Dark metallic gray for Mirror Photons
        glBegin(GL_LINES);
        glVertex3f(position[0], position[1], position[2]);
        glVertex3f(position[0] + direction[0], position[1] + direction[1], position[2] + direction[2]);
        glEnd();
    }

    std::vector<Photon> DiffusePhotons = photonDiffuseTree.toVector();
    for (const auto &[position, direction, color]: DiffusePhotons) {
        glColor3f(0.0f, 0.5f, 0.5f); // Saturated Cyan for Diffuse Photons
        glBegin(GL_LINES);
        glVertex3f(position[0], position[1], position[2]);
        glVertex3f(position[0] + direction[0], position[1] + direction[1], position[2] + direction[2]);
        glEnd();
    }

    glDisable(GL_COLOR_MATERIAL);
}
