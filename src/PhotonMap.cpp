#include "PhotonMap.h"

void PhotonMap::emitPhotons(const std::vector<Light> &lights, const std::vector<Sphere> &spheres, const std::vector<Square> &squares, const std::vector<Mesh> &meshes, int photons) {
    std::mt19937 rng(std::random_device{}());
        std::vector<Photon> photonsGlass, photonsMirror, photonsDiffuse;

        for (int i = 0; i < photons; ++i) {
            if (lights.empty()) {
                std::cerr << "No lights in the scene." << std::endl;
                return;
            }

            Photon photon;
            photon.position = lights[0].pos; // Photon starts at the light position
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

Vec3 PhotonMap::renderCaustics(const Vec3 &position, const Vec3 &normal, const Material &material) const {
    Vec3 causticsColor(0.0f, 0.0f, 0.0f);

    // Find nearby photons for the glass material and add them to the caustics color
    std::vector<Photon> nearbyGlassPhotons = photonGlassTree.findNearestNeighbors(position, 0.4f);
    float totalWeightGlass = 0.0f;
    for (auto &[photonPosition, direction, color] : nearbyGlassPhotons) {
        const float cosTheta = Vec3::dot(normal, direction);
        if (cosTheta > 0) {
            const float distance = (photonPosition - position).length();
            const float weight = 1.0f / (distance * distance + 1e-4f); // Pondération basée sur la distance
            causticsColor += Vec3::compProduct(color, material.diffuse_material) * cosTheta * weight;
            totalWeightGlass += weight;
        }
    }
    if (totalWeightGlass > 0) {
        causticsColor /= totalWeightGlass; // Normalisation
    }

    // Find nearby photons for the mirror material and add them to the caustics color
    std::vector<Photon> nearbyMirrorPhotons = photonMirrorTree.findNearestNeighbors(position, 1.5f);
    for (const auto &[position, direction, color] : nearbyMirrorPhotons) {
        const float cosTheta = Vec3::dot(normal, direction);
        if (cosTheta > 0) {
            causticsColor += Vec3::compProduct(color, material.diffuse_material) * cosTheta / static_cast<float>(nearbyMirrorPhotons.size());
        }
    }

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
    for (const auto &[position, direction, color] : glassPhotons) {
        glColor3f(0.0f, 0.3f, 0.8f);  // Dark blue for Glass Photons
        glBegin(GL_LINES);
        glVertex3f(position[0], position[1], position[2]);
        glVertex3f(position[0] + direction[0], position[1] + direction[1], position[2] + direction[2]);
        glEnd();
    }

    std::vector<Photon> MirrorPhotons = photonMirrorTree.toVector();
    for (const auto &[position, direction, color] : MirrorPhotons) {
        glColor3f(0.3f, 0.3f, 0.3f);  // Dark metallic gray for Mirror Photons
        glBegin(GL_LINES);
        glVertex3f(position[0], position[1], position[2]);
        glVertex3f(position[0] + direction[0], position[1] + direction[1], position[2] + direction[2]);
        glEnd();
    }

    std::vector<Photon> DiffusePhotons = photonDiffuseTree.toVector();
    for (const auto &[position, direction, color] : DiffusePhotons) {
        glColor3f(0.0f, 0.5f, 0.5f);  // Saturated Cyan for Diffuse Photons
        glBegin(GL_LINES);
        glVertex3f(position[0], position[1], position[2]);
        glVertex3f(position[0] + direction[0], position[1] + direction[1], position[2] + direction[2]);
        glEnd();
    }

    glDisable(GL_COLOR_MATERIAL);
}