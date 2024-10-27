#ifndef SCENE_H
#define SCENE_H

#include <cmath>
#include <map>
#include <vector>
#include <string>
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"
#include "Light.h"
#include "PhotonKDTree.h"

struct RaySceneIntersection {
    bool intersectionExists;
    unsigned int typeOfIntersectedObject{};
    unsigned int objectIndex{};
    float t;
    RayTriangleIntersection rayMeshIntersection;
    RaySphereIntersection raySphereIntersection;
    RaySquareIntersection raySquareIntersection;

    RaySceneIntersection() : intersectionExists(false), t(FLT_MAX), rayMeshIntersection(), raySquareIntersection() {
    }
};

class Scene {
    std::vector<Mesh> meshes;
    std::vector<Sphere> spheres;
    std::vector<Square> squares;
    std::vector<Light> lights;
    PhotonKDTree photonGlassTree;
    PhotonKDTree photonMirrorTree;
    PhotonKDTree photonDiffuseTree;

public:
    Scene() = default;

    void draw() const {
        // Iterate over all objects and render them:
        for (const auto &mesh: meshes) mesh.draw();
        for (const auto &sphere: spheres) sphere.draw();
        for (const auto &square: squares) square.draw();

        // Debug : Draw the photons
        std::vector<Photon> glassPhotons = photonGlassTree.toVector();
        for (const auto &[position, direction, color] : glassPhotons) {
            glColor3f(0, 0, 1); // Couleur bleue
            glBegin(GL_LINES);
            glVertex3f(position[0], position[1], position[2]);
            glVertex3f(position[0] + direction[0], position[1] + direction[1], position[2] + direction[2]);
            glEnd();
        }

        std::vector<Photon> MirrorPhotons = photonMirrorTree.toVector();
        for (const auto &[position, direction, color] : MirrorPhotons) {
            glColor3f(0, 1, 0); // Couleur verte
            glBegin(GL_LINES);
            glVertex3f(position[0], position[1], position[2]);
            glVertex3f(position[0] + direction[0], position[1] + direction[1], position[2] + direction[2]);
            glEnd();
        }

        std::vector<Photon> DiffusePhotons = photonDiffuseTree.toVector();
        for (const auto &[position, direction, color] : DiffusePhotons) {
            glColor3f(1, 1, 1); // Couleur blanc
            glBegin(GL_LINES);
            glVertex3f(position[0], position[1], position[2]);
            glVertex3f(position[0] + direction[0], position[1] + direction[1], position[2] + direction[2]);
            glEnd();
        }
    }

    RaySceneIntersection computeIntersection(const Ray &ray, const float z_near) {
        RaySceneIntersection result;
        result.t = FLT_MAX;

        // For each object in the scene, compute the intersection with the ray
        for (unsigned int i = 0; i < spheres.size(); i++) {
            const RaySphereIntersection intersection = spheres[i].intersect(ray);

            // If the intersection exists and is closer than the previous one, keep it
            if (intersection.intersectionExists && intersection.t <= result.t) {
                result.intersectionExists = true;
                result.raySphereIntersection = intersection;
                result.t = intersection.t;
                result.objectIndex = i;
                result.typeOfIntersectedObject = 0;
            }
        }

        for (unsigned int i = 0; i < squares.size(); i++) {
            const RaySquareIntersection intersection = squares[i].intersect(ray);

            // If the intersection exists and is closer than the previous one, keep it
            if (intersection.intersectionExists && intersection.t <= result.t && intersection.t > z_near) {
                result.intersectionExists = true;
                result.raySquareIntersection = intersection;
                result.t = intersection.t;
                result.objectIndex = i;
                result.typeOfIntersectedObject = 1;
            }
        }

        for (unsigned int i = 0; i < meshes.size(); i++) {
            const RayTriangleIntersection intersection = meshes[i].intersect(ray);

            // If the intersection exists and is closer than the previous one, keep it
            if (intersection.intersectionExists && intersection.t <= result.t && intersection.t > z_near) {
                result.intersectionExists = true;
                result.rayMeshIntersection = intersection;
                result.t = intersection.t;
                result.objectIndex = i;
                result.typeOfIntersectedObject = 2;
            }
        }

        return result;
    }

    // Compute the reflected direction using the normal
    static Vec3 computeReflectedDirection(const Vec3 &direction, const Vec3 &normal) {
        return direction - 2.0f * Vec3::dot(direction, normal) * normal;
    }

    // Compute the refracted direction using Snell's Law
    static Vec3 computeRefractedDirection(const Vec3 &incident, const Vec3 &normal, const float &ior) {
        // Snell's Law: etaI * sin(thetaI) = etaT * sin(thetaT)
        float cosI = std::clamp(Vec3::dot(incident, normal), -1.0f, 1.0f);
        float etaI = 1.0f, etaT = ior;
        Vec3 n = normal;

        // If the ray is entering the material, invert the normal and swap the refractive indices
        if (cosI < 0) { cosI = -cosI; }
        else { std::swap(etaI, etaT); n = -n; }
        const float etaRatio = etaI / etaT;
        const float k = 1 - etaRatio * etaRatio * (1 - cosI * cosI);

        // If k is negative, there is total internal reflection, return a zero vector
        return k < 0 ? Vec3(0, 0, 0) : etaRatio * incident + (etaRatio * cosI - sqrtf(k)) * n;
    }

    // Compute the Fresnel effect using Schlick's approximation
    static float computeFresnelEffect(const Vec3 &I, const Vec3 &N, const float &ior) {
        float cosI = std::clamp(-1.0f, 1.0f, Vec3::dot(I, N));
        float etaI = 1.0f, etaT = ior;
        if (cosI > 0) { std::swap(etaI, etaT); }

        // Compute the sine of the angle using Snell's Law
        const float sinT = etaI / etaT * sqrtf(std::max(0.f, 1 - cosI * cosI));

        // Total internal reflection
        if (sinT >= 1) return 1.0f;

        const float cosT = sqrtf(std::max(0.f, 1 - sinT * sinT));
        cosI = fabsf(cosI);
        const float Rs = ((etaT * cosI) - (etaI * cosT)) / ((etaT * cosI) + (etaI * cosT));
        const float Rp = ((etaI * cosI) - (etaT * cosT)) / ((etaI * cosI) + (etaT * cosT));
        return (Rs * Rs + Rp * Rp) / 2;
    }

    static Vec3 samplePointOnQuad(const Light &light, std::mt19937 &rng) {
        std::uniform_real_distribution<float> dist(0.0f, 1.0f);
        const float u = dist(rng);
        const float v = dist(rng);

        return light.quad.vertices[0].position * (1 - u) * (1 - v) +
               light.quad.vertices[1].position * u * (1 - v) +
               light.quad.vertices[2].position * u * v +
               light.quad.vertices[3].position * (1 - u) * v;
    }

    // Compute the Phong components for the light (diffuse and specular only)
    static Vec3 computePhongComponents(const Vec3 &lightDir, const Vec3 &viewDir, const Vec3 &normal, const Material &material, const Light &light) {
        // Compute the diffuse component
        const float diff = std::max(Vec3::dot(normal, lightDir), 0.0f);
        const Vec3 diffuse = Vec3::compProduct(material.diffuse_material, light.material) * diff;

        // Compute the reflection direction and the specular component
        const Vec3 reflectDir = computeReflectedDirection(-lightDir, normal);
        const float spec = static_cast<float>(std::pow(std::max(Vec3::dot(viewDir, reflectDir), 0.0f), material.shininess));
        const Vec3 specular = Vec3::compProduct(material.specular_material, light.material) * spec;

        return diffuse + specular; // Combine diffuse and specular
    }

    // Recursive ray tracing function
    Vec3 rayTraceRecursive(const Ray &ray, const int NRemainingBounces, std::mt19937 &rng, const float z_near = 0.0f) {
        Vec3 color(0.0f, 0.0f, 0.0f);
        const Vec3 ambientLight(0.1f, 0.1f, 0.1f); // Ambient light

        // Compute the intersection with the scene (ray vs. (spheres, squares, meshes))
        const RaySceneIntersection intersection = computeIntersection(ray, z_near);

        // If there isn't any intersection, return the background color
        if (!intersection.intersectionExists || NRemainingBounces == 0) {
            return color;
        }

        auto [intersectionPoint, normal, material] = handleIntersection(intersection);

        // Add the ambient light to the color
        color += Vec3::compProduct(material.ambient_material, ambientLight);

        // Pre-normalize the view vector
        Vec3 viewDir = -ray.direction();
        viewDir.normalize();

        // For each light in the scene
        for (const auto &light: lights) {
            constexpr float epsilon = 1e-5f;
            Vec3 lightDir = (light.pos - intersectionPoint).normalize();;
            float lightDistance = lightDir.length();

            // If the light is spherical (point light) => hard shadows
            if (light.type == LightType_Spherical) {
                // Ray to test the shadow
                Ray shadowRay(intersectionPoint + normal * epsilon, lightDir);
                RaySceneIntersection shadowIntersection = computeIntersection(shadowRay, epsilon);

                // If the point is in shadow, skip the light
                if (shadowIntersection.intersectionExists && shadowIntersection.t < lightDistance) {
                    continue; // The point is in shadow.
                }

                color += computePhongComponents(lightDir, viewDir, normal, material, light);
            }

            // If the light is a quad (area light) => soft shadows
            else if (light.type == LightType_Quad) {
                constexpr int numSamples = 16; // Number of samples for soft shadows
                float shadowFactor = 0.0f;
                constexpr float threshold = numSamples * 0.9f;

                // Sampling for soft shadows
                for (int i = 0; i < numSamples; ++i) {
                    Vec3 samplePoint = samplePointOnQuad(light, rng);
                    Vec3 shadowDir = (samplePoint - intersectionPoint).normalize();
                    float shadowDistance = shadowDir.length();

                    // Ray to test the shadow
                    Ray shadowRay(intersectionPoint + normal * epsilon, shadowDir);
                    RaySceneIntersection shadowIntersection = computeIntersection(shadowRay, epsilon);

                    // If the point is in shadow, increment the shadow factor
                    // If the shadow factor is too high, skip the light
                    if (shadowIntersection.intersectionExists && shadowIntersection.t < shadowDistance) {
                        shadowFactor += 1.0f;
                        if (shadowFactor >= threshold) break; // Early exit if too much shadow
                    }
                }

                float lightVisibility = 1.0f - shadowFactor / numSamples;
                if (lightVisibility > 0.0f) {
                    color += computePhongComponents(lightDir, viewDir, normal, material, light) * lightVisibility;
                }
            }
        }

        // Add caustics effect
        color += renderCaustics(intersectionPoint, normal, material);

        Vec3 rayDirection = ray.direction();
        rayDirection.normalize();
        const Vec3 reflectedDirection = computeReflectedDirection(rayDirection, normal).normalize();
        const Vec3 bias = (Vec3::dot(ray.direction(), normal) < 0) ? normal * 1e-5f : -normal * 1e-5f;

        // If there are remaining bounces, compute the reflected ray
        if (material.type == Material_Mirror && NRemainingBounces > 0) {
            const Ray reflectedRay(intersectionPoint + bias, reflectedDirection);
            return rayTraceRecursive(reflectedRay, NRemainingBounces - 1, rng);
        }

        // If the material is glass, compute the refracted ray
        if (material.type == Material_Glass && NRemainingBounces > 0) {
            float eta = material.index_medium;
            Vec3 refractedDirection = computeRefractedDirection(rayDirection, normal, eta).normalize();
            const Ray refractedRay(intersectionPoint - bias, refractedDirection);

            // Compute Fresnel effect
            float fresnelEffect = computeFresnelEffect(rayDirection, normal, eta);

            // Trace both reflected and refracted rays
            Vec3 reflectedColor = rayTraceRecursive(Ray(intersectionPoint + bias, reflectedDirection), NRemainingBounces - 1, rng);
            Vec3 refractedColor = rayTraceRecursive(refractedRay, NRemainingBounces - 1, rng);

            // Combine the colors based on the Fresnel effect
            return reflectedColor * fresnelEffect + refractedColor * (1.0f - fresnelEffect);
        }

        return color; // Return the accumulated color
    }

    std::tuple<Vec3, Vec3, Material> handleIntersection(const RaySceneIntersection &intersection) const {
        Vec3 intersectionPoint, normal;
        Material material;

        // Determine the intersected object (sphere, square, mesh)
        if (intersection.typeOfIntersectedObject == 0) {
            // Sphere
            intersectionPoint = intersection.raySphereIntersection.intersection;
            normal = intersection.raySphereIntersection.normal;
            material = spheres[intersection.objectIndex].material;
        } else if (intersection.typeOfIntersectedObject == 1) {
            // Square
            intersectionPoint = intersection.raySquareIntersection.intersection;
            normal = intersection.raySquareIntersection.normal;
            material = squares[intersection.objectIndex].material;
        } else {
            // Mesh
            intersectionPoint = intersection.rayMeshIntersection.intersection;
            normal = intersection.rayMeshIntersection.normal;
            material = meshes[intersection.objectIndex].material;
        }

        return {intersectionPoint, normal, material};
    }

    Vec3 rayTrace(Ray const &rayStart, std::mt19937 &rng) {
        // Call the recursive function with the single ray and 1 bounce
        return rayTraceRecursive(rayStart, 5, rng, 4.8f);
    }

    static Vec3 randomDirection(std::mt19937 &rng) {
        std::uniform_real_distribution<float> dist(0.0f, 1.0f);
        const float theta = 2.0f * M_PIf * dist(rng);
        constexpr float coneAngle = M_PIf / 6.f * 2.5f; // Cone angle to limit the spread of the photons
        const float phi = coneAngle * dist(rng); // Limit phi between 0 and coneAngle

        float x = std::sin(phi) * std::cos(theta);
        float y = -std::cos(phi); // Negative to limit the direction to the upper hemisphere
        float z = std::sin(phi) * std::sin(theta);

        return {x, y, z};
    }

    void emitPhotons(const int photons) {
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
                RaySceneIntersection intersection = computeIntersection(ray, 0);

                if (intersection.intersectionExists) {
                    auto [intersectionPoint, normal, material] = handleIntersection(intersection);

                    // Glass material: refraction with transparency
                    if (material.type == Material_Glass) {
                        const Vec3 refractedDirection = computeRefractedDirection(photon.direction, normal, material.index_medium);
                        photon.position = intersectionPoint;
                        photon.direction = refractedDirection;

                        // Add object color to the photon color (with transparency)
                        photon.color = material.transparency * photon.color + (1 - material.transparency) * material.diffuse_material;
                        photonsGlass.emplace_back(photon);
                    }

                    // Mirror material: reflection
                    else if (material.type == Material_Mirror) {
                        photon.position = intersectionPoint + normal * 1e-4f; // small bias to avoid self-intersection
                        photon.direction = computeReflectedDirection(photon.direction, normal).normalize();
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

    Vec3 renderCaustics(const Vec3 &position, const Vec3 &normal, const Material &material) {
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

    // Scene 1
    void setup_single_sphere() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(-5, 5, 5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        } {
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0., 0., 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(0.811, 0.031, 0.129);
            s.material.specular_material = Vec3(0.2, 0.2, 0.2);
            s.material.shininess = 20;
        }
    }

    // Scene 2
    void setup_multiple_spheres() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear(); {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(-5, 5, 5);
            light.radius = 1.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        } {
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1., 0., 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(0.811, 0.031, 0.129);
            s.material.specular_material = Vec3(0.2, 0.2, 0.2);
            s.material.shininess = 20;
        } {
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1., 0., 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(0.129, 0.811, 0.031);
            s.material.specular_material = Vec3(0.2, 0.2, 0.2);
            s.material.shininess = 20;
        }
    }

    // Scene 3
    void setup_single_square() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear(); {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(-5, 5, 5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        } {
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.build_arrays();
            s.material.diffuse_material = Vec3(0.19, 0.40, 0.40); // Teal
            s.material.specular_material = Vec3(0.2, 0.2, 0.2);
            s.material.shininess = 20;
        }
    }

    void setup_cornell_box() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear(); {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(0.0, 1.5, 0.0);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Quad;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;

            // Quad for the light
            Mesh quad;
            quad.vertices.emplace_back(Vec3(-0.5, 1.5, -0.5), Vec3(0, -1, 0));
            quad.vertices.emplace_back(Vec3(0.5, 1.5, -0.5), Vec3(0, -1, 0));
            quad.vertices.emplace_back(Vec3(0.5, 1.5, 0.5), Vec3(0, -1, 0));
            quad.vertices.emplace_back(Vec3(-0.5, 1.5, 0.5), Vec3(0, -1, 0));
            quad.triangles.emplace_back(0, 1, 2);
            quad.triangles.emplace_back(0, 2, 3);
            quad.build_arrays();

            light.quad = quad;
        } {
            // Back Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(0.5, 0.5, 0.5);
            s.material.shininess = 16;
        } {
            // Left Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(194.0f / 255.0f, 49.0f / 255.0f, 44.0f / 255.0f);
            s.material.specular_material = Vec3(194.0f / 255.0f, 49.0f / 255.0f, 44.0f / 255.0f);
            s.material.shininess = 16;
        } {
            // Right Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(22.0f / 255.0f, 34.0f / 255.0f, 101.0f / 255.0f);
            s.material.specular_material = Vec3(22.0f / 255.0f, 34.0f / 255.0f, 101.0f / 255.0f);
            s.material.shininess = 16;
        } {
            // Floor => Checkerboard pattern
            const int numSquares = 8; // Number of squares per row/column
            const float squareSize = 2.0f / numSquares;
            squares.resize(squares.size() + numSquares * numSquares);
            const Vec3 color1(1.0f, 1.0f, 1.0f); // Blanc
            const Vec3 color2(0.0f, 0.0f, 0.0f); // Noir

            for (int i = 0; i < numSquares; ++i) {
                for (int j = 0; j < numSquares; ++j) {
                    Square &s = squares[squares.size() - numSquares * (i + 1) + j];
                    s.setQuad(
                        Vec3(-1 + j * squareSize, -1 + i * squareSize, 0.),
                        Vec3(squareSize, 0, 0),
                        Vec3(0, squareSize, 0),
                        squareSize, squareSize
                    );
                    s.translate(Vec3(0., 0., -2.));
                    s.scale(Vec3(2., 2., 1.));
                    s.rotate_x(-90);
                    s.build_arrays();
                    s.material.diffuse_material = ((i + j) % 2 == 0) ? color1 : color2;
                    s.material.specular_material = Vec3(0.2f, 0.2f, 0.2f);
                    s.material.shininess = 20;
                }
            }
        } {
            // Ceiling
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
            s.material.specular_material = Vec3(0.5, 0.5, 0.5);
            s.material.shininess = 16;
        } {
            // Front Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
            s.material.specular_material = Vec3(0.5, 0.5, 0.5);
            s.material.shininess = 16;
        } {
            // GLASS Sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3(194.0f / 255.0f, 49.0f / 255.0f, 44.0f / 255.0f);
            s.material.specular_material = Vec3(194.0f / 255.0f, 49.0f / 255.0f, 44.0f / 255.0f);
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.5;
        } {
            // MIRRORED Sphere
            spheres.resize(spheres.size() + 1);
            Sphere &s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }
    }
};

#endif
