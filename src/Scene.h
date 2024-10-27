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
#include "Settings.h"
#include "Intersection.h"
#include "Lighting.h"
#include "PhotonMap.h"

class Scene {
    std::vector<Mesh> meshes;
    std::vector<Sphere> spheres;
    std::vector<Square> squares;
    std::vector<Light> lights;
    PhotonMap photonMap;

public:
    Scene() = default;

    void draw() const {
        const Settings &settings = Settings::getInstance();

        // Iterate over all objects and render them:
        for (const auto &mesh: meshes) mesh.draw();
        for (const auto &sphere: spheres) sphere.draw();
        for (const auto &square: squares) square.draw();

        if (settings.drawDebugPhotons) photonMap.debugDrawPhotons();
    }

    // Recursive ray tracing function
    Vec3 rayTraceRecursive(const Ray &ray, const int NRemainingBounces, const Settings &settings, std::mt19937 &rng, const float z_near = 0.0f) {
        Vec3 color(0.0f, 0.0f, 0.0f);
        const Vec3 ambientLight(0.1f, 0.1f, 0.1f); // Ambient light

        // Compute the intersection with the scene (ray vs. (spheres, squares, meshes))
        const RaySceneIntersection intersection = Intersection::computeIntersection(ray, spheres, squares, meshes, z_near);

        // If there isn't any intersection, return the background color
        if (!intersection.intersectionExists || NRemainingBounces == 0) {
            return color;
        }

        auto [intersectionPoint, normal, material] = Intersection::parseIntersection(intersection, spheres, squares, meshes);

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
                RaySceneIntersection shadowIntersection = Intersection::computeIntersection(shadowRay, spheres, squares, meshes, epsilon);

                // If the point is in shadow, skip the light
                if (shadowIntersection.intersectionExists && shadowIntersection.t < lightDistance) {
                    continue; // The point is in shadow.
                }

                color += Lighting::computePhongComponents(lightDir, viewDir, normal, material, light);
            }

            // If the light is a quad (area light) => soft shadows
            else if (light.type == LightType_Quad) {
                constexpr int numSamples = 16; // Number of samples for soft shadows
                float shadowFactor = 0.0f;
                constexpr float threshold = numSamples * 0.9f;

                // Sampling for soft shadows
                for (int i = 0; i < numSamples; ++i) {
                    Vec3 samplePoint = Lighting::samplePointOnQuad(light, rng);
                    Vec3 shadowDir = (samplePoint - intersectionPoint).normalize();
                    float shadowDistance = shadowDir.length();

                    // Ray to test the shadow
                    Ray shadowRay(intersectionPoint + normal * epsilon, shadowDir);
                    RaySceneIntersection shadowIntersection = Intersection::computeIntersection(shadowRay, spheres, squares, meshes, epsilon);

                    // If the point is in shadow, increment the shadow factor
                    // If the shadow factor is too high, skip the light
                    if (shadowIntersection.intersectionExists && shadowIntersection.t < shadowDistance) {
                        shadowFactor += 1.0f;
                        if (shadowFactor >= threshold) break; // Early exit if too much shadow
                    }
                }

                float lightVisibility = 1.0f - shadowFactor / numSamples;
                if (lightVisibility > 0.0f) {
                    color += Lighting::computePhongComponents(lightDir, viewDir, normal, material, light) * lightVisibility;
                }
            }
        }

        // Add caustics effect
        if (settings.caustics) {
            color += photonMap.renderCaustics(intersectionPoint, normal, material);
        }

        Vec3 rayDirection = ray.direction();
        rayDirection.normalize();
        const Vec3 reflectedDirection = Lighting::computeReflectedDirection(rayDirection, normal).normalize();
        const Vec3 bias = (Vec3::dot(ray.direction(), normal) < 0) ? normal * 1e-5f : -normal * 1e-5f;

        // If there are remaining bounces, compute the reflected ray
        if (material.type == Material_Mirror && NRemainingBounces > 0) {
            const Ray reflectedRay(intersectionPoint + bias, reflectedDirection);
            return rayTraceRecursive(reflectedRay, NRemainingBounces - 1, settings, rng);
        }

        // If the material is glass, compute the refracted ray
        if (material.type == Material_Glass && NRemainingBounces > 0) {
            float eta = material.index_medium;
            Vec3 refractedDirection = Lighting::computeRefractedDirection(rayDirection, normal, eta).normalize();
            const Ray refractedRay(intersectionPoint - bias, refractedDirection);

            // Compute Fresnel effect
            float fresnelEffect = Lighting::computeFresnelEffect(rayDirection, normal, eta);

            // Trace both reflected and refracted rays
            Vec3 reflectedColor = rayTraceRecursive(Ray(intersectionPoint + bias, reflectedDirection), NRemainingBounces - 1, settings, rng);
            Vec3 refractedColor = rayTraceRecursive(refractedRay, NRemainingBounces - 1, settings, rng);

            // Combine the colors based on the Fresnel effect
            return reflectedColor * fresnelEffect + refractedColor * (1.0f - fresnelEffect);
        }

        return color; // Return the accumulated color
    }

    Vec3 rayTrace(Ray const &rayStart, std::mt19937 &rng) {
        const Settings &settings = Settings::getInstance();

        // Call the recursive function with the single ray and 1 bounce
        return rayTraceRecursive(rayStart, 5, settings, rng, 4.8f);
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


    /******************************************************************************************************************/
    /********************************************* SCENE SETUP FUNCTIONS **********************************************/
    /******************************************************************************************************************/

    void add_light(Vec3 position, float radius, float powerCorrection, LightType type, Vec3 material, bool isInCamSpace) {
        lights.emplace_back();
        Light &light = lights.back();
        light.pos = position;
        light.radius = radius;
        light.powerCorrection = powerCorrection;
        light.type = type;
        light.material = material;
        light.isInCamSpace = isInCamSpace;
    }

    // Scene 1
    void setup_single_sphere() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        add_light(Vec3(-5, 5, 5), 2.5f, 2.f, LightType_Spherical, Vec3(1, 1, 1), false);
        {
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
        lights.clear();
        add_light(Vec3(-5, 5, 5), 1.5f, 2.f, LightType_Spherical, Vec3(1, 1, 1), false);
        {
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
        lights.clear();
        add_light(Vec3(-5, 5, 5), 2.5f, 2.f, LightType_Spherical, Vec3(1, 1, 1), false);
        {
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
        const Settings &settings = Settings::getInstance();

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
            s.material.specular_material = Vec3(0.05, 0.05, 0.05);
            s.material.shininess = 5;
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
            s.material.specular_material = Vec3(194.0f / 255.0f, 49.0f / 255.0f, 44.0f / 255.0f) * 0.5f;
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
            s.material.specular_material = Vec3(22.0f / 255.0f, 34.0f / 255.0f, 101.0f / 255.0f) * 0.5f;
            s.material.shininess = 16;
        } {
            // Floor
            if (settings.floorType == CHECKERBOARD) {
                constexpr int numSquares = 8; // Number of squares per row/column
                constexpr float squareSize = 2.0f / numSquares;
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
            } else if (settings.floorType == PLAIN) {
                    squares.resize(squares.size() + 1);
                    Square &s = squares[squares.size() - 1];
                    s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
                    s.translate(Vec3(0., 0., -2.));
                    s.scale(Vec3(2., 2., 1.));
                    s.rotate_x(-90);
                    s.build_arrays();
                    s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
                    s.material.specular_material = Vec3(0.5, 0.5, 0.5);
                    s.material.shininess = 16;
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
            s.material.specular_material = Vec3(1.0, 1.0, 1.0) * 0.2f;
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
            s.material.specular_material = Vec3(0.05, 0.05, 0.05);
            s.material.shininess = 5;
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

        if (settings.caustics) {
            photonMap.emitPhotons(lights, spheres, squares, meshes, settings.photons);
        }
    }
};

#endif
