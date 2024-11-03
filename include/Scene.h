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
#include "KDTree.h"
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
    bool photonsEmitted = false;
    float directIlluminationReinhardKey = 0.0f;
    float causticsReinhardKey = 0.0f;
    bool drawCaustics = true;

public:
    Scene() = default;

    void draw() {
        const Settings &settings = Settings::getInstance();

        // Iterate over all objects and render them:
        for (const auto &mesh: meshes) mesh.draw();
        for (const auto &sphere: spheres) sphere.draw();
        for (const auto &square: squares) square.draw();

        if (settings.drawDebugPhotons) photonMap.debugDrawPhotons();
        if (!photonsEmitted && settings.caustics && drawCaustics) {
            std::cout << "Emitting photons, the application will freeze for a few seconds..." << std::endl;
            photonMap.emitPhotons(lights, spheres, squares, meshes, settings.photons);
            photonsEmitted = true;
        }

        if (settings.drawDebugAABBs) {
            for (const auto &sphere: spheres) sphere.aabb.draw();
            for (const auto &square: squares) square.aabb.draw();
            for (const auto &mesh: meshes) mesh.aabb.draw();
        }
    }

    Vec3 rayTrace(Ray const &rayStart, std::mt19937 &rng) {
        const Settings &settings = Settings::getInstance();

        // Call the recursive function with the single ray and 10 bounces
        return rayTraceRecursive(rayStart, 10, settings, rng, 4.8f);
    }

    // Recursive ray tracing function
    Vec3 rayTraceRecursive(const Ray &ray, const int NRemainingBounces, const Settings &settings, std::mt19937 &rng, const float z_near = 0.0f) {
        Vec3 color(0.0f, 0.0f, 0.0f);

        auto reinhardToneMapping = [](const Vec3 &pixelColor, const float keyValue) -> Vec3 {
            Vec3 mappedColor;
            for (int i = 0; i < 3; ++i) {
                mappedColor[i] = (pixelColor[i] * keyValue) / (pixelColor[i] * keyValue + 1.0f);
            }
            return mappedColor;
        };

        // Compute the intersection with the scene (ray vs. (spheres, squares, meshes))
        const RaySceneIntersection intersection = Intersection::computeIntersection(ray, spheres, squares, meshes, z_near);

        // If there isn't any intersection, return the background color
        if (!intersection.intersectionExists) {
            return color; // Background color
        }

        auto [intersectionPoint, normal, material] = Intersection::parseIntersection(intersection, spheres, squares, meshes);

        // Add direct illumination
        if (settings.directIllumination && material.type == Material_Diffuse_Blinn_Phong) {
            const Vec3 directIlluminationColor = computeDirectIllumination(ray, intersectionPoint, normal, material, settings, rng);
            color += directIlluminationReinhardKey > 0 ? reinhardToneMapping(directIlluminationColor, directIlluminationReinhardKey) : directIlluminationColor;
        }

        // Add caustics effect
        if (settings.caustics) {
            const Vec3 causticsColor = photonMap.computeCaustics(intersectionPoint, material);

            // Tone mapping
            color += reinhardToneMapping(causticsColor, causticsReinhardKey);
        }

        // Add reflection
        if (NRemainingBounces > 0 && settings.reflections) {
            const Vec3 reflectionColor = computeReflection(ray, intersectionPoint, normal, material, NRemainingBounces, settings, rng);
            color += reflectionColor;
        }

        // Add refraction
        if (NRemainingBounces > 0 && settings.refractions) {
            color += computeRefraction(ray, intersectionPoint, normal, material, NRemainingBounces, settings, rng);
        }

        return color; // Return the accumulated color
    }

    Vec3 computeDirectIllumination(const Ray &ray, const Vec3 &intersectionPoint, const Vec3 &normal, const Material &material, const Settings &settings, std::mt19937 &rng) {
        Vec3 color(0.0f, 0.0f, 0.0f);
        const Vec3 ambientLight(0.1f, 0.1f, 0.1f); // Ambient light

        // Add the ambient light to the color
        color += Vec3::compProduct(material.ambient_material, ambientLight);

        // Pre-normalize the view vector
        Vec3 viewDir = -ray.direction();
        viewDir.normalize();

        // For each light in the scene
        for (const auto &light: lights) {
            if (light.type == LightType_Spherical) {
                color += computeSphericalLight(intersectionPoint, normal, material, light, viewDir, settings);
            } else if (light.type == LightType_Quad) {
                color += computeQuadLight(intersectionPoint, normal, material, light, viewDir, settings, rng);
            }
        }

        return color;
    }

    [[nodiscard]] Vec3 computeSphericalLight(const Vec3 &intersectionPoint, const Vec3 &normal, const Material &material, const Light &light, const Vec3 &viewDir, const Settings &settings) const {
        constexpr float epsilon = 1e-5f;
        const Vec3 lightDir = (light.position - intersectionPoint).normalize();
        const float lightDistance = lightDir.length();

        // Ray to test the shadow
        const Ray shadowRay(intersectionPoint + normal * epsilon, lightDir);
        const RaySceneIntersection shadowIntersection = Intersection::computeIntersection(shadowRay, spheres, squares, meshes, epsilon);

        // If the point is in shadow, skip the light
        if (shadowIntersection.intersectionExists && shadowIntersection.t < lightDistance) {
            return Vec3(0.0f); // No contribution
        }

        return Lighting::computePhongComponents(lightDir, viewDir, normal, material, light);
    }

    Vec3 computeQuadLight(const Vec3 &intersectionPoint, const Vec3 &normal, const Material &material, const Light &light, const Vec3 &viewDir, const Settings &settings, std::mt19937 &rng) const {
        Vec3 color(0.0f);
        float shadowFactor = 0.0f;
        const float threshold = static_cast<float>(settings.shadowRays) * 0.9f;
        constexpr float epsilon = 1e-5f;
        const Vec3 lightDir = (light.position - intersectionPoint).normalize();

        // Sampling for soft shadows
        for (int i = 0; i < settings.shadowRays; ++i) {
            Vec3 samplePoint = Lighting::samplePointOnQuad(light, rng);
            Vec3 shadowDir = (samplePoint - intersectionPoint).normalize();
            const float shadowDistance = shadowDir.length();

            // Ray to test the shadow
            Ray shadowRay(intersectionPoint + normal * epsilon, shadowDir);
            const RaySceneIntersection shadowIntersection = Intersection::computeIntersection(
                shadowRay, spheres, squares, meshes, epsilon);

            // If the point is in shadow, increment the shadow factor
            if (shadowIntersection.intersectionExists && shadowIntersection.t < shadowDistance) {
                shadowFactor += 1.0f;
                if (shadowFactor >= threshold) break; // Early exit if too much shadow
            }
        }

        const float lightVisibility = 1.0f - shadowFactor / static_cast<float>(settings.shadowRays);
        if (lightVisibility > 0.0f) {
            color += Lighting::computePhongComponents(lightDir, viewDir, normal, material, light) * lightVisibility;
        }

        return color;
    }

    Vec3 computeReflection(const Ray &ray, const Vec3 &intersectionPoint, const Vec3 &normal, const Material &material, const int NRemainingBounces, const Settings &settings, std::mt19937 &rng) {
        Vec3 color(0.0f);

        if (material.type == Material_Mirror) {
            const Vec3 reflectedDirection = Lighting::computeReflectedDirection(ray.direction(), normal).normalize();
            const Vec3 bias = (Vec3::dot(ray.direction(), normal) < 0) ? normal * 1e-5f : -normal * 1e-5f;
            const Ray reflectedRay(intersectionPoint + bias, reflectedDirection);
            color += rayTraceRecursive(reflectedRay, NRemainingBounces - 1, settings, rng);
        }

        return color; // Return the accumulated reflection color
    }

    Vec3 computeRefraction(const Ray &ray, const Vec3 &intersectionPoint, const Vec3 &normal, const Material &material, const int NRemainingBounces, const Settings &settings, std::mt19937 &rng) {
        Vec3 color(0.0f);

        if (material.type == Material_Glass) {
            const float eta = material.index_medium;
            const Vec3 refractedDirection = Lighting::computeRefractedDirection(ray.direction(), normal, eta).normalize();
            const Vec3 reflectedDirection = Lighting::computeReflectedDirection(ray.direction(), normal).normalize();
            const Vec3 bias = (Vec3::dot(ray.direction(), normal) < 0) ? normal * 1e-5f : -normal * 1e-5f;
            const Ray refractedRay(intersectionPoint - bias, refractedDirection);

            // Compute Fresnel effect
            const float fresnelEffect = Lighting::computeFresnelEffect(ray.direction(), normal, eta);

            // Trace both reflected and refracted rays
            const Vec3 reflectedColor = rayTraceRecursive(Ray(intersectionPoint + bias, reflectedDirection), NRemainingBounces - 1, settings, rng);
            const Vec3 refractedColor = rayTraceRecursive(refractedRay, NRemainingBounces - 1, settings, rng);

            // Combine the colors based on the Fresnel effect
            color += reflectedColor * fresnelEffect + refractedColor * (1.0f - fresnelEffect);
        }

        return color; // Return the accumulated refraction color
    }


    /******************************************************************************************************************/
    /********************************************* SCENE SETUP FUNCTIONS **********************************************/
    /******************************************************************************************************************/

    // Scene 1
    void setup_single_sphere() {
        drawCaustics = false;
        addLight(Vec3(-5, 5, 5), LightType_Spherical, Vec3(1, 1, 1), 2.5f, 2.f);
        addSphere(Vec3(0., 0., 0.), 1.f, Material_Diffuse_Blinn_Phong, Vec3(0.811, 0.031, 0.129), Vec3(0.2, 0.2, 0.2), 20);
    }

    void setup_multiple_spheres() {
        drawCaustics = false;
        addLight(Vec3(-5, 5, 5), LightType_Spherical, Vec3(1, 1, 1), 1.5f, 2.f);
        addSphere(Vec3(1., 0., 0.), 1.f, Material_Diffuse_Blinn_Phong, Vec3(0.811, 0.031, 0.129), Vec3(0.2, 0.2, 0.2), 20);
        addSphere(Vec3(-1., 0., 0.), 1.f, Material_Diffuse_Blinn_Phong, Vec3(0.129, 0.811, 0.031), Vec3(0.2, 0.2, 0.2), 20);
    }

    void setup_single_square() {
        drawCaustics = false;
        addLight(Vec3(-5, 5, 5), LightType_Spherical, Vec3(1, 1, 1), 2.5f, 2.f);
        addSquare(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.,
                  Vec3(0.19, 0.40, 0.40), Vec3(0.2, 0.2, 0.2), 20);
    }

    void setup_cornell_box_with_2_spheres() {
        directIlluminationReinhardKey = 0.9f;
        causticsReinhardKey = 0.0008f;
        // 500000 photons => 0.0003f

        setup_cornell_box();
        addSphere(Vec3(1.0, -1.25, 0.5), 0.75f, Material_Glass, Vec3(0.f), Vec3(1.f), 16, 1.0, 1.5);
        addSphere(Vec3(-1.0, -1.25, -0.5), 0.75f, Material_Mirror, Vec3(0.f), Vec3(1.f), 16, 0.f, 0.f);
    }

    void setup_cornell_box_mesh() {
        directIlluminationReinhardKey = 0.9f;
        causticsReinhardKey = 0.0002f;
        // 500000 photons => 0.0001f

        setup_cornell_box();
        addMesh("../data/epcot.off", Vec3(0.0, -0.4, 0.0), Vec3(1.f),
            Material_Mirror, Vec3(0.f), Vec3(1.f), 16, 0.0, 1.5);
    }


    /******************************************************************************************************************/
    /********************************************* SCENE SETUP FUNCTIONS **********************************************/
    /******************************************************************************************************************/

     void addLight(const Vec3 &position, LightType type, const Vec3 &color, float intensity, float size) {
        lights.emplace_back(position, type, color, intensity, size, false);
    }

    void addSphere(const Vec3 &center, const float radius, const MaterialType materialType, const Vec3 &diffuseColor,
                   const Vec3 &specularColor, const float shininess, const float transparency = 0.f, const float indexMedium = 0.f) {
        Sphere s;
        s.m_center = center;
        s.m_radius = radius;
        s.buildArrays();
        s.material.type = materialType;
        s.material.diffuse_material = diffuseColor;
        s.material.specular_material = specularColor;
        s.material.shininess = shininess;
        s.material.transparency = transparency;
        s.material.index_medium = indexMedium;
        spheres.emplace_back(s);
    }

    void addSquare(const Vec3 &quadPos, const Vec3 &quadWidth, const Vec3 &quadHeight, const float scaleWidth, const float scaleHeight,
        const Vec3 &diffuseColor, const Vec3 &specularColor, const float shininess) {
        Square s;
        s.setQuad(quadPos, quadWidth, quadHeight, scaleWidth, scaleHeight);
        s.buildArrays();
        s.material.diffuse_material = diffuseColor;
        s.material.specular_material = specularColor;
        s.material.shininess = shininess;
        squares.emplace_back(s);
    }

    void addMesh(const std::string &filePath, const Vec3 &translation, const Vec3 &scale, const MaterialType materialType,
                 const Vec3 &diffuseColor, const Vec3 &specularColor, const float shininess,
                 const float transparency = 0.f, const float indexMedium = 0.f) {
        meshes.emplace_back();
        Mesh &m = meshes.back();
        m.loadOFF(filePath);
        m.scale(scale);
        m.translate(translation);
        m.buildArrays();
        m.material.type = materialType;
        m.material.diffuse_material = diffuseColor;
        m.material.specular_material = specularColor;
        m.material.shininess = shininess;
        m.material.transparency = transparency;
        m.material.index_medium = indexMedium;
    }

    void setup_cornell_box() {
        const Settings &settings = Settings::getInstance();

        const Light light(Vec3(0.0, 1.5, 0.0), LightType_Quad, Vec3(1, 1, 1), 2.5f, 2.f, false);
        light.quad = Mesh();
        light.quad.vertices.emplace_back(Vec3(-0.5, 1.5, -0.5), Vec3(0, -1, 0));
        light.quad.vertices.emplace_back(Vec3(0.5, 1.5, -0.5), Vec3(0, -1, 0));
        light.quad.vertices.emplace_back(Vec3(0.5, 1.5, 0.5), Vec3(0, -1, 0));
        light.quad.vertices.emplace_back(Vec3(-0.5, 1.5, 0.5), Vec3(0, -1, 0));
        light.quad.triangles.emplace_back(0, 1, 2);
        light.quad.triangles.emplace_back(0, 2, 3);
        light.quad.buildArrays();
        lights.push_back(light);

        {
            // Back Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., -2., -2.), Vec3(2., 0., 0.), Vec3(0., 2., 0.), 4., 4.);
            s.buildArrays();
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(0.05, 0.05, 0.05);
            s.material.shininess = 5;
        } {
            // Left Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., -2., 2.), Vec3(0., 0., -2.), Vec3(0., 2, 0.), 4., 4.);
            s.buildArrays();
            s.material.diffuse_material = Vec3(194.0f / 255.0f, 49.0f / 255.0f, 44.0f / 255.0f);
            s.material.specular_material = Vec3(194.0f / 255.0f, 49.0f / 255.0f, 44.0f / 255.0f) * 0.5f;
            s.material.shininess = 16;
        } {
            // Right Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(2., -2., -2.), Vec3(0., 0., 2.), Vec3(0., 2., 0.), 4., 4.);
            s.buildArrays();
            s.material.diffuse_material = Vec3(22.0f / 255.0f, 34.0f / 255.0f, 101.0f / 255.0f);
            s.material.specular_material = Vec3(22.0f / 255.0f, 34.0f / 255.0f, 101.0f / 255.0f) * 0.5f;
            s.material.shininess = 16;
        } {
            // Floor
            if (settings.floorType == CHECKERBOARD) {
                constexpr int numSquares = 8; // Number of squares per row/column
                constexpr float squareSize = 4.0f / numSquares;
                squares.resize(squares.size() + numSquares * numSquares);
                const Vec3 color1(1.0f, 1.0f, 1.0f); // Blanc
                const Vec3 color2(0.0f, 0.0f, 0.0f); // Noir

                for (int i = 0; i < numSquares; ++i) {
                    for (int j = 0; j < numSquares; ++j) {
                        Square &s = squares[squares.size() - numSquares * (i + 1) + j];
                        s.setQuad(
                            Vec3(2 - static_cast<float>(j) * squareSize, -2, -2 + static_cast<float>(i) * squareSize),
                            Vec3(-squareSize, 0, 0),
                            Vec3(0, 0, squareSize),
                            squareSize, squareSize
                        );
                        s.buildArrays();
                        s.material.diffuse_material = ((i + j) % 2 == 0) ? color1 : color2;
                        s.material.specular_material = Vec3(0.2f, 0.2f, 0.2f);
                        s.material.shininess = 20;
                    }
                }
            } else if (settings.floorType == PLAIN) {
                    squares.resize(squares.size() + 1);
                    Square &s = squares[squares.size() - 1];
                    s.setQuad(Vec3(2., -2., -2.), Vec3(-2., 0., 0.), Vec3(0., 0., 2.), 4., 4.);
                    s.buildArrays();
                    s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
                    s.material.specular_material = Vec3(0.5, 0.5, 0.5);
                    s.material.shininess = 16;
            }
        } {
            // Ceiling
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., 2., -2.), Vec3(2., 0., 0.), Vec3(0., 0., 2.), 4., 4.);
            s.buildArrays();
            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
            s.material.specular_material = Vec3(1.0, 1.0, 1.0) * 0.2f;
            s.material.shininess = 16;
        } {
            // Front Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-2., 2., 2.), Vec3(2., 0., 0.), Vec3(0., -2., 0.), 4., 4.);
            s.buildArrays();
            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
            s.material.specular_material = Vec3(0.05, 0.05, 0.05);
            s.material.shininess = 5;
        }
    }
};

#endif
