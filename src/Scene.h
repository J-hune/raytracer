#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"


#include <GL/glut.h>


enum LightType {
    LightType_Spherical,
    LightType_Quad
};


struct Light {
    Vec3 material;
    bool isInCamSpace{};
    LightType type;

    Vec3 pos;
    float radius{};

    Mesh quad;

    float powerCorrection;

    Light() : type(), powerCorrection(1.0) {
    }
};

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

public:
    Scene() = default;

    void draw() const {
        // On itère sur l'ensemble des objets, et faire leur rendu :
        for (const auto &mesh: meshes) {
            mesh.draw();
        }
        for (const auto &sphere: spheres) {
            sphere.draw();
        }
        for (const auto &square: squares) {
            square.draw();
        }
    }

    RaySceneIntersection computeIntersection(Ray const &ray, const float z_near) {
        RaySceneIntersection result;
        result.t = FLT_MAX;

        // Pour chaque objet de la scène, on calcule l'intersection avec le rayon
        for (unsigned int i = 0; i < spheres.size(); i++) {
            const RaySphereIntersection intersection = spheres[i].intersect(ray);

            // Si l'intersection existe et est plus proche que la précédente, on la garde
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

            // Si l'intersection existe et est plus proche que la précédente, on la garde
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

            // Si l'intersection existe et est plus proche que la précédente, on la garde
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

    Vec3 reflect(const Vec3 &I, const Vec3 &N) {
        return I - 2 * Vec3::dot(I, N) * N;
    }

    static Vec3 samplePointOnQuad(const Light &light) {
        const float u = static_cast<float>(rand()) / RAND_MAX;
        const float v = static_cast<float>(rand()) / RAND_MAX;

        Vec3 pointOnQuad = light.quad.vertices[0].position * (1 - u) * (1 - v) +
                           light.quad.vertices[1].position * u * (1 - v) +
                           light.quad.vertices[2].position * u * v +
                           light.quad.vertices[3].position * (1 - u) * v;

        return pointOnQuad;
    }

    Vec3 rayTraceRecursive(const Ray &ray, int NRemainingBounces, const float z_near) {
        Vec3 color(0.0, 0.0, 0.0);
        const Vec3 ambientLight(0.1, 0.1, 0.1); // Lumière ambiante
        const float epsilon = 1e-5f; // Pour éviter les auto-intersections

        // Calcul de l'intersection avec la scène
        const RaySceneIntersection intersection = computeIntersection(ray, z_near);

        if (intersection.intersectionExists) {
            Vec3 intersectionPoint;
            Vec3 normal;
            Material material;

            // On détermine l'objet intersecté (sphère, carré, mesh).
            if (intersection.typeOfIntersectedObject == 0) {
                intersectionPoint = intersection.raySphereIntersection.intersection;
                normal = intersection.raySphereIntersection.normal;
                material = spheres[intersection.objectIndex].material;
            } else if (intersection.typeOfIntersectedObject == 1) {
                intersectionPoint = intersection.raySquareIntersection.intersection;
                normal = intersection.raySquareIntersection.normal;
                material = squares[intersection.objectIndex].material;
            } else {
                intersectionPoint = intersection.rayMeshIntersection.intersection;
                normal = intersection.rayMeshIntersection.normal;
                material = meshes[intersection.objectIndex].material;
            }

            // Composante ambiante
            color += Vec3::compProduct(material.ambient_material, ambientLight);

            // Pré-normalisation du vecteur de vue
            Vec3 viewDir = -ray.direction();
            viewDir.normalize();

            // Pour chaque lumière dans la scène
            for (const auto &light: lights) {
                Vec3 lightDir = light.pos - intersectionPoint;
                float lightDistance = lightDir.length();
                lightDir.normalize();

                // Si la lumière est de type sphérique
                if (light.type == LightType_Spherical) {
                    // Vérification des ombres
                    Ray shadowRay(intersectionPoint + normal * epsilon, lightDir);
                    RaySceneIntersection shadowIntersection = computeIntersection(shadowRay, epsilon);

                    if (shadowIntersection.intersectionExists && shadowIntersection.t < lightDistance) {
                        continue; // Le point est dans l'ombre.
                    }

                    Vec3 reflectDir = reflect(-lightDir, normal);

                    // Composante diffuse
                    float diff = std::max(Vec3::dot(normal, lightDir), 0.0f);
                    Vec3 diffuse = Vec3::compProduct(material.diffuse_material, light.material) * diff;

                    // Composante spéculaire (approximation si besoin)
                    const float spec = std::pow(std::max(Vec3::dot(viewDir, reflectDir), 0.0f), material.shininess);
                    Vec3 specular = Vec3::compProduct(material.specular_material, light.material) * spec;

                    // Ajout des composantes diffuse et spéculaire à la couleur
                    color += diffuse + specular;
                }

                // Si la lumière est de type quad
                else if (light.type == LightType_Quad) {
                    const int numSamples = 16; // Nombre d'échantillons pour les ombres douces
                    float shadowFactor = 0.0f;
                    const float threshold = 0.8f * numSamples; // Early exit

                    // Échantillonnage pour les ombres douces
                    for (int i = 0; i < numSamples; ++i) {
                        Vec3 samplePoint = samplePointOnQuad(light);
                        Vec3 shadowDir = samplePoint - intersectionPoint;
                        float shadowDistance = shadowDir.length();
                        shadowDir.normalize();

                        // Rayon pour tester l'ombre
                        Ray shadowRay(intersectionPoint + normal * epsilon, shadowDir);
                        RaySceneIntersection shadowIntersection = computeIntersection(shadowRay, epsilon);

                        if (shadowIntersection.intersectionExists && shadowIntersection.t < shadowDistance) {
                            shadowFactor += 1.0f;
                            if (shadowFactor >= threshold) break; // Sortie anticipée si trop d'ombre
                        }
                    }

                    float lightVisibility = 1.0f - shadowFactor / numSamples;
                    if (lightVisibility > 0.0f) {
                        // Composante diffuse
                        float diff = std::max(Vec3::dot(normal, lightDir), 0.0f);
                        Vec3 diffuse = Vec3::compProduct(material.diffuse_material, light.material) * diff *
                                       lightVisibility;

                        // Composante spéculaire
                        Vec3 reflectDir = reflect(-lightDir, normal);
                        const float spec = std::pow(std::max(Vec3::dot(viewDir, reflectDir), 0.0f), material.shininess);
                        Vec3 specular = Vec3::compProduct(material.specular_material, light.material) * spec *
                                        lightVisibility;

                        // Ajout des composantes diffuse et spéculaire à la couleur
                        color += diffuse + specular;
                    }
                }
            }
        } else {
            color += Vec3(1.0, 1.0, 1.0); // Couleur de fond
        }

        return color;
    }


    Vec3 rayTrace(Ray const &rayStart) {
        // On appelle la fonction récursive avec l'unique rayon et 1 rebond
        Vec3 color = rayTraceRecursive(rayStart, 1, 4.8f);
        return color;
    }

    // Scène 1
    void setup_single_sphere() {
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

    // Scène 2
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

    // Scène 3
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
            s.material.diffuse_material = Vec3(0.19, 0.40, 0.40); // Vert sarcelle
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

            // On définit les sommets du quad
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
            //Back Wall
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3(1., 1., 1.);
            s.material.specular_material = Vec3(1., 1., 1.);
            s.material.shininess = 16;
        } {
            //Left Wall
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
            //Right Wall
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
            // Floor
            squares.resize(squares.size() + 1);
            Square &s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
            s.material.shininess = 16;
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
            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
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
            s.material.specular_material = Vec3(1.0, 1.0, 1.0);
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
            s.material.index_medium = 1.4;
        } {
            //MIRRORED Sphere
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
