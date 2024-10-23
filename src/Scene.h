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

    RaySceneIntersection() : intersectionExists(false), t(FLT_MAX), rayMeshIntersection(), raySquareIntersection() {}
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
        for (const auto & mesh : meshes) {
            mesh.draw();
        }
        for (const auto & sphere : spheres) {
            sphere.draw();
        }
        for (const auto & square : squares) {
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

    Vec3 rayTraceRecursive(const Ray &ray, int NRemainingBounces, const float z_near) {
        Vec3 color;

        // On calcule l'intersection avec la scène
        const RaySceneIntersection intersection = computeIntersection(ray, z_near);

        if (intersection.intersectionExists) {
            // Si l'objet intersecté est une sphère
            if (intersection.typeOfIntersectedObject == 0) {
                color += spheres[intersection.objectIndex].material.diffuse_material;
            }

            // Si l'objet intersecté est un carré
            else if (intersection.typeOfIntersectedObject == 1) {
                color += squares[intersection.objectIndex].material.diffuse_material;
            }

            // L'objet intersecté est donc un mesh
            else {
                color += meshes[intersection.objectIndex].material.diffuse_material;
            }
        } else {
            color += Vec3(1.0, 1.0, 1.0);
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
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
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
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
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
