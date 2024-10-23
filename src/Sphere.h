#ifndef Sphere_H
#define Sphere_H

#include "Vec3.h"
#include <cmath>
#include <vector>

struct RaySphereIntersection {
    bool intersectionExists{};
    float t{};
    float theta{}, phi{};
    Vec3 intersection;
    Vec3 secondintersection;
    Vec3 normal;
};

static
Vec3 SphericalCoordinatesToEuclidean(Vec3 ThetaPhiR) {
    return ThetaPhiR[2] * Vec3(
        std::cos(ThetaPhiR[0]) * std::cos(ThetaPhiR[1]),
        std::sin(ThetaPhiR[0]) * std::cos(ThetaPhiR[1]),
        std::sin(ThetaPhiR[1])
    );
}

static
Vec3 SphericalCoordinatesToEuclidean(const float theta, const float phi) {
    return {std::cos(theta) * std::cos(phi), std::sin(theta) * std::cos(phi), std::sin(phi)};
}

static
Vec3 EuclideanCoordinatesToSpherical(Vec3 xyz) {
    const float R = xyz.length();
    const float phi = std::asin(xyz[2] / R);
    const float theta = std::atan2(xyz[1], xyz[0]);
    return {theta, phi, R};
}


class Sphere : public Mesh {
public:
    Vec3 m_center;
    float m_radius{};

    Sphere() : Mesh() {}
    Sphere(const Vec3 &c, const float r) : Mesh(), m_center(c), m_radius(r) {}

    void build_arrays() {
        unsigned int nTheta = 20, nPhi = 20;
        positions_array.resize(3 * nTheta * nPhi);
        normalsArray.resize(3 * nTheta * nPhi);
        uvs_array.resize(2 * nTheta * nPhi);
        for (unsigned int thetaIt = 0; thetaIt < nTheta; ++thetaIt) {
            float u = (float) (thetaIt) / (float) (nTheta - 1);
            float theta = u * 2 * M_PI;
            for (unsigned int phiIt = 0; phiIt < nPhi; ++phiIt) {
                unsigned int vertexIndex = thetaIt + phiIt * nTheta;
                float v = (float) (phiIt) / (float) (nPhi - 1);
                float phi = -M_PI / 2.0 + v * M_PI;
                Vec3 xyz = SphericalCoordinatesToEuclidean(theta, phi);
                positions_array[3 * vertexIndex + 0] = m_center[0] + m_radius * xyz[0];
                positions_array[3 * vertexIndex + 1] = m_center[1] + m_radius * xyz[1];
                positions_array[3 * vertexIndex + 2] = m_center[2] + m_radius * xyz[2];
                normalsArray[3 * vertexIndex + 0] = xyz[0];
                normalsArray[3 * vertexIndex + 1] = xyz[1];
                normalsArray[3 * vertexIndex + 2] = xyz[2];
                uvs_array[2 * vertexIndex + 0] = u;
                uvs_array[2 * vertexIndex + 1] = v;
            }
        }
        triangles_array.clear();
        for (unsigned int thetaIt = 0; thetaIt < nTheta - 1; ++thetaIt) {
            for (unsigned int phiIt = 0; phiIt < nPhi - 1; ++phiIt) {
                unsigned int vertexuv = thetaIt + phiIt * nTheta;
                unsigned int vertexUv = thetaIt + 1 + phiIt * nTheta;
                unsigned int vertexuV = thetaIt + (phiIt + 1) * nTheta;
                unsigned int vertexUV = thetaIt + 1 + (phiIt + 1) * nTheta;
                triangles_array.push_back(vertexuv);
                triangles_array.push_back(vertexUv);
                triangles_array.push_back(vertexUV);
                triangles_array.push_back(vertexuv);
                triangles_array.push_back(vertexUV);
                triangles_array.push_back(vertexuV);
            }
        }
    }


    RaySphereIntersection intersect(const Ray &ray) const {
        RaySphereIntersection intersection;
        // Equation de la sphère ||x-c||² - r² = 0 (avec x étant l'équation du rayon)
        // Equation du rayon t² d.d + 2t d.(o-c) + ||o-c||² - r² = 0
        const Vec3 oc = ray.origin() - m_center; // o-c
        const Vec3 d = ray.direction(); // d

        const float a = Vec3::dot(d, d); // d.d
        const float b = 2.0f * Vec3::dot(d, oc); // 2t d.(o-c)
        const float c = Vec3::dot(oc, oc) - m_radius * m_radius; // ||o-c||² - r²
        const float delta = b * b - 4 * a * c;

        // Si l'équation n'a pas de solution, il n'y a pas d'intersection
        if (delta < 0) {
            intersection.intersectionExists = false;
            return intersection;
        }

        // On sait qu'il y a une intersection, on va chercher la plus proche positive
        intersection.intersectionExists = true;

        const float t1 = (-b - sqrtf(delta)) / (2.0f * a); // Solution 1 : (-b - sqrt(delta)) / 2a    <= Généralement la solution la plus proche
        const float t2 = (-b + sqrtf(delta)) / (2.0f * a); // Solution 2 : (-b + sqrt(delta)) / 2a

        // On garde la plus proche positive
        if (t1 > 0 && t2 > 0) {
            intersection.t = std::min(t1, t2);
        } else if (t1 > 0) {
            intersection.t = t1;
        } else if (t2 > 0) {
            intersection.t = t2;
        } else {
            // Si les deux solutions sont négatives, il n'y a pas d'intersection (finalement, je ne suis pas un margoulin)
            intersection.intersectionExists = false;
            return intersection;
        }

        // On calcule le point d'intersection
        intersection.intersection = ray.origin() + intersection.t * ray.direction();

        // On calcule le deuxième point d'intersection s'il existe
        if (t1 > 0 && t2 > 0) {
            intersection.secondintersection = ray.origin() + t2 * ray.direction();
        }

        // On calcule la normale
        intersection.normal = (intersection.intersection - m_center) / m_radius;
        return intersection;
    }
};
#endif
