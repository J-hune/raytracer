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

static Vec3 SphericalCoordinatesToEuclidean(Vec3 ThetaPhiR) {
    return ThetaPhiR[2] * Vec3(
        std::cos(ThetaPhiR[0]) * std::cos(ThetaPhiR[1]),
        std::sin(ThetaPhiR[0]) * std::cos(ThetaPhiR[1]),
        std::sin(ThetaPhiR[1])
    );
}

static Vec3 SphericalCoordinatesToEuclidean(const float theta, const float phi) {
    return {std::cos(theta) * std::cos(phi), std::sin(theta) * std::cos(phi), std::sin(phi)};
}

static Vec3 EuclideanCoordinatesToSpherical(Vec3 xyz) {
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

    void build_arrays() override {
        unsigned int nTheta = 20, nPhi = 20;
        positions_array.resize(3 * nTheta * nPhi);
        normalsArray.resize(3 * nTheta * nPhi);
        uvs_array.resize(2 * nTheta * nPhi);
        for (unsigned int thetaIt = 0; thetaIt < nTheta; ++thetaIt) {
            const float u = static_cast<float>(thetaIt) / static_cast<float>(nTheta - 1);
            const float theta = u * 2.f * static_cast<float>(M_PI);
            for (unsigned int phiIt = 0; phiIt < nPhi; ++phiIt) {
                const unsigned int vertexIndex = thetaIt + phiIt * nTheta;
                const float v = static_cast<float>(phiIt) / static_cast<float>(nPhi - 1);
                const auto phi = static_cast<float>(-M_PI / 2.0 + v * M_PI);
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
        // Sphere equation ||x-c||² - r² = 0 (with x being the ray equation)
        // Ray equation t² d.d + 2t d.(o-c) + ||o-c||² - r² = 0
        const Vec3 oc = ray.origin() - m_center; // o-c
        const Vec3 d = ray.direction(); // d

        const float a = Vec3::dot(d, d); // d.d
        const float b = 2.0f * Vec3::dot(d, oc); // 2t d.(o-c)
        const float c = Vec3::dot(oc, oc) - m_radius * m_radius; // ||o-c||² - r²
        const float delta = b * b - 4 * a * c;

        // If the equation has no solution, there is no intersection
        if (delta < 0) {
            intersection.intersectionExists = false;
            return intersection;
        }

        // We know there is an intersection, we will look for the closest positive one
        intersection.intersectionExists = true;

        const float t1 = (-b - sqrtf(delta)) / (2.0f * a); // Solution 1: (-b - sqrt(delta)) / 2a    <= Generally the closest solution
        const float t2 = (-b + sqrtf(delta)) / (2.0f * a); // Solution 2: (-b + sqrt(delta)) / 2a

        // If both solutions are positive, we take the closest one else we take the positive one
        if (t1 > 0 && t2 > 0) {
            intersection.t = std::min(t1, t2);
        } else if (t1 > 0) {
            intersection.t = t1;
        } else if (t2 > 0) {
            intersection.t = t2;
        } else {
            // If both solutions are negative, there is no intersection (finally, I'm not a crook)
            intersection.intersectionExists = false;
            return intersection;
        }

        // Calculate the intersection point
        intersection.intersection = ray.origin() + intersection.t * ray.direction();

        // Calculate the second intersection point if the ray is not tangent to the sphere
        if (t1 > 0 && t2 > 0) {
            intersection.secondintersection = ray.origin() + t2 * ray.direction();
        }

        // Calculate the normal at the intersection point
        intersection.normal = (intersection.intersection - m_center) / m_radius;
        return intersection;
    }
};
#endif
