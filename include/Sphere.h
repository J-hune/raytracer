#ifndef SPHERE_H
#define SPHERE_H

#include "Vec3.h"
#include "Ray.h"
#include <cmath>
#include <vector>

// Struct to hold the results of a ray-sphere intersection
struct RaySphereIntersection : RayIntersection {
    float theta{}, phi{};
    Vec3 secondIntersection;
};

// Converts spherical coordinates (theta, phi, radius) to Euclidean coordinates
[[maybe_unused]] static Vec3 SphericalCoordinatesToEuclidean(const Vec3 &thetaPhiR) {
    return thetaPhiR[2] * Vec3(
        std::cos(thetaPhiR[0]) * std::cos(thetaPhiR[1]),
        std::sin(thetaPhiR[0]) * std::cos(thetaPhiR[1]),
        std::sin(thetaPhiR[1])
    );
}

// Converts spherical coordinates (theta, phi) to Euclidean coordinates
[[maybe_unused]] static Vec3 SphericalCoordinatesToEuclidean(const float theta, const float phi) {
    return {std::cos(theta) * std::cos(phi), std::sin(theta) * std::cos(phi), std::sin(phi)};
}

// Converts Euclidean coordinates to spherical coordinates
[[maybe_unused]] static Vec3 EuclideanCoordinatesToSpherical(const Vec3 &xyz) {
    const float R = xyz.length();
    const float phi = std::asin(xyz[2] / R);
    const float theta = std::atan2(xyz[1], xyz[0]);
    return {theta, phi, R};
}

// Sphere class representing a 3D sphere mesh
class Sphere : public Mesh {
public:
    Vec3 m_center; // Center of the sphere
    float m_radius{}; // Radius of the sphere

    Sphere() : Mesh() {}
    Sphere(const Vec3 &c, const float r) : Mesh(), m_center(c), m_radius(r) {}

    // Builds arrays for vertices, normals, and texture coordinates
    void buildArrays() override {
        constexpr unsigned int nTheta = 20, nPhi = 20; // Number of subdivisions
        positions_array.resize(3 * nTheta * nPhi);
        normals_array.resize(3 * nTheta * nPhi);
        uvs_array.resize(2 * nTheta * nPhi);

        for (unsigned int thetaIt = 0; thetaIt < nTheta; ++thetaIt) {
            const float u = static_cast<float>(thetaIt) / (nTheta - 1);
            const float theta = u * 2.f * static_cast<float>(M_PI);
            for (unsigned int phiIt = 0; phiIt < nPhi; ++phiIt) {
                const unsigned int vertexIndex = thetaIt + phiIt * nTheta;
                const float v = static_cast<float>(phiIt) / (nPhi - 1);
                const auto phi = static_cast<float>(-M_PI / 2.0 + v * M_PI);
                Vec3 xyz = SphericalCoordinatesToEuclidean(theta, phi);

                // Store positions and normals
                positions_array[3 * vertexIndex + 0] = m_center[0] + m_radius * xyz[0];
                positions_array[3 * vertexIndex + 1] = m_center[1] + m_radius * xyz[1];
                positions_array[3 * vertexIndex + 2] = m_center[2] + m_radius * xyz[2];
                normals_array[3 * vertexIndex + 0] = xyz[0];
                normals_array[3 * vertexIndex + 1] = xyz[1];
                normals_array[3 * vertexIndex + 2] = xyz[2];
                uvs_array[2 * vertexIndex + 0] = u;
                uvs_array[2 * vertexIndex + 1] = v;
            }
        }

        // Build the triangle indices
        triangles_array.clear();
        for (unsigned int thetaIt = 0; thetaIt < nTheta - 1; ++thetaIt) {
            for (unsigned int phiIt = 0; phiIt < nPhi - 1; ++phiIt) {
                unsigned int vertexuv = thetaIt + phiIt * nTheta;
                unsigned int vertexUv = thetaIt + 1 + phiIt * nTheta;
                unsigned int vertexuV = thetaIt + (phiIt + 1) * nTheta;
                unsigned int vertexUV = thetaIt + 1 + (phiIt + 1) * nTheta;

                // Add two triangles for each quad
                triangles_array.push_back(vertexuv);
                triangles_array.push_back(vertexUv);
                triangles_array.push_back(vertexUV);
                triangles_array.push_back(vertexuv);
                triangles_array.push_back(vertexUV);
                triangles_array.push_back(vertexuV);
            }
        }

        Sphere::computeAABB();
    }

    void computeAABB() override {
        // Calculate the AABB of the sphere
        const Vec3 min = m_center - Vec3(m_radius, m_radius, m_radius);
        const Vec3 max = m_center + Vec3(m_radius, m_radius, m_radius);
        aabb = AABB(min, max);
    }

    [[nodiscard]] bool intersectAABB(const Ray &ray) const override {
        if (aabb.min == aabb.max) {
            exit(1);
        }
        // Check if the ray intersects the sphere's AABB
        return ray.intersectAABB(aabb);
    }

    // Checks for intersection between a ray and the sphere
    [[nodiscard]] RayIntersection intersect(const Ray &ray) const override {
        RaySphereIntersection intersection;
        const Vec3 oc = ray.origin() - m_center; // Vector from sphere center to ray origin
        const Vec3 d = ray.direction(); // Ray direction

        const float a = Vec3::dot(d, d); // d·d
        const float b = 2.0f * Vec3::dot(d, oc); // 2·t·d·(o-c)
        const float c = Vec3::dot(oc, oc) - m_radius * m_radius; // ||o-c||² - r²
        const float delta = b * b - 4 * a * c; // Discriminant

        // If the equation has no solution, there is no intersection
        if (delta < 0) {
            intersection.intersectionExists = false;
            return intersection;
        }

        // There is an intersection; we will look for the closest positive one
        intersection.intersectionExists = true;

        const float t1 = (-b - std::sqrt(delta)) / (2.0f * a); // First solution
        const float t2 = (-b + std::sqrt(delta)) / (2.0f * a); // Second solution

        // Determine the closest positive solution
        if (t1 > 0 && t2 > 0) {
            intersection.t = std::min(t1, t2);
        } else if (t1 > 0) {
            intersection.t = t1;
        } else if (t2 > 0) {
            intersection.t = t2;
        } else {
            intersection.intersectionExists = false; // Both solutions are negative
            return intersection;
        }

        // Calculate the intersection point
        intersection.intersection = ray.origin() + intersection.t * ray.direction();

        // Calculate the second intersection point if the ray is not tangent to the sphere
        if (t1 > 0 && t2 > 0) {
            intersection.secondIntersection = ray.origin() + t2 * ray.direction();
        }

        // Calculate the normal at the intersection point
        intersection.normal = (intersection.intersection - m_center) / m_radius;
        return intersection;
    }
};

#endif // SPHERE_H
