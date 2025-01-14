#ifndef SPHERE_H
#define SPHERE_H

#include "Vec3.h"
#include "Ray.h"
#include <cmath>
#include <vector>

// -------------------------------------------
// Sphere Class
// -------------------------------------------

/**
 * Sphere class representing a 3D sphere mesh.
 */
class Sphere final : public Mesh {
public:
    Vec3 m_center;        ///< Center of the sphere.
    float m_radius{};     ///< Radius of the sphere.

    Sphere() : Mesh() {}
    Sphere(const Vec3 &c, const float r) : Mesh(), m_center(c), m_radius(r) {}

    /**
     * Builds the arrays for the sphere mesh.
     */
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
                uvs_array[2 * vertexIndex + 0] = std::fmod(u - 0.5f, 1.0f); // Add 0.5 to shift the texture to match OpenGL
                uvs_array[2 * vertexIndex + 1] = 1 - v; // Flip the V coordinate to match OpenGL
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

        computeAABB();
    }

    /**
     * Computes the Axis-Aligned Bounding Box (AABB) for the sphere.
     */
    void computeAABB() override {
        const Vec3 min = m_center - Vec3(m_radius, m_radius, m_radius);
        const Vec3 max = m_center + Vec3(m_radius, m_radius, m_radius);
        aabb = AABB(min, max);
    }

    /**
     * Checks if the ray intersects with the AABB of the sphere.
     * @param ray The ray to check for intersection.
     * @return True if the ray intersects with the AABB, false otherwise.
     */
    [[nodiscard]] bool intersectAABB(const Ray &ray) const override {
        return ray.intersectAABB(aabb);
    }

    /**
     * Computes the intersection of the ray with the sphere.
     * @param ray The ray to test for intersection.
     * @return A RayIntersection object containing the intersection details.
     */
    [[nodiscard]] RayIntersection intersect(const Ray &ray) const override {
        RayIntersection intersection;
        const Vec3 oc = ray.origin() - m_center; // Vector from sphere center to ray origin
        const Vec3 d = ray.direction(); // Ray direction

        const float a = Vec3::dot(d, d);
        const float b = 2.0f * Vec3::dot(d, oc);
        const float c = Vec3::dot(oc, oc) - m_radius * m_radius;
        const float delta = b * b - 4 * a * c; // Discriminant

        // No intersection
        if (delta < 0) {
            intersection.intersectionExists = false;
            return intersection;
        }

        intersection.intersectionExists = true;
        const float t1 = (-b - std::sqrt(delta)) / (2.0f * a);
        const float t2 = (-b + std::sqrt(delta)) / (2.0f * a);

        // Closest positive intersection
        if (t1 > 0 && t2 > 0) {
            intersection.t = std::min(t1, t2);
        } else if (t1 > 0) {
            intersection.t = t1;
        } else if (t2 > 0) {
            intersection.t = t2;
        } else {
            intersection.intersectionExists = false;
            return intersection;
        }

        // Calculate intersection point
        intersection.intersection = ray.origin() + intersection.t * ray.direction();

        // Calculate second intersection point if applicable NOTE: This is not used in the current implementation
        //intersection.secondIntersection = ray.origin() + t2 * ray.direction();

        intersection.u = 0.5f + std::atan2(intersection.intersection[1], intersection.intersection[0]) / (2.0f * M_PIf);
        intersection.v = 0.5f - std::asin((intersection.intersection[2]) / m_radius) / M_PIf;

        // Calculate the normal at the intersection point
        intersection.normal = (intersection.intersection - m_center) / m_radius;
        return intersection;
    }

private:
    /**
     * Converts spherical coordinates (theta, phi, radius) to Euclidean coordinates.
     * @param thetaPhiR Spherical coordinates.
     * @return Euclidean coordinates.
     */
    [[maybe_unused]] static Vec3 SphericalCoordinatesToEuclidean(const Vec3 &thetaPhiR) {
        return thetaPhiR[2] * Vec3(
            std::cos(thetaPhiR[0]) * std::cos(thetaPhiR[1]),
            std::sin(thetaPhiR[0]) * std::cos(thetaPhiR[1]),
            std::sin(thetaPhiR[1])
        );
    }

    /**
     * Converts spherical coordinates (theta, phi) to Euclidean coordinates.
     * @param theta Theta angle.
     * @param phi Phi angle.
     * @return Euclidean coordinates.
     */
    [[maybe_unused]] static Vec3 SphericalCoordinatesToEuclidean(const float theta, const float phi) {
        return {std::cos(theta) * std::cos(phi), std::sin(theta) * std::cos(phi), std::sin(phi)};
    }

    /**
     * Converts Euclidean coordinates to spherical coordinates.
     * @param xyz Euclidean coordinates.
     * @return Spherical coordinates.
     */
    [[maybe_unused]] static Vec3 EuclideanCoordinatesToSpherical(const Vec3 &xyz) {
        const float R = xyz.length();
        const float phi = std::asin(xyz[2] / R);
        const float theta = std::atan2(xyz[1], xyz[0]);
        return {theta, phi, R};
    }
};

#endif // SPHERE_H