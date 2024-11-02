#ifndef SQUARE_H
#define SQUARE_H

#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

// Structure to hold the results of a ray-square intersection
struct RaySquareIntersection : RayIntersection {
    float u{}, v{};
};

// Class representing a square in 3D space
class Square : public Mesh {
public:
    Vec3 m_normal; // Normal vector of the square
    Vec3 m_bottom_left; // Bottom-left corner of the square
    Vec3 m_right_vector; // Right vector of the square
    Vec3 m_up_vector; // Up vector of the square

    // Default constructor
    Square() : Mesh() {}
    Square(const Vec3 &bottomLeft, const Vec3 &rightVector, const Vec3 &upVector,
           const float width = 1.0f, const float height = 1.0f,
           const float uMin = 0.0f, const float uMax = 1.0f,
           const float vMin = 0.0f, const float vMax = 1.0f) : Mesh() {
        setQuad(bottomLeft, rightVector, upVector, width, height, uMin, uMax, vMin, vMax);
        Square::computeAABB();
    }

    // Set up the square with given parameters
    void setQuad(const Vec3 &bottomLeft, const Vec3 &rightVector, const Vec3 &upVector,
                 float const width = 1.0f, const float height = 1.0f,
                 float uMin = 0.0f, float uMax = 1.0f,
                 float vMin = 0.0f, float vMax = 1.0f) {
        m_right_vector = rightVector.normalized() * width; // Normalize and scale right vector
        m_up_vector = upVector.normalized() * height; // Normalize and scale up vector
        m_bottom_left = bottomLeft; // Set bottom-left corner
        m_normal = Vec3::cross(m_right_vector, m_up_vector).normalized(); // Calculate normal

        // Initialize vertices and texture coordinates
        vertices.clear();
        vertices.resize(4);
        vertices[0] = MeshVertex(bottomLeft, m_normal, uMin, vMin);
        vertices[1] = MeshVertex(bottomLeft + m_right_vector, m_normal, uMax, vMin);
        vertices[2] = MeshVertex(bottomLeft + m_right_vector + m_up_vector, m_normal, uMax, vMax);
        vertices[3] = MeshVertex(bottomLeft + m_up_vector, m_normal, uMin, vMax);

        // Set up triangle indices for rendering
        triangles.clear();
        triangles.resize(2);
        triangles[0] = {0, 1, 2}; // First triangle
        triangles[1] = {0, 2, 3}; // Second triangle
        Square::computeAABB();
    }

    void computeAABB() override {
        // Calculate the AABB of the square
        const Vec3 min = m_bottom_left;
        const Vec3 max = m_bottom_left + m_right_vector + m_up_vector;

        aabb = AABB(min, max);
    }

    [[nodiscard]] bool intersectAABB(const Ray &ray) const override {
        // Check if the ray intersects the square's AABB
        return ray.intersectAABB(aabb);
    }

    // Check for intersection between the ray and the square
    [[nodiscard]] RayIntersection intersect(const Ray &ray) const override {
        RaySquareIntersection intersection;
        intersection.intersectionExists = false; // Initialize to no intersection

        // Retrieve square data
        const Vec3 &m_bottom_left = vertices[0].position;
        const Vec3 m_right_vector = vertices[1].position - vertices[0].position;
        const Vec3 m_up_vector = vertices[3].position - vertices[0].position;

        // Calculate the normal of the square
        const Vec3 m_normal = Vec3::cross(m_right_vector, m_up_vector).normalized();

        // Calculate intersection with the square's plane
        const float denominator = Vec3::dot(m_normal, ray.direction());
        const float numerator = Vec3::dot(m_bottom_left - ray.origin(), m_normal);

        // If the ray is parallel to the plane, there's no intersection
        if (std::fabs(denominator) < 1e-6) {
            return intersection; // No intersection
        }

        // Calculate the distance t at which the ray intersects the plane
        const float t = numerator / denominator;

        // If t < 0, the intersection is behind the ray's origin
        if (t < 0) {
            return intersection; // No intersection
        }

        // Calculate the intersection point with the plane
        const Vec3 PointIntersection = ray.origin() + t * ray.direction();

        // Transform point into the square's local space
        const Vec3 localP = PointIntersection - m_bottom_left;

        // Calculate projections on the right and up vectors
        const float u = Vec3::dot(localP, m_right_vector) / m_right_vector.squareLength();
        const float v = Vec3::dot(localP, m_up_vector) / m_up_vector.squareLength();

        // Check if the intersection is within the square's bounds
        if (u >= 0 && u <= 1 && v >= 0 && v <= 1) {
            intersection.intersectionExists = true;
            intersection.t = t;
            intersection.u = u;
            intersection.v = v;
            intersection.intersection = PointIntersection;
            intersection.normal = m_normal; // Plane normal
        }

        return intersection; // Return the intersection result
    }
};

#endif // SQUARE_H
