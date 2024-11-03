#ifndef SQUARE_H
#define SQUARE_H

#include "Vec3.h"
#include "Mesh.h"
#include <vector>
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
        computeAABB();
    }

    // Set up the square with given parameters
    void setQuad(const Vec3 &bottomLeft, const Vec3 &rightVector, const Vec3 &upVector,
                 float width = 1.0f, float height = 1.0f,
                 float uMin = 0.0f, float uMax = 1.0f,
                 float vMin = 0.0f, float vMax = 1.0f) {
        m_right_vector = rightVector.normalized() * width;
        m_up_vector = upVector.normalized() * height;
        m_bottom_left = bottomLeft;
        m_normal = Vec3::cross(m_right_vector, m_up_vector).normalized();

        initializeVertices(uMin, uMax, vMin, vMax);
        initializeTriangles();
        computeAABB();
    }

    void computeAABB() override {
        const Vec3 min = m_bottom_left;
        const Vec3 max = m_bottom_left + m_right_vector + m_up_vector;
        aabb = AABB(min, max);
    }

    [[nodiscard]] bool intersectAABB(const Ray &ray) const override {
        return ray.intersectAABB(aabb);
    }

    // Check for intersection between the ray and the square
    [[nodiscard]] RayIntersection intersect(const Ray &ray) const override {
        RaySquareIntersection intersection;
        intersection.intersectionExists = false;

        const Vec3 &bottomLeft = vertices[0].position;
        const Vec3 rightVector = vertices[1].position - bottomLeft;
        const Vec3 upVector = vertices[3].position - bottomLeft;

        const Vec3 normal = Vec3::cross(rightVector, upVector).normalized();
        const float denominator = Vec3::dot(normal, ray.direction());
        const float numerator = Vec3::dot(bottomLeft - ray.origin(), normal);

        // Check for parallelism
        if (std::fabs(denominator) < 1e-6) {
            return intersection; // No intersection
        }

        const float t = numerator / denominator;
        if (t < 0) {
            return intersection; // No intersection
        }

        const Vec3 PointIntersection = ray.origin() + t * ray.direction();
        const Vec3 localP = PointIntersection - bottomLeft;

        const float u = Vec3::dot(localP, rightVector) / rightVector.squareLength();
        const float v = Vec3::dot(localP, upVector) / upVector.squareLength();

        // Check bounds of the square
        if (u >= 0 && u <= 1 && v >= 0 && v <= 1) {
            intersection.intersectionExists = true;
            intersection.t = t;
            intersection.u = u;
            intersection.v = v;
            intersection.intersection = PointIntersection;
            intersection.normal = normal; // Plane normal
        }

        return intersection; // Return the intersection result
    }

private:
    // Initialize vertices of the square
    void initializeVertices(float uMin, float uMax, float vMin, float vMax) {
        vertices.clear();
        vertices.resize(4);
        vertices[0] = MeshVertex(m_bottom_left, m_normal, uMin, vMin);
        vertices[1] = MeshVertex(m_bottom_left + m_right_vector, m_normal, uMax, vMin);
        vertices[2] = MeshVertex(m_bottom_left + m_right_vector + m_up_vector, m_normal, uMax, vMax);
        vertices[3] = MeshVertex(m_bottom_left + m_up_vector, m_normal, uMin, vMax);
    }

    // Initialize triangles for rendering
    void initializeTriangles() {
        triangles.clear();
        triangles.resize(2);
        triangles[0] = {0, 1, 2}; // First triangle
        triangles[1] = {0, 2, 3}; // Second triangle
    }
};

#endif // SQUARE_H
