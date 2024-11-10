#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "Vec3.h"
#include "Ray.h"

/**
 * Structure for ray-triangle intersection results.
 */
struct RayTriangleIntersection : RayIntersection {
    unsigned int tIndex{};  ///< Triangle index.
};

/**
 * Class representing a triangle.
 */
class Triangle {
private:
    Vec3 m_c[3];           ///< Vertices of the triangle.
    Vec3 m_normal;         ///< Normal vector of the triangle.
    float area{};          ///< Area of the triangle.

public:
    Triangle() = default;
    Triangle(Vec3 const &c0, Vec3 const &c1, Vec3 const &c2) {
        m_c[0] = c0;
        m_c[1] = c1;
        m_c[2] = c2;
        updateAreaAndNormal(); ///< Update area and normal on construction.
    }

    /**
     * Get a vertex of the triangle.
     * @param i Index of the vertex.
     * @return The vertex at the specified index.
     */
    [[nodiscard]] Vec3 const &getVertex(const unsigned int i) const {
        return m_c[i];
    }

    /**
     * Update the area and normal vector of the triangle.
     */
    void updateAreaAndNormal() {
        const Vec3 nNotNormalized = Vec3::cross(m_c[1] - m_c[0], m_c[2] - m_c[0]);
        const float norm = nNotNormalized.length();
        m_normal = nNotNormalized / norm;
        area = norm / 2.f;
    }

    /**
     * Set the first vertex of the triangle.
     * @param c0 Position of the first vertex.
     */
    void setC0(Vec3 const &c0) { m_c[0] = c0; updateAreaAndNormal(); }

    /**
     * Set the second vertex of the triangle.
     * @param c1 Position of the second vertex.
     */
    void setC1(Vec3 const &c1) { m_c[1] = c1; updateAreaAndNormal(); }

    /**
     * Set the third vertex of the triangle.
     * @param c2 Position of the third vertex.
     */
    void setC2(Vec3 const &c2) { m_c[2] = c2; updateAreaAndNormal(); }

    /**
     * Get the normal vector of the triangle.
     * @return The normal vector of the triangle.
     */
    [[nodiscard]] Vec3 const &normal() const {
        return m_normal;
    }

    /**
     * Get the Axis-Aligned Bounding Box (AABB) of the triangle.
     * @return The AABB of the triangle.
     */
    [[nodiscard]] AABB getAABB() const {
        const Vec3 min = Vec3::min(m_c[0], Vec3::min(m_c[1], m_c[2]));
        const Vec3 max = Vec3::max(m_c[0], Vec3::max(m_c[1], m_c[2]));
        return {min, max};
    }

    /**
     * Get the combined AABB of a vector of triangles.
     * @param triangles Vector of triangles.
     * @return The combined AABB of the triangles.
     */
    static AABB getAABB(std::vector<Triangle> const &triangles) {
        Vec3 min(FLT_MAX);
        Vec3 max(-FLT_MAX);

        for (const Triangle &triangle : triangles) {
            const AABB aabb = triangle.getAABB();
            min = Vec3::min(min, aabb.getMin());
            max = Vec3::max(max, aabb.getMax());
        }

        return {min, max};
    }

    /**
     * Get the area of the triangle.
     * @return The area of the triangle.
     */
    [[nodiscard]] float getArea() const {
        return area;
    }

    /**
     * Check for intersection with a ray.
     * @param ray The ray to check for intersection.
     * @return The intersection result.
     */
    [[nodiscard]] RayTriangleIntersection intersect(const Ray& ray) const {
        RayTriangleIntersection intersection;
        intersection.intersectionExists = false;

        const Vec3& v0 = m_c[0];
        const Vec3& v1 = m_c[1];
        const Vec3& v2 = m_c[2];

        // Möller-Trumbore intersection algorithm
        const Vec3 edge1 = v1 - v0;
        const Vec3 edge2 = v2 - v0;
        const Vec3 h = Vec3::cross(ray.direction(), edge2);
        const float a = Vec3::dot(edge1, h);

        if (a > -1e-6 && a < 1e-6) {
            return intersection; ///< Ray is parallel to the triangle.
        }

        const float f = 1.0f / a;
        const Vec3 s = ray.origin() - v0;
        const float u = f * Vec3::dot(s, h);

        if (u < 0.0 || u > 1.0) {
            return intersection; ///< Intersection out of bounds.
        }

        const Vec3 q = Vec3::cross(s, edge1);
        const float v = f * Vec3::dot(ray.direction(), q);

        if (v < 0.0 || u + v > 1.0) {
            return intersection; ///< Intersection out of bounds.
        }

        // Calculate t to find the intersection point.
        const float t = f * Vec3::dot(edge2, q);

        if (t > 1e-6) { ///< Valid intersection with the ray.
            intersection.intersectionExists = true;
            intersection.t = t;
            intersection.intersection = ray.origin() + ray.direction() * t;
            intersection.normal = Vec3::cross(edge1, edge2).normalized();
            intersection.u = u;
            intersection.v = v;
        }

        return intersection;
    }

    /**
     * Get the centroid of the triangle.
     * @return The centroid of the triangle.
     */
    [[nodiscard]] Vec3 getCentroid() const {
        return (m_c[0] + m_c[1] + m_c[2]) / 3.f;
    }
};

#endif // TRIANGLE_H