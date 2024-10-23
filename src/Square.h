#ifndef SQUARE_H
#define SQUARE_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySquareIntersection{
    bool intersectionExists;
    float t;
    float u,v;
    Vec3 intersection;
    Vec3 normal;
};


class Square : public Mesh {
public:
    Vec3 m_normal;
    Vec3 m_bottom_left;
    Vec3 m_right_vector;
    Vec3 m_up_vector;

    Square() : Mesh() {}
    Square(Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
           float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) : Mesh() {
        setQuad(bottomLeft, rightVector, upVector, width, height, uMin, uMax, vMin, vMax);
    }

    void setQuad( Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
                  float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) {
        m_right_vector = rightVector;
        m_up_vector = upVector;
        m_normal = Vec3::cross(rightVector , upVector);
        m_bottom_left = bottomLeft;

        m_normal.normalize();
        m_right_vector.normalize();
        m_up_vector.normalize();

        m_right_vector = m_right_vector*width;
        m_up_vector = m_up_vector*height;

        vertices.clear();
        vertices.resize(4);
        vertices[0].position = bottomLeft;                                      vertices[0].u = uMin; vertices[0].v = vMin;
        vertices[1].position = bottomLeft + m_right_vector;                     vertices[1].u = uMax; vertices[1].v = vMin;
        vertices[2].position = bottomLeft + m_right_vector + m_up_vector;       vertices[2].u = uMax; vertices[2].v = vMax;
        vertices[3].position = bottomLeft + m_up_vector;                        vertices[3].u = uMin; vertices[3].v = vMax;
        vertices[0].normal = vertices[1].normal = vertices[2].normal = vertices[3].normal = m_normal;
        triangles.clear();
        triangles.resize(2);
        triangles[0][0] = 0;
        triangles[0][1] = 1;
        triangles[0][2] = 2;
        triangles[1][0] = 0;
        triangles[1][1] = 2;
        triangles[1][2] = 3;
    }

    RaySquareIntersection intersect(const Ray &ray) const {
    RaySquareIntersection intersection;
    intersection.intersectionExists = false;

    // Retrieve square data
    Vec3 m_bottom_left = vertices[0].position;
    Vec3 m_right_vector = vertices[1].position - vertices[0].position;
    Vec3 m_up_vector = vertices[3].position - vertices[0].position;
    Vec3 m_normal = Vec3::cross(m_right_vector, m_up_vector);
    m_normal.normalize();

    // Calculate intersection with the square's plane (if the ray is parallel, there is no intersection)
    const float denominator = Vec3::dot(m_normal, ray.direction());
    const float numerator = Vec3::dot(m_bottom_left - ray.origin(), m_normal);

    if (std::fabs(denominator) < 1e-6) {
        return intersection; // No intersection because the ray is parallel to the plane
    }

    // Calculate the distance t at which the ray intersects the plane
    const float t = numerator / denominator;

    // If t < 0, the intersection is behind the ray's origin
    if (t < 0) {
        return intersection; // No intersection
    }

    // Calculate the intersection point with the plane
    const Vec3 PointIntersection = ray.origin() + t * ray.direction();

    // Transform P into the square's local space
    const Vec3 localP = PointIntersection - m_bottom_left;

    // Calculate projections on the right (right_vector) and up (up_vector) axes
    const float u = Vec3::dot(localP, m_right_vector) / m_right_vector.squareLength();
    const float v = Vec3::dot(localP, m_up_vector) / m_up_vector.squareLength();

    // Check if the intersection is within the square's bounds (u, v must be between 0 and 1)
    if (u >= 0 && u <= 1 && v >= 0 && v <= 1) {
        intersection.intersectionExists = true;
        intersection.t = t;
        intersection.u = u;
        intersection.v = v;
        intersection.intersection = PointIntersection;
        intersection.normal = m_normal;  // Plane normal
    }

    return intersection;
}
};
#endif // SQUARE_H
