#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "Vec3.h"
#include "Ray.h"
#include "Plane.h"

struct RayTriangleIntersection : RayIntersection {
    float w0{},w1{},w2{};
    unsigned int tIndex{};
};

class Triangle {
private:
    Vec3 m_c[3] , m_normal;
    float area;
public:
    Triangle() {}

    Triangle(Vec3 const &c0, Vec3 const &c1, Vec3 const &c2) {
        m_c[0] = c0;
        m_c[1] = c1;
        m_c[2] = c2;
        updateAreaAndNormal();
    }

    void updateAreaAndNormal() {
        Vec3 nNotNormalized = Vec3::cross(m_c[1] - m_c[0], m_c[2] - m_c[0]);
        float norm = nNotNormalized.length();
        m_normal = nNotNormalized / norm;
        area = norm / 2.f;
    }

    void setC0(Vec3 const &c0) { m_c[0] = c0; } // remember to update the area and normal afterwards!
    void setC1(Vec3 const &c1) { m_c[1] = c1; } // remember to update the area and normal afterwards!
    void setC2(Vec3 const &c2) { m_c[2] = c2; } // remember to update the area and normal afterwards!

    [[nodiscard]] Vec3 const &normal() const {
        return m_normal;
    }

    [[nodiscard]] RayTriangleIntersection intersect(const Ray& ray) const {
        RayTriangleIntersection intersection;
        intersection.intersectionExists = false;

        const Vec3& v0 = m_c[0];
        const Vec3& v1 = m_c[1];
        const Vec3& v2 = m_c[2];

        // Algorithme d'intersection de Möller-Trumbore
        const Vec3 edge1 = v1 - v0;
        const Vec3 edge2 = v2 - v0;
        const Vec3 h = Vec3::cross(ray.direction(), edge2);
        const float a = Vec3::dot(edge1, h);

        if (a > -1e-6 && a < 1e-6) {
            return intersection; // Ce rayon est parallèle à ce triangle.
        }

        const float f = 1.0f / a;
        const Vec3 s = ray.origin() - v0;
        const float u = f * Vec3::dot(s, h);

        if (u < 0.0 || u > 1.0) {
            return intersection;
        }

        const Vec3 q = Vec3::cross(s, edge1);
        const float v = f * Vec3::dot(ray.direction(), q);

        if (v < 0.0 || u + v > 1.0) {
            return intersection;
        }

        // À ce stade, nous pouvons calculer t pour savoir où se trouve le point d'intersection sur la ligne.
        const float t = f * Vec3::dot(edge2, q);

        if (t > 1e-6) { // intersection du rayon
            intersection.intersectionExists = true;
            intersection.t = t;
            intersection.intersection = ray.origin() + ray.direction() * t;
            intersection.normal = Vec3::cross(edge1, edge2).normalized();
        }

        return intersection;
    }
};
#endif
