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

        // Récupération des données du carré
        Vec3 m_bottom_left = vertices[0].position;
        Vec3 m_right_vector = vertices[1].position - vertices[0].position;
        Vec3 m_up_vector = vertices[3].position - vertices[0].position;
        Vec3 m_normal = Vec3::cross(m_right_vector, m_up_vector);
        m_normal.normalize();

        // Calcul de l'intersection avec le plan du carré (si le rayon est parallèle, il n'y a pas d'intersection)
        float denominator = Vec3::dot(m_normal, ray.direction());
        float numerator = Vec3::dot(m_bottom_left - ray.origin(), m_normal);

        if (std::fabs(denominator) < 1e-6) {
            return intersection; // Pas d'intersection car le rayon est parallèle au plan
        }

        // Calcul de la distance t pour laquelle le rayon intersecte le plan du carré
        float t = numerator / denominator;

        // Si t < 0, l'intersection est derrière l'origine du rayon
        if (t < 0) {
            return intersection; // Pas d'intersection
        }

        // Calculer le point d'intersection avec le plan
        const Vec3 PointIntersection = ray.origin() + t * ray.direction();

        // Transformer P dans l'espace local du carré
        const Vec3 localP = PointIntersection - m_bottom_left;

        // Calculer les projections sur les axes droit (right_vector) et haut (up_vector)
        const float u = Vec3::dot(localP, m_right_vector) / m_right_vector.squareLength();
        const float v = Vec3::dot(localP, m_up_vector) / m_up_vector.squareLength();

        // Vérifier si l'intersection est à l'intérieur des limites du carré (u, v doivent être entre 0 et 1)
        if (u >= 0 && u <= 1 && v >= 0 && v <= 1) {
            intersection.intersectionExists = true;
            intersection.t = t;
            intersection.u = u;
            intersection.v = v;
            intersection.intersection = PointIntersection;
            intersection.normal = m_normal;  // Normale du plan
        }

        return intersection;
    }
};
#endif // SQUARE_H
