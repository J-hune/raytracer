#ifndef MESH_H
#define MESH_H

#include <random>
#include <vector>
#include <string>
#include <stdexcept>

#include "Vec3.h"
#include "Material.h"
#include "AABB.h"
#include "Ray.h"
#include "Triangle.h"

// -------------------------------------------
// MeshVertex Structure
// -------------------------------------------
struct MeshVertex {
    Vec3 position;
    Vec3 normal;
    float u, v;

    MeshVertex() : u(0), v(0) {}
    MeshVertex(const Vec3& _p, const Vec3& _n) : position(_p), normal(_n), u(0), v(0) {}
    MeshVertex(const Vec3& _p, const Vec3& _n, const float _u, const float _v) : position(_p), normal(_n), u(_u), v(_v) {}
    MeshVertex(const MeshVertex& vertex) = default;
    MeshVertex& operator=(const MeshVertex& vertex) = default;
};

// -------------------------------------------
// MeshTriangle Structure
// -------------------------------------------
struct MeshTriangle {
    unsigned int v[3]; // vertex indices

    MeshTriangle() : v{0, 0, 0} {}
    MeshTriangle(const unsigned int v0, const unsigned int v1, const unsigned int v2) : v{v0, v1, v2} {}
    MeshTriangle(const MeshTriangle& t) = default;
    MeshTriangle& operator=(const MeshTriangle& t) = default;
    unsigned int& operator[](const unsigned int iv) { return v[iv]; }
    unsigned int operator[](const unsigned int iv) const { return v[iv]; }
};

// -------------------------------------------
// Mesh Class
// -------------------------------------------
class Mesh {
protected:
    std::vector<float> positions_array;
    std::vector<float> normals_array;
    std::vector<float> uvs_array;
    std::vector<unsigned int> triangles_array;

public:
    virtual ~Mesh() = default;

    std::vector<MeshVertex> vertices;
    std::vector<MeshTriangle> triangles;
    Material material;
    AABB aabb;

    [[nodiscard]] Vec3 getPosition() const {
        return aabb.center();
    }

    [[nodiscard]] std::vector<Triangle> getTriangles() const {
        std::vector<Triangle> meshTriangles;
        meshTriangles.reserve(triangles.size());

        // Convert mesh triangles to triangles
        for (const MeshTriangle& triangle : triangles) {
            const MeshVertex& v0 = vertices[triangle[0]];
            const MeshVertex& v1 = vertices[triangle[1]];
            const MeshVertex& v2 = vertices[triangle[2]];

            meshTriangles.emplace_back(v0.position, v1.position, v2.position);
        }

        return meshTriangles;
    }

    void loadOFF(const std::string& filename);

    MeshVertex getRandomPointOnSurface(std::mt19937 &rng) const {
        if (triangles.empty() || vertices.empty()) {
            throw std::runtime_error("Mesh has no triangles or vertices.");
        }

        std::uniform_int_distribution<size_t> triangleDist(0, triangles.size() - 1);
        const MeshTriangle& triangle = triangles[triangleDist(rng)];

        const MeshVertex& v0 = vertices[triangle[0]];
        const MeshVertex& v1 = vertices[triangle[1]];
        const MeshVertex& v2 = vertices[triangle[2]];

        std::uniform_real_distribution<float> dist(0.0f, 1.0f);
        float u = dist(rng);
        float v = dist(rng);

        if (u + v > 1.0f) {
            u = 1.0f - u;
            v = 1.0f - v;
        }

        Vec3 position = v0.position * (1.0f - u - v) + v1.position * u + v2.position * v;
        Vec3 normal = v0.normal * (1.0f - u - v) + v1.normal * u + v2.normal * v;

        return MeshVertex(position, normal);
    }

    void recomputeNormals();
    void centerAndScaleToUnit();
    void scaleUnit();

    void translate(const Vec3& translation);
    void applyTransformationMatrix(const Mat3& transform);
    void scale(const Vec3& scale);
    void rotateX(float angle);
    void rotateY(float angle);
    void rotateZ(float angle);
    void draw() const;

    virtual void computeAABB();
    [[nodiscard]] virtual bool intersectAABB(const Ray& ray) const;
    [[nodiscard]] virtual RayIntersection intersect(const Ray& ray) const;

    virtual void buildArrays();

private:
    void buildPositionsArray();
    void buildNormalsArray();
    void buildUVsArray();
    void buildTrianglesArray();
};

#endif // MESH_H