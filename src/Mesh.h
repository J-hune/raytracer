#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include "Vec3.h"
#include "Ray.h"
#include "Triangle.h"
#include "Material.h"

// -------------------------------------------
// Basic Mesh class
// -------------------------------------------

struct MeshVertex {
    Vec3 position; // position
    Vec3 normal;   // normal
    float u, v;    // UV coordinates

    MeshVertex() : u(0), v(0) {}
    MeshVertex(const Vec3& _p, const Vec3& _n) : position(_p), normal(_n), u(0), v(0) {}
    MeshVertex(const Vec3& _p, const Vec3& _n, const float _u, const float _v) : position(_p), normal(_n), u(_u), v(_v) {}
    MeshVertex(const MeshVertex& vertex) = default;
    MeshVertex& operator=(const MeshVertex& vertex) = default;
};

struct MeshTriangle {
    unsigned int v[3]; // vertex indices

    MeshTriangle() : v{0, 0, 0} {}
    MeshTriangle(const unsigned int v0, const unsigned int v1, const unsigned int v2) : v{v0, v1, v2} {}
    MeshTriangle(const MeshTriangle& t) = default;
    MeshTriangle& operator=(const MeshTriangle& t) = default;
    unsigned int& operator[](const unsigned int iv) { return v[iv]; }
    unsigned int operator[](const unsigned int iv) const { return v[iv]; }
};

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

    void loadOFF(const std::string& filename);
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

    [[nodiscard]] virtual RayTriangleIntersection intersect(const Ray& ray) const;

    virtual void buildArrays();

private:
    void buildPositionsArray();
    void buildNormalsArray();
    void buildUVsArray();
    void buildTrianglesArray();
};

#endif // MESH_H