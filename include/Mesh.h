#ifndef MESH_H
#define MESH_H

#include <random>
#include <vector>
#include <string>

#include "Vec3.h"
#include "Material.h"
#include "AABB.h"
#include "Ray.h"
#include "Triangle.h"

// -------------------------------------------
// MeshVertex Structure
// -------------------------------------------

/**
 * Structure representing a vertex in a mesh.
 */
struct MeshVertex {
    Vec3 position;  ///< Position of the vertex.
    Vec3 normal;    ///< Normal vector at the vertex.
    float u, v;     ///< Texture coordinates.

    MeshVertex() : u(0), v(0) {}
    MeshVertex(const Vec3& _p, const Vec3& _n) : position(_p), normal(_n), u(0), v(0) {}
    MeshVertex(const Vec3& _p, const Vec3& _n, const float _u, const float _v) : position(_p), normal(_n), u(_u), v(_v) {}
    MeshVertex(const MeshVertex& vertex) = default;

    /**
     * Assignment operator for MeshVertex.
     * @param vertex The vertex to assign.
     * @return Reference to the assigned vertex.
     */
    MeshVertex& operator=(const MeshVertex& vertex) = default;
};

// -------------------------------------------
// MeshTriangle Structure
// -------------------------------------------

/**
 * Structure representing a triangle in a mesh.
 */
struct MeshTriangle {
    unsigned int v[3]; ///< Indices of the vertices forming the triangle.

    MeshTriangle() : v{0, 0, 0} {}
    MeshTriangle(const unsigned int v0, const unsigned int v1, const unsigned int v2) : v{v0, v1, v2} {}
    MeshTriangle(const MeshTriangle& t) = default;

    /**
     * Assignment operator for MeshTriangle.
     * @param t The triangle to assign.
     * @return Reference to the assigned triangle.
     */
    MeshTriangle& operator=(const MeshTriangle& t) = default;

    /**
     * Access operator for vertex indices.
     * @param iv Index of the vertex to access.
     * @return Reference to the vertex index.
     */
    unsigned int& operator[](const unsigned int iv) { return v[iv]; }

    /**
     * Access operator for vertex indices (const version).
     * @param iv Index of the vertex to access.
     * @return The vertex index.
     */
    unsigned int operator[](const unsigned int iv) const { return v[iv]; }
};

// -------------------------------------------
// Mesh Class
// -------------------------------------------

/**
 * Class representing a mesh.
 */
class Mesh {
protected:
    std::vector<float> positions_array;         ///< Array of vertex positions.
    std::vector<float> normals_array;           ///< Array of vertex normals.
    std::vector<float> uvs_array;               ///< Array of texture coordinates.
    std::vector<unsigned int> triangles_array;  ///< Array of triangle vertex indices.
    std::vector<uint8_t> texture_array;         ///< Array of texture colors.

public:
    virtual ~Mesh() = default;

    std::vector<MeshVertex> vertices;       ///< List of vertices in the mesh.
    std::vector<MeshTriangle> triangles;    ///< List of triangles in the mesh.
    Material material;                      ///< Material properties of the mesh.
    AABB aabb;                              ///< Axis-aligned bounding box of the mesh.
    GLuint textureID = 0;                   ///< OpenGL texture ID.
    int textureWidth = 0;                   ///< Width of the texture.
    int textureHeight = 0;                  ///< Height of the texture.

    /**
     * Gets the position of the mesh.
     * @return The center of the axis-aligned bounding box.
     */
    [[nodiscard]] Vec3 getPosition() const { return aabb.center(); }

    /**
     * Gets the list of triangles in the mesh.
     * @return The list of triangles.
     */
    [[nodiscard]] std::vector<Triangle> getTriangles() const;

    /**
     * Loads a mesh from an OFF file.
     * @param filename The name of the OFF file.
     */
    void loadOFF(const std::string &filename);

    /**
     * Loads a texture from a PPM file.
     * @param filename The name of the PPM file.
     * @return The OpenGL texture ID.
     */
    void loadTexture(const std::string &filename);

    /**
     * Samples the texture at the specified UV coordinates.
     * @param u The U coordinate.
     * @param v The V coordinate.
     * @return The color of the texture at the specified coordinates.
     */
    [[nodiscard]] Vec3 sampleTexture(float u, float v) const;

    /**
     * Gets a random point on the surface of the mesh.
     * @param rng A random number generator.
     * @return A random point on the surface.
     */
    MeshVertex getRandomPointOnSurface(std::mt19937 &rng) const;

    /**
     * Recomputes the normals of the mesh.
     */
    void recomputeNormals();

    /**
     * Centers the mesh and scales it to unit size.
     */
    void centerAndScaleToUnit();

    /**
     * Translates the mesh.
     * @param translation The translation vector.
     */
    void translate(const Vec3& translation);

    /**
     * Applies a transformation matrix to the mesh.
     * @param transform The transformation matrix.
     */
    void applyTransformationMatrix(const Mat3& transform);

    /**
     * Scales the mesh.
     * @param scale The scale vector.
     */
    void scale(const Vec3& scale);

    /**
     * Rotates the mesh around the X axis.
     * @param angle The rotation angle in radians.
     */
    void rotateX(float angle);

    /**
     * Rotates the mesh around the Y axis.
     * @param angle The rotation angle in radians.
     */
    void rotateY(float angle);

    /**
     * Rotates the mesh around the Z axis.
     * @param angle The rotation angle in radians.
     */
    void rotateZ(float angle);

    /**
     * Builds the arrays for rendering the mesh.
     */
    virtual void buildArrays();

    /**
     * Draws the mesh.
     */
    void draw() const;

    /**
     * Computes the axis-aligned bounding box of the mesh.
     */
    virtual void computeAABB();

    /**
     * Checks if a ray intersects the axis-aligned bounding box of the mesh.
     * @param ray The ray to check.
     * @return True if the ray intersects the bounding box, false otherwise.
     */
    [[nodiscard]] virtual bool intersectAABB(const Ray& ray) const;

    /**
     * Checks if a ray intersects the mesh.
     * @param ray The ray to check.
     * @return The intersection information.
     */
    [[nodiscard]] virtual RayIntersection intersect(const Ray& ray) const;

private:
    /**
     * Builds the array of vertex positions.
     */
    void buildPositionsArray();

    /**
     * Builds the array of vertex normals.
     */
    void buildNormalsArray();

    /**
     * Builds the array of texture coordinates.
     */
    void buildUVsArray();

    /**
     * Builds the array of triangle vertex indices.
     */
    void buildTrianglesArray();
};

#endif // MESH_H