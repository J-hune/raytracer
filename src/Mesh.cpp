#include "Mesh.h"

#include <cfloat>
#include <cmath>
#include <fstream>
#include <GL/gl.h>

// Load mesh from OFF file
void Mesh::loadOFF(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string offString;
    unsigned int sizeV, sizeT, tmp;
    in >> offString >> sizeV >> sizeT >> tmp;

    vertices.resize(sizeV);
    triangles.resize(sizeT);

    for (unsigned int i = 0; i < sizeV; i++) {
        in >> vertices[i].position;
    }

    for (unsigned int i = 0; i < sizeT; i++) {
        unsigned int s;
        in >> s;
        for (unsigned int & j : triangles[i].v) {
            in >> j;
        }
    }
    in.close();
}

// Recompute normals for each vertex
void Mesh::recomputeNormals() {
    for (auto &vertice: vertices) {
        vertice.normal = Vec3(0.0, 0.0, 0.0);
    }

    for (const auto &triangle: triangles) {
        Vec3 e01 = vertices[triangle.v[1]].position - vertices[triangle.v[0]].position;
        Vec3 e02 = vertices[triangle.v[2]].position - vertices[triangle.v[0]].position;
        const Vec3 n = Vec3::cross(e01, e02).normalized();
        for (const unsigned int j: triangle.v) {
            vertices[j].normal += n;
        }
    }

    for (auto &vertice: vertices)
        vertice.normal.normalize();
}

// Center and scale mesh to unit size
void Mesh::centerAndScaleToUnit() {
    Vec3 centroid(0, 0, 0);
    for (const auto& vertex : vertices) {
        centroid += vertex.position;
    }
    centroid /= static_cast<float>(vertices.size());

    float maxD = 0.0f;
    for (const auto& vertex : vertices) {
        maxD = std::max(maxD, (vertex.position - centroid).length());
    }

    for (auto& vertex : vertices) {
        vertex.position = (vertex.position - centroid) / maxD;
    }
}

// Build arrays for OpenGL rendering
void Mesh::buildArrays() {
    recomputeNormals();
    buildPositionsArray();
    buildNormalsArray();
    buildUVsArray();
    buildTrianglesArray();
}

// Implementations for building arrays (optimized)
void Mesh::buildPositionsArray() {
    positions_array.resize(3 * vertices.size());
    for (unsigned int v = 0; v < vertices.size(); ++v) {
        positions_array[3 * v + 0] = vertices[v].position[0];
        positions_array[3 * v + 1] = vertices[v].position[1];
        positions_array[3 * v + 2] = vertices[v].position[2];
    }
}

void Mesh::buildNormalsArray() {
    normals_array.resize(3 * vertices.size());
    for (unsigned int v = 0; v < vertices.size(); ++v) {
        normals_array[3 * v + 0] = vertices[v].normal[0];
        normals_array[3 * v + 1] = vertices[v].normal[1];
        normals_array[3 * v + 2] = vertices[v].normal[2];
    }
}

void Mesh::buildUVsArray() {
    uvs_array.resize(2 * vertices.size());
    for (unsigned int vert = 0; vert < vertices.size(); ++vert) {
        uvs_array[2 * vert + 0] = vertices[vert].u;
        uvs_array[2 * vert + 1] = vertices[vert].v;
    }
}

void Mesh::buildTrianglesArray() {
    triangles_array.resize(3 * triangles.size());
    for (unsigned int t = 0; t < triangles.size(); ++t) {
        triangles_array[3 * t + 0] = triangles[t].v[0];
        triangles_array[3 * t + 1] = triangles[t].v[1];
        triangles_array[3 * t + 2] = triangles[t].v[2];
    }
}

// Transformation methods
void Mesh::translate(const Vec3& translation) {
    for (auto& vertex : vertices) {
        vertex.position += translation;
    }
}

void Mesh::applyTransformationMatrix(const Mat3& transform) {
    for (auto& vertex : vertices) {
        vertex.position = transform * vertex.position;
    }
}

void Mesh::scale(const Vec3 &scale) {
    const Mat3 scale_matrix(
        scale[0], 0.0f, 0.0f,
        0.0f, scale[1], 0.0f,
        0.0f, 0.0f, scale[2]
    );
    applyTransformationMatrix(scale_matrix);
}

void Mesh::rotateX(const float angle) {
    const float radians = angle * M_PIf / 180.0f;
    const Mat3 rotation(
        1.0f, 0.0f, 0.0f,
        0.0f, std::cos(radians), -std::sin(radians),
        0.0f, std::sin(radians), std::cos(radians)
    );
    applyTransformationMatrix(rotation);
}

void Mesh::rotateY(const float angle) {
    const float radians = angle * M_PIf / 180.0f;
    const Mat3 rotation(
        std::cos(radians), 0.0f, std::sin(radians),
        0.0f, 1.0f, 0.0f,
        -std::sin(radians), 0.0f, std::cos(radians)
    );
    applyTransformationMatrix(rotation);
}

void Mesh::rotateZ(const float angle) {
    const float radians = angle * M_PIf / 180.0f;
    const Mat3 rotation(
        std::cos(radians), -std::sin(radians), 0.0f,
        std::sin(radians), std::cos(radians), 0.0f,
        0.0f, 0.0f, 1.0f
    );
    applyTransformationMatrix(rotation);
}

// Draw the mesh using OpenGL
void Mesh::draw() const {
    if( triangles_array.size() == 0 ) return;
    GLfloat material_color[4] = {
        material.diffuse_material[0],
        material.diffuse_material[1],
        material.diffuse_material[2],
        1.0
    };

    GLfloat material_specular[4] = {
        material.specular_material[0],
        material.specular_material[1],
        material.specular_material[2],
        1.0
    };

    GLfloat material_ambient[4] = {
        material.ambient_material[0],
        material.ambient_material[1],
        material.ambient_material[2],
        1.0
    };

    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_color);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, material.shininess);

    glEnableClientState(GL_VERTEX_ARRAY) ;
    glEnableClientState (GL_NORMAL_ARRAY);
    glNormalPointer (GL_FLOAT, 3*sizeof (float), (GLvoid*)(normals_array.data()));
    glVertexPointer (3, GL_FLOAT, 3*sizeof (float) , (GLvoid*)(positions_array.data()));
    glDrawElements(GL_TRIANGLES, triangles_array.size(), GL_UNSIGNED_INT, (GLvoid*)(triangles_array.data()));
}

// Intersection with a ray
RayTriangleIntersection Mesh::intersect(const Ray& ray) const {
    RayTriangleIntersection closestIntersection;
    closestIntersection.intersectionExists = false;
    closestIntersection.t = FLT_MAX;

    for (const auto& triangle : triangles) {
        const Vec3& v0 = vertices[triangle.v[0]].position;
        const Vec3& v1 = vertices[triangle.v[1]].position;
        const Vec3& v2 = vertices[triangle.v[2]].position;

        // Möller-Trumbore intersection algorithm
        Vec3 edge1 = v1 - v0;
        Vec3 edge2 = v2 - v0;
        Vec3 h = Vec3::cross(ray.direction(), edge2);
        const float a = Vec3::dot(edge1, h);

        if (a > -1e-6 && a < 1e-6) {
            continue; // This ray is parallel to this triangle.
        }

        const float f = 1.0f / a;
        Vec3 s = ray.origin() - v0;
        const float u = f * Vec3::dot(s, h);

        if (u < 0.0 || u > 1.0) {
            continue;
        }

        Vec3 q = Vec3::cross(s, edge1);
        const float v = f * Vec3::dot(ray.direction(), q);

        if (v < 0.0 || u + v > 1.0) {
            continue;
        }

        // At this stage we can compute t to find out where the intersection point is on the line.
        const float t = f * Vec3::dot(edge2, q);

        if (t > 1e-6 && t < closestIntersection.t) { // ray intersection
            closestIntersection.intersectionExists = true;
            closestIntersection.t = t;
            closestIntersection.intersection = ray.origin() + ray.direction() * t;
            closestIntersection.normal = Vec3::cross(edge1, edge2).normalized();
        }
    }

    return closestIntersection;
}