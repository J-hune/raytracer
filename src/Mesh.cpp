#include "../include/Mesh.h"
#include "../include/Triangle.h"

#include <cfloat>
#include <cmath>
#include <fstream>
#include <GL/gl.h>

std::vector<Triangle> Mesh::getTriangles() const {
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

    computeAABB();
}

void Mesh::applyTransformationMatrix(const Mat3& transform) {
    for (auto& vertex : vertices) {
        vertex.position = transform * vertex.position;
    }
    computeAABB();
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
    if (triangles_array.empty()) return;
    const GLfloat material_color[4] = {
        material.diffuse_material[0],
        material.diffuse_material[1],
        material.diffuse_material[2],
        1.0
    };

    if (textureID > 0) {
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, textureID);
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
        glTexCoordPointer(2, GL_FLOAT, 2 * sizeof(float), uvs_array.data());
    }

    const GLfloat material_specular[4] = {
        material.specular_material[0],
        material.specular_material[1],
        material.specular_material[2],
        1.0
    };

    const GLfloat material_ambient[4] = {
        material.ambient_material[0],
        material.ambient_material[1],
        material.ambient_material[2],
        1.0
    };

    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_color);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, static_cast<GLfloat>(material.shininess));

    glEnableClientState(GL_VERTEX_ARRAY) ;
    glEnableClientState (GL_NORMAL_ARRAY);
    glNormalPointer(GL_FLOAT, 3 * sizeof(float), normals_array.data());
    glVertexPointer(3, GL_FLOAT, 3 * sizeof(float), positions_array.data());
    glDrawElements(GL_TRIANGLES, static_cast<int>(triangles_array.size()), GL_UNSIGNED_INT, triangles_array.data());

    if (textureID > 0) {
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        glDisable(GL_TEXTURE_2D);
    }
}

MeshVertex Mesh::getRandomPointOnSurface(std::mt19937 &rng) const {
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

    const Vec3 position = v0.position * (1.0f - u - v) + v1.position * u + v2.position * v;
    const Vec3 normal = v0.normal * (1.0f - u - v) + v1.normal * u + v2.normal * v;

    return {position, normal};
}

void Mesh::computeAABB() {
    Vec3 min = vertices[0].position;
    Vec3 max = vertices[0].position;

    for (const auto& vertex : vertices) {
        for (int i = 0; i < 3; ++i) {
            if (vertex.position[i] < min[i]) min[i] = vertex.position[i];
            if (vertex.position[i] > max[i]) max[i] = vertex.position[i];
        }
    }

    aabb = AABB(min, max);
}

[[nodiscard]] bool Mesh::intersectAABB(const Ray& ray) const {
    float tmin = (aabb.getMin()[0] - ray.origin()[0]) / ray.direction()[0];
    float tmax = (aabb.getMax()[0] - ray.origin()[0]) / ray.direction()[0];

    if (tmin > tmax) std::swap(tmin, tmax);

    for (int i = 1; i < 3; ++i) {
        float t1 = (aabb.getMin()[i] - ray.origin()[i]) / ray.direction()[i];
        float t2 = (aabb.getMax()[i] - ray.origin()[i]) / ray.direction()[i];

        if (t1 > t2) std::swap(t1, t2);

        tmin = std::max(tmin, t1);
        tmax = std::min(tmax, t2);

        if (tmin > tmax) return false;
    }

    return true;
}

RayIntersection Mesh::intersect(const Ray& ray) const {
    RayIntersection closestIntersection;
    closestIntersection.intersectionExists = false;
    closestIntersection.t = FLT_MAX;

    for (const auto& triangle : triangles) {
        Triangle tri(vertices[triangle.v[0]].position, vertices[triangle.v[1]].position, vertices[triangle.v[2]].position);
        const RayIntersection intersection = tri.intersect(ray);

        if (intersection.intersectionExists && intersection.t < closestIntersection.t) {
            closestIntersection.intersectionExists = true;
            closestIntersection.t = intersection.t;
            closestIntersection.intersection = intersection.intersection;
            closestIntersection.normal = intersection.normal;
        }
    }

    return closestIntersection;
}

// Load mesh from OFF file
void Mesh::loadOFF(const std::string &filename) {
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
        for (unsigned int &j: triangles[i].v) {
            in >> j;
        }
    }
    in.close();

    // Compute AABB of the mesh
    computeAABB();
}

void Mesh::loadTexture(const std::string &filename) {
    // Open PPM file in binary mode
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Failed to open PPM file: " << filename << std::endl;
        return;
    }

    std::string header;
    int maxColor = 0;

    // Read header and skip comments
    file >> header;
    while (file.peek() == '#' || file.peek() == '\n') file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    file >> textureWidth >> textureHeight >> maxColor;
    file.ignore();

    if (header != "P6" || maxColor != 255) {
        std::cerr << "Invalid PPM file: " << filename << std::endl;
        return;
    }

    texture_array.resize(textureWidth * textureHeight * 3);
    if (!file.read(reinterpret_cast<char*>(texture_array.data()), static_cast<std::streamsize>(texture_array.size()))) {
        std::cerr << "Failed to read pixel data from file: " << filename << std::endl;
        return;
    }

    // Generate OpenGL texture
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, textureWidth, textureHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, texture_array.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
}


Vec3 Mesh::sampleTexture(float u, float v) const {
    // Clamp UV coordinates to [0, 1]
    u = std::clamp(u, 0.0f, 1.0f);
    v = std::clamp(v, 0.0f, 1.0f);

    // Convert UV coordinates to pixel coordinates
    const int x = static_cast<int>(u * (static_cast<float>(textureWidth) - 1.0f));
    const int y = static_cast<int>(v * (static_cast<float>(textureHeight) - 1.0f));

    // Calculate the index in the texture data
    const int index = (y * textureWidth + x) * 3; // Assuming RGB format

    // Get the color values (ensure index is within bounds)
    if (index < 0 || index + 2 >= static_cast<int>(texture_array.size())) {
        return {-1.0f, -1.0f, -1.0f};
    }

    const float r = static_cast<float>(texture_array[index]) / 255.0f;
    const float g = static_cast<float>(texture_array[index + 1]) / 255.0f;
    const float b = static_cast<float>(texture_array[index + 2]) / 255.0f;

    return {r, g, b};
}
