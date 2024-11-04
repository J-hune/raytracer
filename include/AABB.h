#ifndef AABB_H
#define AABB_H

#include <cfloat>
#include <GL/gl.h>
#include "Vec3.h"

/**
 * Class representing an axis-aligned bounding box (AABB).
 */
class AABB {
private:
    Vec3 min;   ///< Minimum vector
    Vec3 max;   ///< Maximum vector

public:
    AABB() : min(0.0f), max(0.0f) {}
    AABB(const Vec3& _min, const Vec3& _max) : min(_min), max(_max) {}
    AABB(const AABB& aabb) = default;

    /**
     * Gets the minimum vector of the AABB.
     * @return The minimum vector.
     */
    [[nodiscard]] Vec3 getMin() const { return min; }

    /**
     * Sets the minimum vector of the AABB.
     * @param _min The minimum vector.
     */
    void setMin(const Vec3& _min) { min = _min; }

    /**
     * Sets the minimum vector of the AABB.
     * @param _min The new minimum vector.
     */
    [[nodiscard]] Vec3 getMax() const { return max; }

    /**
     * Sets the maximum vector of the AABB.
     * @param _max The new maximum vector.
     */
    void setMax(const Vec3& _max) { max = _max; }

    /**
     * Calculates the center of the AABB.
     * @return The center vector.
     */
    [[nodiscard]] Vec3 center() const { return (min + max) * 0.5f; }

    /**
     * Merges two AABB objects into one.
     * @param a The first AABB.
     * @param b The second AABB.
     * @return A new AABB that encompasses both input AABBs.
     */
    static AABB merge(const AABB& a, const AABB& b) {
        return {Vec3::min(a.min, b.min), Vec3::max(a.max, b.max)};
    }

    /**
     * Draws the AABB using OpenGL.
     * This method disables lighting and enables color material to draw the AABB as a wireframe.
     * It uses a line strip to draw the edges of the AABB.
     */
    void draw() const {
        glDisable(GL_LIGHTING);
        glEnable(GL_COLOR_MATERIAL);
        glLineWidth(1.0f);
        glBegin(GL_LINE_STRIP);
        glColor3f(0.8f, 0.8f, 0.2f);

        // Bottom face
        glVertex3f(min[0], min[1], min[2]);
        glVertex3f(max[0], min[1], min[2]);
        glVertex3f(max[0], min[1], max[2]);
        glVertex3f(min[0], min[1], max[2]);
        glVertex3f(min[0], min[1], min[2]);

        // Vertical edges and top face
        glVertex3f(min[0], max[1], min[2]);
        glVertex3f(max[0], max[1], min[2]);
        glVertex3f(max[0], min[1], min[2]);
        glVertex3f(max[0], min[1], max[2]);
        glVertex3f(max[0], max[1], max[2]);
        glVertex3f(min[0], max[1], max[2]);
        glVertex3f(min[0], min[1], max[2]);
        glVertex3f(min[0], max[1], max[2]);
        glVertex3f(min[0], max[1], min[2]);
        glVertex3f(max[0], max[1], min[2]);
        glVertex3f(max[0], max[1], max[2]);

        glEnd();
        glDisable(GL_COLOR_MATERIAL);
        glEnable(GL_LIGHTING);
    }
};

#endif // AABB_H