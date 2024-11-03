#ifndef AABB_H
#define AABB_H

#include "Vec3.h"
#include <GL/gl.h>

class AABB {
private:
    Vec3 min;
    Vec3 max;

public:
    AABB() : min(0.0f), max(0.0f) {}
    AABB(const Vec3& _min, const Vec3& _max) : min(_min), max(_max) {}
    AABB(const AABB& aabb) = default;

    [[nodiscard]] Vec3 getMin() const { return min; }
    void setMin(const Vec3& _min) { min = _min; }

    [[nodiscard]] Vec3 getMax() const { return max; }
    void setMax(const Vec3& _max) { max = _max; }

    [[nodiscard]] Vec3 center() const {
        return (min + max) * 0.5f;
    }

    void draw() const {
        glLineWidth(1.0f);
        glBegin(GL_LINE_STRIP);

        // Bottom face
        glVertex3f(min[0], min[1], min[2]);
        glVertex3f(max[0], min[1], min[2]);
        glVertex3f(max[0], min[1], max[2]);
        glVertex3f(min[0], min[1], max[2]);
        glVertex3f(min[0], min[1], min[2]);

        // Vertical edges
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
    }
};

#endif //AABB_H