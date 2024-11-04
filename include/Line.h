#ifndef LINE_H
#define LINE_H

#include "Vec3.h"
#include <iostream>

/**
 * Class representing a line in 3D space.
 */
class Line {
private:
    Vec3 m_origin; ///< Origin point of the line.
    Vec3 m_direction; ///< Direction vector of the line.
public:
    Line() = default;
    Line(const Vec3 &o, const Vec3 &d) {
        m_origin = o;
        m_direction = d;
        m_direction.normalize();
    }

    /**
     * Gets the origin point of the line.
     * @return Reference to the origin point.
     */
    Vec3 &origin() { return m_origin; }

    /**
     * Gets the origin point of the line (const version).
     * @return Constant reference to the origin point.
     */
    [[nodiscard]] const Vec3 &origin() const { return m_origin; }

    /**
     * Gets the direction vector of the line (const version).
     * @return Constant reference to the direction vector.
     */
    [[nodiscard]] const Vec3 &direction() const { return m_direction; }

    /**
     * Gets the direction vector of the line.
     * @return Reference to the direction vector.
     */
    Vec3 &direction() { return m_direction; }
};

/**
 * Overloaded output stream operator for Line.
 * @param s The output stream.
 * @param l The line to output.
 * @return Reference to the output stream.
 */
static inline std::ostream &operator<<(std::ostream &s, const Line &l) {
    s << l.origin() << " " << l.direction();
    return s;
}

#endif
