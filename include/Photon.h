#ifndef PHOTON_H
#define PHOTON_H

#include "Material.h"
#include "Vec3.h"

// -------------------------------------------
// Photon Structure
// -------------------------------------------

/**
 * Structure representing a photon.
 */
struct Photon {
    Vec3 position;              ///< Position of the photon.
    Vec3 direction;             ///< Direction of the photon.
    Vec3 color;                 ///< Color of the photon.
    Vec3 debugColor;            ///< Debug color of the photon.

    /**
     * Gets the position of the photon.
     * @return The position of the photon.
     */
    [[nodiscard]] Vec3 getPosition() const { return position; }
};

#endif //PHOTON_H