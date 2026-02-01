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
    Vec3 position;      ///< Position of the photon.
    Vec3 direction;     ///< Direction of the photon.
    Vec3 color;         ///< Color of the photon.

    /// True once the photon has been reflected by a mirror or refracted by a transparent surface.
    /// The caustics map only stores photons carrying this flag: a caustic is by definition light
    /// concentrated by a specular or refractive element before reaching a diffuse surface.
    bool hasSpecularBounce = false;

    /**
     * Gets the position of the photon.
     * @return The position of the photon.
     */
    [[nodiscard]] Vec3 getPosition() const { return position; }
};

#endif //PHOTON_H