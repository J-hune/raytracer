#ifndef MATERIAL_H
#define MATERIAL_H

#include "Vec3.h"

/**
 * Enum representing the type of material.
 */
enum MaterialType {
    Material_Diffuse_Blinn_Phong,   ///< Diffuse Blinn-Phong material.
    Material_Glass,                 ///< Glass material.
    Material_Mirror                 ///< Mirror material.
};

/**
 * Structure representing material properties.
 */
struct Material {
    Vec3 ambient_material;      ///< Ambient material properties.
    Vec3 diffuse_material;      ///< Diffuse material properties.
    Vec3 specular_material;     ///< Specular material properties.
    double shininess{};         ///< Shininess coefficient.

    float index_medium;         ///< Index of refraction for the medium.
    float transparency;         ///< Transparency factor.

    MaterialType type;          ///< Type of the material.

    Material() {
        type = Material_Diffuse_Blinn_Phong;
        transparency = 0.0;
        index_medium = 1.0;
        ambient_material = Vec3(0., 0., 0.);
    }
};

#endif // MATERIAL_H