#ifndef SETTINGS_H
#define SETTINGS_H

#pragma once

// -------------------------------------------
// TypeOfFloor Enum
// -------------------------------------------

/**
 * Enum representing the type of floor pattern.
 */
enum TypeOfFloor {
    CHECKERBOARD, ///< Checkerboard pattern.
    PLAIN ///< Plain pattern.
};

// -------------------------------------------
// Settings Class
// -------------------------------------------

/**
 * Class representing the settings for the application.
 * This class follows the Singleton design pattern.
 */
class Settings {
public:
    /**
     * Gets the singleton instance of the Settings class.
     * @return The singleton instance of the Settings class.
     */
    static Settings& getInstance() {
        static Settings instance; ///< Singleton instance of the Settings class.
        return instance;
    }

    // Delete copy constructor and assignment operator to prevent copying
    Settings(const Settings&) = delete;
    Settings& operator=(const Settings&) = delete;

    mutable int width = 480;                    ///< Width of the window in pixels.
    mutable int height = 480;                   ///< Height of the window in pixels.
    mutable int samples = 100;                  ///< Number of samples per pixel for antialiasing to reduce noise.
    mutable int shadowRays = 16;                ///< Number of rays used to calculate shadows for soft shadow effects.
    int globalPhotons = 20000;                  ///< Number of photons to emit for global illumination.
    int causticsPhotons = 20000;                ///< Number of photons to emit for caustics effects.
    bool directIllumination = true;             ///< Enable or disable direct illumination using the Phong model.
    bool reflections = true;                    ///< Enable or disable reflections for mirror-like surfaces.
    bool refractions = true;                    ///< Enable or disable refractions for glass-like surfaces.
    bool caustics = true;                       ///< Enable or disable caustics effects using photon mapping.
    mutable short drawDebugPhotons = 0;         ///< 0: Do not draw photons, 1: Draw photons, 2: Draw photons and light paths.
    mutable bool drawDebugAABBs = false;        ///< If true, draw the Axis-Aligned Bounding Boxes (AABBs) of the meshes for debugging.
    bool useKDTree = true;                      ///< Use KD tree for mesh intersection to improve performance, otherwise use a vector (slower).
    int maxKdTreeDepth = 20;                    ///< Maximum depth of the KD tree for mesh intersection.
    TypeOfFloor floorType = CHECKERBOARD;       ///< Type of floor pattern, either CHECKERBOARD or PLAIN.

private:
    Settings() = default;
    ~Settings() = default;
};

#endif //SETTINGS_H