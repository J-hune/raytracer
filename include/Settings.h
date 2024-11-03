#ifndef SETTINGS_H
#define SETTINGS_H

#pragma once

enum TypeOfFloor {
    CHECKERBOARD,
    PLAIN
};

class Settings {
public:
    static Settings& getInstance() {
        static Settings instance;
        return instance;
    }
    Settings(const Settings&) = delete;
    Settings& operator=(const Settings&) = delete;

    mutable int width = 480;
    mutable int height = 480;
    mutable int samples = 100;
    mutable int shadowRays = 16;
    int photons = 20000;
    bool directIllumination = true;
    bool reflections = true;
    bool refractions = true;
    bool caustics = true;
    mutable bool drawDebugPhotons = false;
    mutable bool drawDebugAABBs = false;
    TypeOfFloor floorType = CHECKERBOARD;

private:
    Settings() = default;
    ~Settings() = default;
};

#endif //SETTINGS_H
