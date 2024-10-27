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

    int width = 480;
    int height = 480;
    int samples = 100;
    int photon = 20000;
    bool caustics = true;
    bool drawDebugPhotons = false;
    TypeOfFloor floorType = CHECKERBOARD;

private:
    Settings() = default;
    ~Settings() = default;
};

#endif //SETTINGS_H
