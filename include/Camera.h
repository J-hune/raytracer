#ifndef CAMERA_H
#define CAMERA_H

#include "Vec3.h"

class Camera final {
public:
    Camera();
    virtual ~Camera() = default;

    /**
     * Gets the field of view angle of the camera.
     * @return The field of view angle.
     */
    [[nodiscard]] constexpr float getFovAngle() const noexcept { return fovAngle; }

    /**
     * Sets the field of view angle of the camera.
     * @param newFovAngle The new field of view angle.
     */
    constexpr void setFovAngle(const float newFovAngle) noexcept { fovAngle = newFovAngle; }

    /**
     * Gets the aspect ratio of the camera.
     * @return The aspect ratio.
     */
    [[nodiscard]] constexpr float getAspectRatio() const noexcept { return aspectRatio; }

    /**
     * Gets the near plane distance of the camera.
     * @return The near plane distance.
     */
    [[nodiscard]] constexpr float getNearPlane() const noexcept { return nearPlane; }

    /**
     * Sets the near plane distance of the camera.
     * @param newNearPlane The new near plane distance.
     */
    constexpr void setNearPlane(const float newNearPlane) noexcept { nearPlane = newNearPlane; }

    /**
     * Gets the far plane distance of the camera.
     * @return The far plane distance.
     */
    [[nodiscard]] constexpr float getFarPlane() const noexcept { return farPlane; }

    /**
     * Sets the far plane distance of the camera.
     * @param newFarPlane The new far plane distance.
     */
    constexpr void setFarPlane(const float newFarPlane) noexcept { farPlane = newFarPlane; }

    /**
     * Gets the screen width.
     * @return The screen width.
     */
    [[nodiscard]] constexpr unsigned int getScreenWidth() const noexcept { return width; }

    /**
     * Gets the screen height.
     * @return The screen height.
     */
    [[nodiscard]] constexpr unsigned int getScreenHeight() const noexcept { return height; }

    /**
     * Resizes the camera's viewport.
     * @param width The new width of the viewport.
     * @param height The new height of the viewport.
     */
    void resize(int width, int height);

    /**
     * Initializes the camera's position.
     */
    void initPos();

    /**
     * Moves the camera by the specified deltas.
     * @param dx The delta movement along the x-axis.
     * @param dy The delta movement along the y-axis.
     * @param dz The delta movement along the z-axis.
     */
    void move(float dx, float dy, float dz) noexcept;

    /**
     * Begins the rotation of the camera.
     * @param u The initial horizontal coordinate.
     * @param v The initial vertical coordinate.
     */
    void beginRotate(int u, int v) noexcept;

    /**
     * Rotates the camera based on the given coordinates.
     * @param u The current horizontal coordinate.
     * @param v The current vertical coordinate.
     */
    void rotate(int u, int v);

    /**
     * Ends the rotation of the camera.
     */
    void endRotate() noexcept { moving = false; }

    /**
     * Zooms the camera by the specified amount.
     * @param z The zoom level to add.
     */
    void zoom(const float z) noexcept { zoomLevel += z; }

    /**
     * Applies the camera transformations.
     */
    void apply();

    /**
     * Gets the current position of the camera.
     * @param x Reference to store the x-coordinate of the camera.
     * @param y Reference to store the y-coordinate of the camera.
     * @param z Reference to store the z-coordinate of the camera.
     */
    void getPos(float &x, float &y, float &z) noexcept;

    /**
     * Gets the current position of the camera.
     * @param pos Reference to a Vec3 object to store the position of the camera.
     */
    void getPos(Vec3 &pos) noexcept { getPos(pos[0], pos[1], pos[2]); }

private:
    float fovAngle{45.0f};              ///< Field of view angle
    float aspectRatio{1.0f};            ///< Aspect ratio
    float nearPlane{4.1f};              ///< Near plane distance
    float farPlane{10000.0f};           ///< Far plane distance

    int spinning{0}, moving{0};         ///< Flags for spinning and moving the camera
    int startU{0}, startV{0};           ///< Initial coordinates for rotation
    int height{0}, width{0};            ///< Screen height and width
    float currentQuat[4];               ///< Current quaternion
    float lastQuat[4];                  ///< Last quaternion
    float x{0.0f}, y{0.0f}, z{0.0f};    ///< Camera position
    float zoomLevel{3.0f};              ///< Zoom level
    bool initialized{false};            ///< Flag to check if the camera has been initialized

    /**
     * Resets the quaternion.
     */
    void resetQuat();
};

#endif // CAMERA_H
