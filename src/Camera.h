#ifndef CAMERA_H
#define CAMERA_H

#include "Vec3.h"

class Camera {
public:
  Camera();
  virtual ~Camera() = default;

  [[nodiscard]] constexpr float getFovAngle() const noexcept { return fovAngle; }
  constexpr void setFovAngle(const float newFovAngle) noexcept { fovAngle = newFovAngle; }
  [[nodiscard]] constexpr float getAspectRatio() const noexcept { return aspectRatio; }
  [[nodiscard]] constexpr float getNearPlane() const noexcept { return nearPlane; }
  constexpr void setNearPlane(const float newNearPlane) noexcept { nearPlane = newNearPlane; }
  [[nodiscard]] constexpr float getFarPlane() const noexcept { return farPlane; }
  constexpr void setFarPlane(const float newFarPlane) noexcept { farPlane = newFarPlane; }
  [[nodiscard]] constexpr unsigned int getScreenWidth() const noexcept { return width; }
  [[nodiscard]] constexpr unsigned int getScreenHeight() const noexcept { return height; }

  void resize(int width, int height);
  void initPos();
  void move(float dx, float dy, float dz) noexcept;
  void beginRotate(int u, int v) noexcept;
  void rotate(int u, int v);
  void endRotate() noexcept { moving = false; }
  void zoom(const float z) noexcept { zoomLevel += z; }
  void apply();

  void getPos(float &x, float &y, float &z) noexcept;
  void getPos(Vec3 &pos) noexcept { getPos(pos[0], pos[1], pos[2]); }

private:
  float fovAngle{45.0f};
  float aspectRatio{1.0f};
  float nearPlane{4.1f};
  float farPlane{10000.0f};

  int spinning{0}, moving{0};
  int startU{0}, startV{0};
  int height{0}, width{0};
  float currentQuat[4];
  float lastQuat[4];
  float x{0.0f}, y{0.0f}, z{0.0f};
  float zoomLevel{3.0f};
  bool initialized{false};

  void resetQuat();
};

#endif // CAMERA_H