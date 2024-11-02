#include "../include/Camera.h"
#include "../include/Trackball.h"
#include <GL/gl.h>
#include <GL/glu.h>

Camera::Camera() {
    resetQuat();
}

void Camera::resize(const int width, const int height) {
    this->width = width;
    this->height = height;
    aspectRatio = static_cast<float>(width) / static_cast<float>(height);
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fovAngle, aspectRatio, nearPlane, farPlane);
    glMatrixMode(GL_MODELVIEW);
}

void Camera::initPos() {
    if (!initialized) {
        spinning = 0;
        moving = 0;
        resetQuat();
        initialized = true;
    }
}

void Camera::resetQuat() {
    trackball(currentQuat, 0.0, 0.0, 0.0, 0.0);
}

void Camera::move(const float dx, const float dy, const float dz) noexcept {
    x += dx;
    y += dy;
    z += dz;
}

void Camera::beginRotate(const int u, const int v) noexcept {
    startU = u;
    startV = v;
    moving = true;
    spinning = false;
}

void Camera::rotate(const int u, const int v) {
    if (moving) {
        trackball(
            lastQuat,
            (2.0f * static_cast<float>(startU - width)) / static_cast<float>(width),
            (static_cast<float>(height) - 2.0f * static_cast<float>(startV)) / static_cast<float>(height),
            (2.0f * static_cast<float>(u - width)) / static_cast<float>(width),
            (static_cast<float>(height) - 2.0f * static_cast<float>(v)) / static_cast<float>(height)
            );
        startU = u;
        startV = v;
        spinning = true;
        add_quats(lastQuat, currentQuat, currentQuat);
    }
}

void Camera::apply() {
    glLoadIdentity();
    glTranslatef(x, y, z);
    GLfloat m[4][4];
    build_rotmatrix(m, currentQuat);
    glTranslatef(0.0f, 0.0f, -zoomLevel);
    glMultMatrixf(&m[0][0]);
}

void Camera::getPos(float& x, float& y, float& z) noexcept {
    GLfloat m[4][4];
    build_rotmatrix(m, currentQuat);
    x = m[0][0] * -x + m[0][1] * -y + m[0][2] * -z;
    y = m[1][0] * -x + m[1][1] * -y + m[1][2] * -z;
    z = m[2][0] * -x + m[2][1] * -y + m[2][2] * -z;
}

