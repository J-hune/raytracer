// -------------------------------------------
// gMini : a minimal OpenGL/GLUT application
// for 3D graphics.
// Copyright (C) 2006-2008 Tamy Boubekeur
// All rights reserved.
// -------------------------------------------

// -------------------------------------------
// Disclaimer: this code is dirty in the
// meaning that there is no attention paid to
// proper class attribute access, memory
// management or optimisation of any kind. It
// is designed for quick-and-dirty testing
// purpose.
// -------------------------------------------


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <iomanip>
#include <filesystem>

#include "src/Vec3.h"
#include "src/Camera.h"
#include "src/Scene.h"
#include <GL/glut.h>

#include "src/matrixUtilities.h"
#include "src/imageLoader.h"
#include "src/Material.h"

using namespace std;

// -------------------------------------------
// OpenGL/GLUT application code.
// -------------------------------------------

static GLint window;
static int SCREENWIDTH = 480;
static int SCREENHEIGHT = 480;
static Camera camera;
static bool mouseRotatePressed = false;
static bool mouseMovePressed = false;
static bool mouseZoomPressed = false;
static int lastX = 0, lastY = 0, lastZoom = 0;
static unsigned int FPS = 0;
static bool fullScreen = false;

std::vector<Scene> scenes;
unsigned int selected_scene;

std::vector<std::pair<Vec3, Vec3>> rays;

void printUsage() {
    cerr << endl
            << "gMini: a minimal OpenGL/GLUT application for 3D graphics." << endl
            << "Author : Tamy Boubekeur (https://www.labri.fr/~boubek)" << endl << endl
            << "Usage : ./gmini [<file.off>]" << endl
            << "Keyboard commands" << endl
            << "------------------" << endl
            << " ?: Print help" << endl
            << " w: Toggle Wireframe Mode" << endl
            << " g: Toggle Gouraud Shading Mode" << endl
            << " f: Toggle full screen mode" << endl
            << " <drag>+<left button>: rotate model" << endl
            << " <drag>+<right button>: move model" << endl
            << " <drag>+<middle button>: zoom" << endl
            << " q, <esc>: Quit" << endl << endl;
}

void usage() {
    printUsage();
    exit(EXIT_FAILURE);
}


// ------------------------------------
void initLight() {
    constexpr GLfloat light_position[4] = {0.0, 1.5, 0.0, 1.0};
    constexpr GLfloat color[4] = {1.0, 1.0, 1.0, 1.0};
    constexpr GLfloat ambient[4] = {1.0, 1.0, 1.0, 1.0};

    glLightfv(GL_LIGHT1, GL_POSITION, light_position);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, color);
    glLightfv(GL_LIGHT1, GL_SPECULAR, color);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHTING);
}

void init() {
    camera.resize(SCREENWIDTH, SCREENHEIGHT);
    initLight();
    //glCullFace (GL_BACK);
    glDisable(GL_CULL_FACE);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.2f, 0.2f, 0.3f, 1.0f);
}


// ------------------------------------
// Replace the code of these
// functions for cleaning memory,
// closing sockets, etc.
// ------------------------------------

void clear() {
}

// ------------------------------------
// Replace the code of these
// functions for alternative rendering.
// ------------------------------------


void draw() {
    glEnable(GL_LIGHTING);
    scenes[selected_scene].draw();

    // draw rays : (for debug)
    //  std::cout << rays.size() << std::endl;
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);
    glLineWidth(6);
    glColor3f(1, 0, 0);
    glBegin(GL_LINES);

    // draw rays
    for (unsigned int r = 0; r < rays.size(); ++r) {
        glVertex3f(rays[r].first[0], rays[r].first[1], rays[r].first[2]);
        glVertex3f(rays[r].second[0], rays[r].second[1], rays[r].second[2]);
    }
    glEnd();
}

void display() {
    glLoadIdentity();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply();
    draw();
    glFlush();
    glutSwapBuffers();
}

void idle() {
    static int lastTime = glutGet(GLUT_ELAPSED_TIME);
    static unsigned int counter = 0;
    counter++;

    // Update FPS every second
    const int currentTime = glutGet(GLUT_ELAPSED_TIME);
    if (currentTime - lastTime >= 1000) {
        FPS = counter;
        counter = 0;
        static char winTitle[64];
        sprintf(winTitle, "Raytracer - FPS: %d", FPS);
        glutSetWindowTitle(winTitle);
        lastTime = currentTime;
    }
    glutPostRedisplay();
}


void ray_trace_from_camera() {
    int w = glutGet(GLUT_WINDOW_WIDTH), h = glutGet(GLUT_WINDOW_HEIGHT);
    std::cout << "Ray tracing a " << w << " x " << h << " image..." << std::endl;
    camera.apply();
    Vec3 pos, dir;
    //    unsigned int samples = 100;
    unsigned int nsamples = 50;
    std::vector<Vec3> image(w * h, Vec3(0, 0, 0));

    // On récupère le temps actuel
    auto start = std::chrono::high_resolution_clock::now();

    // For each pixel in the image, we cast a ray and accumulate the color
    //#pragma omp parallel for
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            for (unsigned int s = 0; s < nsamples; ++s) {
                float u = (static_cast<float>(x) + static_cast<float>(rand()) / static_cast<float>(RAND_MAX)) / static_cast<float>(w);
                float v = (static_cast<float>(y) + static_cast<float>(rand()) / static_cast<float>(RAND_MAX)) / static_cast<float>(h);
                // this is a random uv that belongs to the pixel xy.
                screen_space_to_world_space_ray(u, v, pos, dir);
                Vec3 color = scenes[selected_scene].rayTrace(Ray(pos, dir));
                image[x + y * w] += color;
            }
            image[x + y * w] /= static_cast<float>(nsamples);
        }
    }

    // On calcule le temps écoulé en secondes
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start);
    std::cout << "Ray tracing réalisé en " << duration.count() << " secondes avec " << nsamples << " samples par pixel." << std::endl;

    // Get current time
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    // Format time to string
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
    std::string time_str = oss.str();

    // Create the directory if it does not exist
    std::filesystem::create_directories("rendus");

    // Create filename with date and time
    std::string filename = "rendus/rendu_" + time_str + ".ppm";
    std::ofstream f(filename.c_str(), std::ios::binary);
    if (f.fail()) {
        std::cout << "Could not open file: " << filename << std::endl;
        return;
    }
    f << "P3\n" << w << " " << h << "\n255\n";
    for (int i = 0; i < w * h; i++) {
        f << static_cast<int>(255.f * std::min(1.f, image[i][0])) << " "
          << static_cast<int>(255.f * std::min(1.f, image[i][1])) << " "
          << static_cast<int>(255.f * std::min(1.f, image[i][2])) << " ";
    }
    f << std::endl;
    std::cout << "Image exportée dans " << filename << std::endl;
    f.close();
}


void key(const unsigned char keyPressed, int _x, int _y) {
    Vec3 pos, dir;
    switch (keyPressed) {
        case 'f':
            if (fullScreen == true) {
                glutReshapeWindow(SCREENWIDTH, SCREENHEIGHT);
                fullScreen = false;
            } else {
                glutFullScreen();
                fullScreen = true;
            }
            break;
        case 'q':
        case 27:
            clear();
            exit(0);
        case 'w':
            GLint polygonMode[2];
            glGetIntegerv(GL_POLYGON_MODE, polygonMode);
            if (polygonMode[0] != GL_FILL)
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            else
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            break;

        case 'r':
            camera.apply();
            rays.clear();
            ray_trace_from_camera();
            break;
        case '+':
            selected_scene = (selected_scene + 1) % scenes.size();
            break;
        default:
            printUsage();
            break;
    }
    idle();
}

void mouse(const int button, const int state, const int x, const int y) {
    if (state == GLUT_UP) {
        mouseMovePressed = false;
        mouseRotatePressed = false;
        mouseZoomPressed = false;
    } else {
        if (button == GLUT_LEFT_BUTTON) {
            camera.beginRotate(x, y);
            mouseMovePressed = false;
            mouseRotatePressed = true;
            mouseZoomPressed = false;
        } else if (button == GLUT_RIGHT_BUTTON) {
            lastX = x;
            lastY = y;
            mouseMovePressed = true;
            mouseRotatePressed = false;
            mouseZoomPressed = false;
        } else if (button == GLUT_MIDDLE_BUTTON) {
            if (mouseZoomPressed == false) {
                lastZoom = y;
                mouseMovePressed = false;
                mouseRotatePressed = false;
                mouseZoomPressed = true;
            }
        }
    }
    idle();
}

void motion(int x, int y) {
    if (mouseRotatePressed == true) {
        camera.rotate(x, y);
    } else if (mouseMovePressed == true) {
        camera.move(
            static_cast<float>(x - lastX) / static_cast<float>(SCREENWIDTH),
            static_cast<float>(lastY - y) / static_cast<float>(SCREENHEIGHT),
            0.0
        );
        lastX = x;
        lastY = y;
    } else if (mouseZoomPressed == true) {
        camera.zoom(static_cast<float>(y - lastZoom) / static_cast<float>(SCREENHEIGHT));
        lastZoom = y;
    }
}


void reshape(const int w, const int h) {
    camera.resize(w, h);
}


int main(int argc, char **argv) {
    if (argc > 2) {
        printUsage();
        exit(EXIT_FAILURE);
    }
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(SCREENWIDTH, SCREENHEIGHT);
    window = glutCreateWindow("gMini");

    init();
    glutIdleFunc(idle);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutReshapeFunc(reshape);
    glutMotionFunc(motion);
    glutMouseFunc(mouse);
    key('?', 0, 0);


    camera.move(0., 0., -3.1);
    selected_scene = 0;
    scenes.resize(4);
    scenes[0].setup_single_sphere();
    scenes[1].setup_multiple_spheres();
    scenes[2].setup_single_square();
    scenes[3].setup_cornell_box();

    glutMainLoop();
    return EXIT_SUCCESS;
}
