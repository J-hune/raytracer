#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <mutex>

#include <algorithm>
#include <condition_variable>
#include <iomanip>
#include <filesystem>
#include <queue>
#include <random>

#include "src/Vec3.h"
#include "src/Camera.h"
#include "src/Scene.h"
#include "src/Settings.h"
#include <GL/glut.h>

#include "src/matrixUtilities.h"

using namespace std;

// -------------------------------------------
// OpenGL/GLUT application code.
// -------------------------------------------
struct FastRenderSettings {
    int width;
    int height;
    int samples;
    int shadowRays;
};

struct RenderTask {
    int startY;
    int endY;
};

static GLint window;
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
std::atomic totalProcessedRows(0);
std::atomic lastPrintedProgress(-1);
std::mutex progressMutex;
std::queue<RenderTask> taskQueue;
std::mutex queueMutex;
std::condition_variable queueCondition;


// Print the usage information for the application.
void printUsage() {
    cerr << endl
        << "=========================================" << endl
        << "         Path Tracer - Usage Guide       " << endl
        << "=========================================" << endl
        << endl
        << "Keyboard Commands:" << endl
        << "------------------" << endl
        << "  ?   : Display help" << endl
        << "  w   : Toggle Wireframe Mode" << endl
        << "  f   : Toggle Full Screen Mode" << endl
        << "  p   : Quick Render" << endl
        << "  r   : Render Render" << endl
        << "  +   : Switch Scene" << endl
        << "  d   : Debug draw mode" << endl
        << endl
        << "Mouse Controls:" << endl
        << "---------------" << endl
        << "  Drag + Left Button   : Rotate Model" << endl
        << "  Drag + Right Button  : Move Model" << endl
        << "  Drag + Middle Button : Zoom" << endl
        << endl
        << "  q or esc : Quit Application" << endl
        << endl;
}

// Display usage information and exit.
void usage() {
    printUsage();
    exit(EXIT_FAILURE);
}

// ------------------------------------
// Initialize the lighting settings.
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

// Initialize the OpenGL context and settings.
void init() {
    const Settings &settings = Settings::getInstance();
    camera.resize(settings.width, settings.height);
    initLight();
    // glCullFace(GL_BACK); // Uncomment if back face culling is needed.
    glDisable(GL_CULL_FACE);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.2f, 0.2f, 0.3f, 1.0f);
}

void draw() {
    glEnable(GL_LIGHTING);
    scenes[selected_scene].draw();

    // Draw rays for debugging
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);
    glLineWidth(6);
    glColor3f(1, 0, 0);
    glBegin(GL_LINES);

    // Draw rays
    for (unsigned int r = 0; r < rays.size(); ++r) {
        glVertex3f(rays[r].first[0], rays[r].first[1], rays[r].first[2]);
        glVertex3f(rays[r].second[0], rays[r].second[1], rays[r].second[2]);
    }
    glEnd();
}

// Display callback function.
void display() {
    glLoadIdentity();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply();
    draw();
    glFlush();
    glutSwapBuffers();
}

// Idle callback function to update the frame rate.
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

// Print a progress bar to the console.
void printProgressBar(const int progress) {
    // Add the new progress bar with the format [■■■■■■□□□□□□□□□□] 26%
    std::cout << "[";
    for (int i = 0; i < 50; ++i) {
        std::cout << (i < progress / 2 ? "■" : "□");
    }
    std::cout << "] " << progress << "%" << std::flush;
    if (progress == 100) std::cout << std::endl;
}

void export_partial_image(const std::vector<Vec3>& image, const int w, const int h, const std::string& filename) {
    std::ofstream f(filename, std::ios::binary);
    if (!f) {
        std::cout << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    f << "P3\n" << w << " " << h << "\n255\n";
    for (const auto& pixel : image) {
        f << static_cast<int>(255.f * std::min(1.f, pixel[0])) << " "
          << static_cast<int>(255.f * std::min(1.f, pixel[1])) << " "
          << static_cast<int>(255.f * std::min(1.f, pixel[2])) << " ";
    }

    f.close();
}

void ray_trace_section(const int w, const int h, const unsigned int nsamples,
    std::vector<Vec3> &image, std::mt19937 &rng, std::uniform_real_distribution<float> &dist, const int totalRows)
{
    Vec3 pos, dir;

    // Print the launch of ray tracing and thread ID
    std::ostringstream oss;
    oss << std::this_thread::get_id();
    printf("Ray tracing thread %s started\n", oss.str().c_str());

    while (true) {
        RenderTask task{};

        // Get the next task from the queue
        {
            std::unique_lock lock(queueMutex);
            queueCondition.wait(lock, [&] { return !taskQueue.empty(); });
            task = taskQueue.front();
            taskQueue.pop();
        }

        // If the task is invalid, exit the thread
        if (task.startY == -1) break;

        // Ray tracing for the specified section
        for (int y = task.startY; y < task.endY; y++) {
            for (int x = 0; x < w; x++) {
                for (unsigned int s = 0; s < nsamples; ++s) {
                    const float u = (static_cast<float>(x) + dist(rng)) / static_cast<float>(w);
                    const float v = (static_cast<float>(y) + dist(rng)) / static_cast<float>(h);

                    screenSpaceToWorldSpaceRay(u, v, pos, dir);
                    const Vec3 color = scenes[selected_scene].rayTrace(Ray(pos, dir), rng);
                    image[x + y * w] += color;
                }
                image[x + y * w] /= static_cast<float>(nsamples);
            }

            // Update the total processed rows
            const int processedRows = ++totalProcessedRows;
            const int progress = static_cast<int>(std::round(100.0 * processedRows / totalRows));
            if (progress % 2 == 0 && progress > lastPrintedProgress) {
                std::lock_guard lock(progressMutex);

                // Remove the previous progress bar
                for (int i = 0; i < 100; i++) {
                    std::cout << "\b";
                }

                std::cout << "Thread " << oss.str() << " processed row " << task.startY << " - " << task.endY << " (" << progress << "%)" << std::endl;
                if (progress > lastPrintedProgress) {
                    lastPrintedProgress = progress;
                    printProgressBar(progress);
                }

                if (progress % 10 == 0) {
                    std::ostringstream partialFilename;
                    partialFilename << "renders/render_partial.ppm";
                    export_partial_image(image, w, h, partialFilename.str());
                }
            }
        }
    }
}

// Main ray tracing function from the camera perspective.
void ray_trace_from_camera(const Settings &settings) {
    int w = settings.width; int h = settings.height;
    camera.apply();
    std::vector<Vec3> image(w * h, Vec3(0, 0, 0));

    auto start = std::chrono::high_resolution_clock::now();
    initializeMatrices();

    // Random number generator for anti-aliasing
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    // Number of threads to use
    unsigned int numThreads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
    int totalRows = h;
    int bandHeight = 10;
    totalProcessedRows = 0; // Reset the counter
    lastPrintedProgress = -1; // Reset the progress bar
    std::cout << "Ray tracing image of size " << w << "x" << h << " with " << numThreads << " threads and " << bandHeight << " rows per band" << std::endl;

    // Fill the task queue with bands of rows
    for (int y = 0; y < h; y += bandHeight) {
        taskQueue.push({y, std::min(y + bandHeight, h)});
    }

    // Start the threads
    for (unsigned int i = 0; i < numThreads; ++i) {
        threads.emplace_back(ray_trace_section, w, h, settings.samples, std::ref(image), std::ref(rng), std::ref(dist), totalRows);
    }

    // Wait for all threads to finish
    for (unsigned int i = 0; i < numThreads; ++i) {
        taskQueue.push({-1, -1}); // Signal the thread to exit
    }
    queueCondition.notify_all(); // Notify all threads to exit

    // Wait for all threads to finish
    for (auto &thread: threads) {
        thread.join();
    }

    // Calculate elapsed time in seconds
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
    std::cout << "Ray tracing completed in " << static_cast<float>(duration.count()) / 1000.f << " seconds" << std::endl;

    // Export the image after ray tracing
    std::ostringstream oss;
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    oss << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
    std::string time_str = oss.str();

    std::filesystem::create_directories("renders");
    std::string filename = "renders/render_" + time_str + ".ppm";

    std::ofstream f(filename, std::ios::binary);
    if (!f) {
        std::cout << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    f << "P3\n" << w << " " << h << "\n255\n";
    for (const auto &pixel: image) {
        f << static_cast<int>(255.f * std::min(1.f, pixel[0])) << " "
          << static_cast<int>(255.f * std::min(1.f, pixel[1])) << " "
          << static_cast<int>(255.f * std::min(1.f, pixel[2])) << " ";
    }

    f.close();
    std::cout << "Image saved to " << filename << std::endl;
}

void fast_ray_trace_from_camera(const Settings &settings) {
    // Copy the original settings
    FastRenderSettings originalSettings{};
    originalSettings.width = settings.width;
    originalSettings.height = settings.height;
    originalSettings.samples = settings.samples;
    originalSettings.shadowRays = settings.shadowRays;

    // Edit the settings for a faster render and ray trace
    settings.width = 240;
    settings.height = 240;
    settings.samples = 8;
    settings.shadowRays = 6;
    camera.resize(settings.width, settings.height);
    ray_trace_from_camera(settings);

    // Restore the original settings
    settings.width = originalSettings.width;
    settings.height = originalSettings.height;
    settings.samples = originalSettings.samples;
    settings.shadowRays = originalSettings.shadowRays;
    camera.resize(settings.width, settings.height);
}

// Function to handle keyboard inputs.
void keyboard(const unsigned char key, int _x, int _y) {
    const Settings &settings = Settings::getInstance();
    switch (key) {
        case 'f':
            if (fullScreen == true) {
                glutReshapeWindow(settings.width, settings.height);
                fullScreen = false;
            } else {
                glutFullScreen();
                fullScreen = true;
            }
            break;
        case 'd':
            settings.drawDebugPhotons = !settings.drawDebugPhotons;
            break;
        case 'q':
        case 27:
            exit(0);
        case 'w':
            GLint polygonMode[2];
            glGetIntegerv(GL_POLYGON_MODE, polygonMode);
            if (polygonMode[0] != GL_FILL)
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            else
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            break;
        case 'p':
            camera.apply();
            rays.clear();
            fast_ray_trace_from_camera(settings);
            break;
        case 'r':
            camera.apply();
            rays.clear();
            ray_trace_from_camera(settings);
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

// Function to handle mouse motion events.
void motion(const int x, const int y) {
    const Settings &settings = Settings::getInstance();

    if (mouseRotatePressed == true) {
        camera.rotate(x, y);
    } else if (mouseMovePressed == true) {
        camera.move(
            static_cast<float>(x - lastX) / static_cast<float>(settings.width),
            static_cast<float>(lastY - y) / static_cast<float>(settings.height),
            0.0
        );
        lastX = x;
        lastY = y;
    } else if (mouseZoomPressed == true) {
        camera.zoom(static_cast<float>(y - lastZoom) / static_cast<float>(settings.height));
        lastZoom = y;
    }
}


void reshape(const int w, const int h) {
    camera.resize(w, h);
}


int main(int argc, char **argv) {
    Settings &settings = Settings::getInstance();
    settings.width = 480;
    settings.height = 480;
    settings.samples = 80;
    settings.shadowRays = 16;
    settings.photons = 500000;
    settings.directIllumination = true;
    settings.caustics = true;
    settings.reflections = true;
    settings.refractions = true;
    settings.drawDebugPhotons = true;
    settings.floorType = CHECKERBOARD; //PLAIN, CHECKERBOARD (checkerboard is a lot slower)

    if (argc > 2) {
        printUsage();
        exit(EXIT_FAILURE);
    }
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(settings.width, settings.height);
    window = glutCreateWindow("gMini");

    init();
    glutIdleFunc(idle);
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutReshapeFunc(reshape);
    glutMotionFunc(motion);
    glutMouseFunc(mouse);
    keyboard('?', 0, 0);


    camera.move(0., 0., -3.1);
    selected_scene = 3;
    scenes.resize(4);
    scenes[0].setup_single_sphere();
    scenes[1].setup_multiple_spheres();
    scenes[2].setup_single_square();
    scenes[3].setup_cornell_box();
    glutMainLoop();

    return EXIT_SUCCESS;
}
