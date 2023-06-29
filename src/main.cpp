#include <algorithm>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <vector>

#define NANOSVG_IMPLEMENTATION

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include <nanosvg.h>
#define GLAD_GL_IMPLEMENTATION
#include <glad/gl.h>

#include "simple_math.h"

#define SCANLINE_DEBUG

NSVGimage *vgImage = nullptr;

static unsigned char bgColor[4] = {205, 202, 200, 255};
static unsigned char lineColor[4] = {0, 160, 192, 255};

void checkGlError(const char *flag) {
    for (GLint error = glGetError(); error; error = glGetError()) {
        std::cout << "Error " << error << " after " << flag << std::endl;
    }
}

static float distPtSeg(float x, float y, float px, float py, float qx, float qy) {
    float pqx, pqy, dx, dy, d, t;
    pqx = qx - px;
    pqy = qy - py;
    dx = x - px;
    dy = y - py;
    d = pqx * pqx + pqy * pqy;
    t = pqx * dx + pqy * dy;

    if (d > 0) t /= d;

    if (t < 0)
        t = 0;
    else if (t > 1)
        t = 1;

    dx = px + t * pqx - x;
    dy = py + t * pqy - y;

    return dx * dx + dy * dy;
}

static void cubicBez(
    float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4, float tol, int level) {
    float x12, y12, x23, y23, x34, y34, x123, y123, x234, y234, x1234, y1234;
    float d;

    // Reach recursion limit.
    if (level > 12) return;

    x12 = (x1 + x2) * 0.5f;
    y12 = (y1 + y2) * 0.5f;
    x23 = (x2 + x3) * 0.5f;
    y23 = (y2 + y3) * 0.5f;
    x34 = (x3 + x4) * 0.5f;
    y34 = (y3 + y4) * 0.5f;
    x123 = (x12 + x23) * 0.5f;
    y123 = (y12 + y23) * 0.5f;
    x234 = (x23 + x34) * 0.5f;
    y234 = (y23 + y34) * 0.5f;
    x1234 = (x123 + x234) * 0.5f;
    y1234 = (y123 + y234) * 0.5f;

    // Distance.
    d = distPtSeg(x1234, y1234, x1, y1, x4, y4);

    // Recursion.
    if (d > tol * tol) {
        cubicBez(x1, y1, x12, y12, x123, y123, x1234, y1234, tol, level + 1);
        cubicBez(x1234, y1234, x234, y234, x34, y34, x4, y4, tol, level + 1);
    } else {
        glVertex2f(x4, y4);
    }
}

void drawPath(float *pts, int npts, char closed, float tol) {
    glBegin(GL_LINE_STRIP);
    glColor4ubv(lineColor);
    glVertex2f(pts[0], pts[1]);

    for (int i = 0; i < npts - 1; i += 3) {
        float *p = &pts[i * 2];
        cubicBez(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], tol, 0);
    }

    if (closed) {
        glVertex2f(pts[0], pts[1]);
    }

    glEnd();
}

void drawControlPts(float *pts, int npts) {
    // Control lines.
    glColor4ubv(lineColor);
    glBegin(GL_LINES);
    for (int i = 0; i < npts - 1; i += 3) {
        float *p = &pts[i * 2];
        glVertex2f(p[0], p[1]);
        glVertex2f(p[2], p[3]);
        glVertex2f(p[4], p[5]);
        glVertex2f(p[6], p[7]);
    }
    glEnd();

    // Central points
    glPointSize(6.0f);
    glColor4ubv(lineColor);

    glBegin(GL_POINTS);
    glVertex2f(pts[0], pts[1]);
    for (int i = 0; i < npts - 1; i += 3) {
        float *p = &pts[i * 2];
        glVertex2f(p[6], p[7]);
    }
    glEnd();

    // Points
    glPointSize(3.0f);

    glBegin(GL_POINTS);
    glColor4ubv(bgColor);
    glVertex2f(pts[0], pts[1]);
    for (int i = 0; i < npts - 1; i += 3) {
        float *p = &pts[i * 2];
        glColor4ubv(lineColor);
        glVertex2f(p[2], p[3]);
        glVertex2f(p[4], p[5]);

        // Make the central point hollow.
        glColor4ubv(bgColor);
        glVertex2f(p[6], p[7]);
    }
    glEnd();
}

void drawFrame(GLFWwindow *window, double delta) {
    static int firstFrame = 1;

    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

    // Get framebuffer size.
    int width = 0, height = 0;
    glfwGetFramebufferSize(window, &width, &height);

    glViewport(0, 0, width, height);
    checkGlError("glViewport");

    glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_BLEND);

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_TEXTURE_2D);

    // Specify which matrix is the current matrix
    glMatrixMode(GL_PROJECTION);

    // Replace the current matrix with the identity matrix
    glLoadIdentity();

    // Fit view to bounds
    float cx = vgImage->width * 0.5f;
    float cy = vgImage->height * 0.5f;
    float hw = vgImage->width * 0.5f;
    float hh = vgImage->height * 0.5f;

    float view[4], aspect;
    if ((float)width / hw < (float)height / hh) {
        aspect = (float)height / (float)width;
        view[0] = cx - hw * 1.2f;
        view[2] = cx + hw * 1.2f;
        view[1] = cy - hw * 1.2f * aspect;
        view[3] = cy + hw * 1.2f * aspect;
    } else {
        aspect = (float)width / (float)height;
        view[0] = cx - hh * 1.2f * aspect;
        view[2] = cx + hh * 1.2f * aspect;
        view[1] = cy - hh * 1.2f;
        view[3] = cy + hh * 1.2f;
    }

    // Multiply the current matrix with an orthographic matrix.
    glOrtho(view[0], view[2], view[3], view[1], -1, 1);

    // Specify which matrix is the current matrix.
    glMatrixMode(GL_MODELVIEW);

    // Replace the current matrix with the identity matrix.
    glLoadIdentity();

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Draw the image's bounding box.
    glColor4ub(255, 0, 0, 255);
    glBegin(GL_LINE_LOOP);
    glVertex2f(0, 0);
    glVertex2f(vgImage->width, 0);
    glVertex2f(vgImage->width, vgImage->height);
    glVertex2f(0, vgImage->height);
    glEnd();

    // Traverse the scanlines.
    for (int scanline_y = 0; scanline_y < (int)vgImage->height; scanline_y++) {
        // Scanline formula.
        float lx[] = {0.0f, vgImage->width};
        float ly[] = {(float)scanline_y, (float)scanline_y};

#ifdef SCANLINE_DEBUG
        // Draw scanline
        glColor4ub(220, 220, 220, 50);
        glLineWidth(1);
        glBegin(GL_LINES);
        glVertex2f(lx[0], ly[0]);
        glVertex2f(lx[1], ly[1]);
        glEnd();
#endif

        // Traverse shapes.
        for (NSVGshape *shape = vgImage->shapes; shape != nullptr; shape = shape->next) {
            // Intersection points.
            std::vector<std::vector<float>> intersection_points;

            // Get the shape fill color.
            unsigned int cr = shape->fill.color & 0xff;
            unsigned int cg = (shape->fill.color >> 8) & 0xff;
            unsigned int cb = (shape->fill.color >> 16) & 0xff;
            unsigned int ca = (shape->fill.color >> 24) & 0xff;

            // Traverse paths in the shape.
            for (NSVGpath *path = shape->paths; path != nullptr; path = path->next) {
                // Coordinates of the intersection points.
                float I[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

                // Traverse bezier curves in the path.
                for (int b_curve = 0; b_curve < path->npts / 3; b_curve += 1) {
                    float *p = &path->pts[b_curve * 2 * 3];

                    // Calculate intersection points.
                    compute_intersections(p, lx, ly, I);

                    for (int i = 0; i < 3; i++) {
                        // Store valid intersection points and discard invalid ones.
                        if (abs(I[i * 2]) < 10000.0f && abs(I[i * 2 + 1]) < 10000.0f) {
                            // std::cout << i << " - " << I[i * 2] << "," << I[i * 2 + 1] << std::endl;
                            std::vector<float> ip{I[i * 2], I[i * 2 + 1]};
                            intersection_points.push_back(ip);
                        }
                    }
                }
            }

            // Sort intersection points by X
            sort(intersection_points.begin(),
                 intersection_points.end(),
                 [](const std::vector<float> &a, const std::vector<float> &b) { return a[0] < b[0]; });

            // Draw image line by line
            bool to_fill = false;
            bool is_first_point = true;
            for (int i = 0; i < intersection_points.size(); i++) {
                if (!is_first_point) {
                    if (to_fill) {
                        glColor4ub(cr, cg, cb, ca);
                        glBegin(GL_LINES);
                        glVertex2f(intersection_points[i - 1][0], scanline_y);
                        glVertex2f(intersection_points[i][0], scanline_y);
                        glEnd();
                    }
                } else {
                    is_first_point = false;
                }
                to_fill = !to_fill;
            }

#ifdef SCANLINE_DEBUG
            // Draw intersection points.
            for (auto &point : intersection_points) {
                glColor4ub(255, 0, 0, 255);
                glBegin(GL_POINTS);
                glVertex2f(point[0], point[1]);
                // printf(" (%.1f, %.1f)", intersection_points[i][0], intersection_points[i][1]);
                glEnd();
            }
#endif
        }
    }

    if (firstFrame > 0) firstFrame = 0;

    glfwSwapBuffers(window);
}

void resize_cb(GLFWwindow *window, int width, int height) {
    // Update and render
    NSVG_NOTUSED(width);
    NSVG_NOTUSED(height);
    drawFrame(window, 0.0);
}

int main() {
    if (!glfwInit()) {
        return -1;
    }

    const GLFWvidmode *mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    int screenWidth = mode->width;
    int screenHeight = mode->height;

    GLFWwindow *window =
        glfwCreateWindow(screenWidth - 800, screenHeight - 400, "Scanline VG Renderer", nullptr, nullptr);
    if (!window) {
        printf("Could not create window!\n");
        glfwTerminate();
        return -1;
    }

    glfwSetFramebufferSizeCallback(window, resize_cb);
    glfwMakeContextCurrent(window);

    // GLAD: load all OpenGL function pointers.
    if (!gladLoadGL(glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    vgImage = nsvgParseFromFile("../assets/features.svg", "px", 96.0f);
    if (vgImage == nullptr) {
        printf("Could not open the SVG file!\n");
        glfwTerminate();
        return -1;
    }

    auto timeLastFrame = std::chrono::steady_clock::now();
    auto timeLastPrint = timeLastFrame;
    std::chrono::duration<double> elapsedSeconds{};

    unsigned int frame = 0;

    while (!glfwWindowShouldClose(window)) {
        auto now = std::chrono::steady_clock::now();

        // In seconds
        elapsedSeconds = now - timeLastFrame;
        double delta = elapsedSeconds.count();

        elapsedSeconds = now - timeLastPrint;
        if (elapsedSeconds.count() > 1) {
            timeLastPrint = now;
            std::cout << "Frame time: " << round(delta * 1000.0) << " ms." << std::endl;
        }

        timeLastFrame = std::chrono::steady_clock::now();
        frame++;

        drawFrame(window, delta);

        glfwPollEvents();
    }

    nsvgDelete(vgImage);

    glfwTerminate();

    return 0;
}
