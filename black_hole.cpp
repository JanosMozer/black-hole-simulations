#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <chrono>
#include <fstream>
#include <sstream>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace glm;
using namespace std;
using Clock = std::chrono::high_resolution_clock;

// VARS
double lastPrintTime = 0.0;
int    framesCount   = 0;
double c = 299792458.0;
double G = 6.67430e-11;
struct Ray;
bool Gravity = false;

struct Camera {
    // Center the camera orbit on the black hole at (0, 0, 0)
    vec3 target = vec3(0.0f, 0.0f, 0.0f); // Always look at the black hole center
    float radius = 6.34194e10f;
    float minRadius = 1e10f, maxRadius = 1e12f;

    float azimuth = 0.0f;
    float elevation = M_PI / 2.0f;

    float orbitSpeed = 0.01f;
    float panSpeed = 0.01f;
    double zoomSpeed = 25e9f;

    bool dragging = false;
    bool panning = false;
    bool moving = false; // For compute shader optimization
    double lastX = 0.0, lastY = 0.0;

    // Calculate camera position in world space
    vec3 position() const {
        float clampedElevation = glm::clamp(elevation, 0.01f, float(M_PI) - 0.01f);
        // Orbit around (0,0,0) always
        return vec3(
            radius * sin(clampedElevation) * cos(azimuth),
            radius * cos(clampedElevation),
            radius * sin(clampedElevation) * sin(azimuth)
        );
    }
    void update() {
        // Always keep target at black hole center
        target = vec3(0.0f, 0.0f, 0.0f);
        if(dragging | panning) {
            moving = true;
        } else {
            moving = false;
        }
    }

    void processMouseMove(double x, double y) {
        float dx = float(x - lastX);
        float dy = float(y - lastY);

        if (dragging && panning) {
            // Pan: Shift + Left or Middle Mouse
            // Disable panning to keep camera centered on black hole
        }
        else if (dragging && !panning) {
            // Orbit: Left mouse only
            azimuth   += dx * orbitSpeed;
            elevation -= dy * orbitSpeed;
            elevation = glm::clamp(elevation, 0.01f, float(M_PI) - 0.01f);
        }

        lastX = x;
        lastY = y;
        update();
    }
    void processMouseButton(int button, int action, int mods, GLFWwindow* win) {
        if (button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_MIDDLE) {
            if (action == GLFW_PRESS) {
                dragging = true;
                // Disable panning so camera always orbits center
                panning = false;
                glfwGetCursorPos(win, &lastX, &lastY);
            } else if (action == GLFW_RELEASE) {
                dragging = false;
                panning = false;
            }
        }
        if (button == GLFW_MOUSE_BUTTON_RIGHT) {
            if (action == GLFW_PRESS) {
                Gravity = true;
            } else if (action == GLFW_RELEASE) {
                Gravity = false;
            }
        }
    }
    void processScroll(double xoffset, double yoffset) {
        radius -= yoffset * zoomSpeed;
        radius = glm::clamp(radius, minRadius, maxRadius);
        update();
    }
    void processKey(int key, int scancode, int action, int mods) {
        if (action == GLFW_PRESS && key == GLFW_KEY_G) {
            Gravity = !Gravity;
            cout << "[INFO] Gravity turned " << (Gravity ? "ON" : "OFF") << endl;
        }
    }
};
Camera camera;

struct BlackHole {
    vec3 position;
    double mass;
    double radius;
    double r_s;

    BlackHole(vec3 pos, float m) : position(pos), mass(m) {r_s = 2.0 * G * mass / (c*c);}
    bool Intercept(float px, float py, float pz) const {
        double dx = double(px) - double(position.x);
        double dy = double(py) - double(position.y);
        double dz = double(pz) - double(position.z);
        double dist2 = dx * dx + dy * dy + dz * dz;
        return dist2 < r_s * r_s;
    }
};
BlackHole SagA(vec3(0.0f, 0.0f, 0.0f), 8.54e36); // Sagittarius A black hole
struct ObjectData {
    vec4 posRadius; // xyz = position, w = radius
    vec4 color;     // rgb = color, a = unused
    float  mass;
    vec3 velocity = vec3(0.0f, 0.0f, 0.0f); // Initial velocity
};
vector<ObjectData> objects = {
    { vec4(4e11f, 0.0f, 0.0f, 4e10f)   , vec4(1,1,0,1), 1.98892e30 },
    { vec4(0.0f, 0.0f, 4e11f, 4e10f)   , vec4(1,0,0,1), 1.98892e30 },
    { vec4(0.0f, 0.0f, 0.0f, SagA.r_s) , vec4(0,0,0,1), static_cast<float>(SagA.mass)  },
    //{ vec4(6e10f, 0.0f, 0.0f, 5e10f), vec4(0,1,0,1) }
};

void raytraceCPU(vector<unsigned char>& pixels, int W, int H, const Camera& cam, const vector<ObjectData>& objs);

struct Ray;
void rk4Step(Ray& ray, double dλ, double rs);

// -- Physics functions ported from CPU-geodesic.cpp -- //
struct Ray {
    // cartesian coords
    double x, y, z;
    // spherical coords
    double r, theta, phi;
    double dr, dtheta, dphi;
    double E, L; // conserved quantities

    Ray(vec3 pos, vec3 dir, double rs) : x(pos.x), y(pos.y), z(pos.z) {
        // Convert to spherical coordinates
        r = sqrt(x*x + y*y + z*z);
        theta = acos(z / r);
        phi = atan2(y, x);

        // Seed velocities
        double dx = dir.x, dy = dir.y, dz = dir.z;
        dr     = sin(theta)*cos(phi)*dx + sin(theta)*sin(phi)*dy + cos(theta)*dz;
        dtheta = (cos(theta)*cos(phi)*dx + cos(theta)*sin(phi)*dy - sin(theta)*dz) / r;
        dphi   = (-sin(phi)*dx + cos(phi)*dy) / (r * sin(theta));

        // Store conserved quantities
        L = r * r * sin(theta) * dphi;
        double f = 1.0 - rs / r;
        if (f > 0) {
            double dt_dλ = sqrt((dr*dr)/f + r*r*dtheta*dtheta + r*r*sin(theta)*sin(theta)*dphi*dphi);
            E = f * dt_dλ;
        } else {
            E = 0;
        }
    }

    void step(double dλ, double rs) {
        if (r <= rs) return;
        rk4Step(*this, dλ, rs);
        // convert back to cartesian
        x = r * sin(theta) * cos(phi);
        y = r * sin(theta) * sin(phi);
        z = r * cos(theta);
    }
};

void geodesicRHS(const Ray& ray, double rhs[6], double rs) {
    double r = ray.r, theta = ray.theta, dr = ray.dr, dtheta = ray.dtheta, dphi = ray.dphi, E = ray.E;
    double f = 1.0 - rs / r;
    double dt_dlambda = (f > 0) ? E / f : 0;

    // First derivatives
    rhs[0] = dr;
    rhs[1] = dtheta;
    rhs[2] = dphi;
    // Second derivatives
    rhs[3] = - (rs / (2 * r * r)) * f * dt_dlambda * dt_dlambda + (rs / (2 * r * r * f)) * dr * dr + r * (dtheta * dtheta + sin(theta) * sin(theta) * dphi * dphi);
    rhs[4] = - (2.0 / r) * dr * dtheta + sin(theta) * cos(theta) * dphi * dphi;
    rhs[5] = - (2.0 / r) * dr * dphi - 2.0 * cos(theta) / sin(theta) * dtheta * dphi;
}

void addState(const double a[6], const double b[6], double factor, double out[6]) {
    for (int i = 0; i < 6; i++) out[i] = a[i] + b[i] * factor;
}

void rk4Step(Ray& ray, double dλ, double rs) {
    double y0[6] = { ray.r, ray.theta, ray.phi, ray.dr, ray.dtheta, ray.dphi };
    double k1[6], k2[6], k3[6], k4[6], temp[6];

    geodesicRHS(ray, k1, rs);
    addState(y0, k1, dλ/2.0, temp);
    Ray r2 = ray; r2.r=temp[0]; r2.theta=temp[1]; r2.phi=temp[2]; r2.dr=temp[3]; r2.dtheta=temp[4]; r2.dphi=temp[5];
    geodesicRHS(r2, k2, rs);

    addState(y0, k2, dλ/2.0, temp);
    Ray r3 = ray; r3.r=temp[0]; r3.theta=temp[1]; r3.phi=temp[2]; r3.dr=temp[3]; r3.dtheta=temp[4]; r3.dphi=temp[5];
    geodesicRHS(r3, k3, rs);

    addState(y0, k3, dλ, temp);
    Ray r4 = ray; r4.r=temp[0]; r4.theta=temp[1]; r4.phi=temp[2]; r4.dr=temp[3]; r4.dtheta=temp[4]; r4.dphi=temp[5];
    geodesicRHS(r4, k4, rs);

    ray.r      += (dλ/6.0)*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
    ray.theta  += (dλ/6.0)*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    ray.phi    += (dλ/6.0)*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);
    ray.dr     += (dλ/6.0)*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3]);
    ray.dtheta += (dλ/6.0)*(k1[4] + 2*k2[4] + 2*k3[4] + k4[4]);
    ray.dphi   += (dλ/6.0)*(k1[5] + 2*k2[5] + 2*k3[5] + k4[5]);
}

void raytraceCPU(vector<unsigned char>& pixels, int W, int H, const Camera& cam, const vector<ObjectData>& objs) {
    pixels.resize(W * H * 3);

    // build camera basis
    vec3 forward = normalize(cam.target - cam.position());
    vec3 right   = normalize(cross(forward, vec3(0,1,0)));
    vec3 up      = cross(right, forward);
    float aspect = float(W) / float(H);
    float tanHalfFov = tan(radians(60.0f) * 0.5f);

    ObjectData blackHole;
    for(const auto& obj : objs) {
        if(obj.mass > 1.98892e30 * 100) {
            blackHole = obj;
            break;
        }
    }
    double rs = 2.0 * G * blackHole.mass / (c * c);

    #pragma omp parallel for schedule(dynamic, 4)
    for(int y = 0; y < H; ++y) {
        for(int x = 0; x < W; ++x) {
            float u = (2.0f * (x + 0.5f) / float(W)  - 1.0f) * aspect * tanHalfFov;
            float v = (1.0f - 2.0f * (y + 0.5f) / float(H)) * tanHalfFov;
            vec3 dir = normalize(u*right + v*up + forward);

            Ray ray(cam.position(), dir, rs);

            const int MAX_STEPS = 1000;
            const double D_LAMBDA = 5e8;
            const double ESCAPE_R = 6e12;

            vec4 color(0.0f, 0.0f, 0.05, 1.0); // Background color
            
            for(int i = 0; i < MAX_STEPS; ++i) {
                ray.step(D_LAMBDA, rs);
                vec3 rayPos(ray.x, ray.y, ray.z);

                for(const auto& obj : objs) {
                    if (distance(rayPos, vec3(obj.posRadius)) <= obj.posRadius.w) {
                        color = obj.color;
                        goto pixel_done;
                    }
                }

                if (ray.r > ESCAPE_R) break;
            }
            pixel_done:;

            int idx = (y * W + x) * 3;
            pixels[idx+0] = (unsigned char)(color.r * 255);
            pixels[idx+1] = (unsigned char)(color.g * 255);
            pixels[idx+2] = (unsigned char)(color.b * 255);
        }
    }
}


struct Engine {
    GLuint gridShaderProgram;
    // -- Quad & Texture render -- //
    GLFWwindow* window;
    GLuint quadVAO;
    GLuint texture;
    GLuint shaderProgram;
    GLuint computeProgram = 0;
    GLuint objectShaderProgram = 0;

    // -- UBOs -- //
    GLuint cameraUBO = 0;
    GLuint diskUBO = 0;
    GLuint objectsUBO = 0;
    // -- grid mess vars -- //
    GLuint gridVAO = 0;
    GLuint gridVBO = 0;
    GLuint gridEBO = 0;
    int gridIndexCount = 0;
    int sphereIndexCount = 0;
    GLuint sphereVAO = 0;

    int WIDTH = 800;  // Window width
    int HEIGHT = 600; // Window height
    int COMPUTE_WIDTH  = 200;   // Compute resolution width
    int COMPUTE_HEIGHT = 150;  // Compute resolution height
    float width = 100000000000.0f; // Width of the viewport in meters
    float height = 75000000000.0f; // Height of the viewport in meters
    
    Engine() {
        cout << "[INFO] Initializing Black Hole Engine..." << endl;
        
        if (!glfwInit()) {
            cerr << "[ERROR] GLFW init failed" << endl;
            exit(EXIT_FAILURE);
        }
        cout << "[INFO] GLFW initialized successfully" << endl;
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
        glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GLFW_TRUE);
        cout << "[INFO] Creating GLFW window with OpenGL 4.1..." << endl;
        window = glfwCreateWindow(WIDTH, HEIGHT, "Black Hole", nullptr, nullptr);
        if (!window) {
            cerr << "[ERROR] Failed to create GLFW window" << endl;
            glfwTerminate();
            exit(EXIT_FAILURE);
        }
        cout << "[INFO] GLFW window created successfully" << endl;
        glfwMakeContextCurrent(window);
        cout << "[INFO] OpenGL context made current" << endl;
        
        glewExperimental = GL_TRUE;
        GLenum glewErr = glewInit();
        if (glewErr != GLEW_OK) {
            cerr << "[ERROR] Failed to initialize GLEW: "
                << (const char*)glewGetErrorString(glewErr)
                << endl;
            glfwTerminate();
            exit(EXIT_FAILURE);
        }
        cout << "[INFO] GLEW initialized successfully" << endl;
        
        cout << "[INFO] OpenGL " << glGetString(GL_VERSION) << endl;
        cout << "[INFO] OpenGL Renderer: " << glGetString(GL_RENDERER) << endl;
        cout << "[INFO] OpenGL Vendor: " << glGetString(GL_VENDOR) << endl;
        
        cout << "[INFO] Creating shader programs..." << endl;
        this->shaderProgram = CreateShaderProgram();
        if (this->shaderProgram == 0) {
            cerr << "[ERROR] Failed to create main shader program" << endl;
            exit(EXIT_FAILURE);
        }
        cout << "[INFO] Main shader program created successfully" << endl;
        
        gridShaderProgram = CreateShaderProgram("grid.vert", "grid.frag");
        if (gridShaderProgram == 0) {
            cerr << "[ERROR] Failed to create grid shader program" << endl;
            exit(EXIT_FAILURE);
        }
        cout << "[INFO] Grid shader program created successfully" << endl;

        cout << "[INFO] Creating object shader program..." << endl;
        const char* objectVertSrc = R"(
        #version 410 core
        layout (location = 0) in vec3 aPos;
        uniform mat4 model;
        uniform mat4 viewProj;
        void main() {
            gl_Position = viewProj * model * vec4(aPos, 1.0);
        }
        )";
        const char* objectFragSrc = R"(
        #version 410 core
        out vec4 FragColor;
        uniform vec4 objColor;
        void main() {
            FragColor = objColor;
        }
        )";
        objectShaderProgram = CreateShaderProgramFromSource(objectVertSrc, objectFragSrc);
        if (objectShaderProgram == 0) {
            cerr << "[ERROR] Failed to create object shader program" << endl;
            exit(EXIT_FAILURE);
        }
        cout << "[INFO] Object shader program created successfully" << endl;

        cout << "[INFO] Creating compute program..." << endl;
        computeProgram = CreateComputeProgram("geodesic.comp");
        if (computeProgram == 0) {
            cout << "[WARNING] Compute shaders not supported, using fallback rendering" << endl;
        } else {
            cout << "[INFO] Compute program created successfully" << endl;
        }
        cout << "[INFO] Creating uniform buffers..." << endl;
        glGenBuffers(1, &cameraUBO);
        glBindBuffer(GL_UNIFORM_BUFFER, cameraUBO);
        glBufferData(GL_UNIFORM_BUFFER, 128, nullptr, GL_DYNAMIC_DRAW);
        glBindBufferBase(GL_UNIFORM_BUFFER, 1, cameraUBO);
        cout << "[INFO] Camera UBO created" << endl;

        glGenBuffers(1, &diskUBO);
        glBindBuffer(GL_UNIFORM_BUFFER, diskUBO);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(float) * 4, nullptr, GL_DYNAMIC_DRAW);
        glBindBufferBase(GL_UNIFORM_BUFFER, 2, diskUBO);
        cout << "[INFO] Disk UBO created" << endl;

        glGenBuffers(1, &objectsUBO);
        glBindBuffer(GL_UNIFORM_BUFFER, objectsUBO);
        GLsizeiptr objUBOSize = sizeof(int) + 3 * sizeof(float)
            + 16 * (sizeof(vec4) + sizeof(vec4))
            + 16 * sizeof(float);
        glBufferData(GL_UNIFORM_BUFFER, objUBOSize, nullptr, GL_DYNAMIC_DRAW);
        glBindBufferBase(GL_UNIFORM_BUFFER, 3, objectsUBO);
        cout << "[INFO] Objects UBO created" << endl;

        cout << "[INFO] Creating quad VAO and texture..." << endl;
        auto result = QuadVAO();
        this->quadVAO = result[0];
        this->texture = result[1];
        cout << "[INFO] Quad VAO and texture created successfully" << endl;

        cout << "[INFO] Creating sphere mesh..." << endl;
        createSphere();
        cout << "[INFO] Sphere mesh created" << endl;
        
        cout << "[INFO] Black Hole Engine initialization complete!" << endl;
    }
    void generateGrid(const vector<ObjectData>& objects) {
        const int gridSize = 25;
        const float spacing = 1e10f;  // tweak this

        vector<vec3> vertices;
        vector<GLuint> indices;

        for (int z = 0; z <= gridSize; ++z) {
            for (int x = 0; x <= gridSize; ++x) {
                float worldX = (x - gridSize / 2) * spacing;
                float worldZ = (z - gridSize / 2) * spacing;

                float y = 0.0f;

                // Warp grid using Schwarzschild geometry
                for (const auto& obj : objects) {
                    vec3 objPos = vec3(obj.posRadius);
                    double mass = obj.mass;
                    double radius = obj.posRadius.w;

                    double r_s = 2.0 * G * mass / (c * c);
                    double dx = worldX - objPos.x;
                    double dz = worldZ - objPos.z;
                    double dist = sqrt(dx * dx + dz * dz);

                    // prevent sqrt of negative or divide-by-zero (inside or at the black hole center)
                    if (dist > r_s) {
                        double deltaY = 2.0 * sqrt(r_s * (dist - r_s));
                        y += static_cast<float>(deltaY) - 3e10f;
                    } else {
                        // For points inside or at r_s: make it dip down sharply
                        y += 2.0f * static_cast<float>(sqrt(r_s * r_s)) - 3e10f;  // or add a deep pit
                    }
                }

                vertices.emplace_back(worldX, y, worldZ);
            }
        }

        // Add indices for GL_LINE rendering
        for (int z = 0; z < gridSize; ++z) {
            for (int x = 0; x < gridSize; ++x) {
                int i = z * (gridSize + 1) + x;
                indices.push_back(i);
                indices.push_back(i + 1);

                indices.push_back(i);
                indices.push_back(i + gridSize + 1);
            }
        }

        // Upload to GPU
        if (gridVAO == 0) glGenVertexArrays(1, &gridVAO);
        if (gridVBO == 0) glGenBuffers(1, &gridVBO);
        if (gridEBO == 0) glGenBuffers(1, &gridEBO);

        glBindVertexArray(gridVAO);

        glBindBuffer(GL_ARRAY_BUFFER, gridVBO);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vec3), vertices.data(), GL_DYNAMIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gridEBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), indices.data(), GL_STATIC_DRAW);

        glEnableVertexAttribArray(0); // location = 0
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (void*)0);

        gridIndexCount = indices.size();

        glBindVertexArray(0);
    }
    void drawGrid(const mat4& viewProj) {
        glUseProgram(gridShaderProgram);
        glUniformMatrix4fv(glGetUniformLocation(gridShaderProgram, "viewProj"),
                        1, GL_FALSE, glm::value_ptr(viewProj));
        glBindVertexArray(gridVAO);

        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glDrawElements(GL_LINES, gridIndexCount, GL_UNSIGNED_INT, 0);

        glBindVertexArray(0);
        glEnable(GL_DEPTH_TEST);
    }
    void drawFullScreenQuad() {
        glUseProgram(shaderProgram); // fragment + vertex shader
        glBindVertexArray(quadVAO);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture);
        glUniform1i(glGetUniformLocation(shaderProgram, "screenTexture"), 0);

        glDisable(GL_DEPTH_TEST);  // draw as background
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 6);  // 2 triangles
        glEnable(GL_DEPTH_TEST);
    }
    GLuint CreateShaderProgram(){
        const char* vertexShaderSource = R"(
        #version 330 core
        layout (location = 0) in vec2 aPos;  // Changed to vec2
        layout (location = 1) in vec2 aTexCoord;
        out vec2 TexCoord;
        void main() {
            gl_Position = vec4(aPos, 0.0, 1.0);  // Explicit z=0
            TexCoord = aTexCoord;
        })";

        const char* fragmentShaderSource = R"(
        #version 330 core
        in vec2 TexCoord;
        out vec4 FragColor;
        uniform sampler2D screenTexture;
        void main() {
            FragColor = texture(screenTexture, TexCoord);
        })";

        // vertex shader
        GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertexShader, 1, &vertexShaderSource, nullptr);
        glCompileShader(vertexShader);

        // fragment shader
        GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
        glCompileShader(fragmentShader);

        GLuint shaderProgram = glCreateProgram();
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
        glLinkProgram(shaderProgram);

        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);

        return shaderProgram;
    };
    GLuint CreateShaderProgram(const char* vertPath, const char* fragPath) {
        auto loadShader = [](const char* path, GLenum type) -> GLuint {
            std::ifstream in(path);
            if (!in.is_open()) {
                std::cerr << "Failed to open shader: " << path << "\n";
                exit(EXIT_FAILURE);
            }
            std::stringstream ss;
            ss << in.rdbuf();
            std::string srcStr = ss.str();
            const char* src = srcStr.c_str();

            GLuint shader = glCreateShader(type);
            glShaderSource(shader, 1, &src, nullptr);
            glCompileShader(shader);

            GLint success;
            glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
            if (!success) {
                GLint logLen;
                glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLen);
                std::vector<char> log(logLen);
                glGetShaderInfoLog(shader, logLen, nullptr, log.data());
                std::cerr << "Shader compile error (" << path << "):\n" << log.data() << "\n";
                exit(EXIT_FAILURE);
            }
            return shader;
        };

        GLuint vertShader = loadShader(vertPath, GL_VERTEX_SHADER);
        GLuint fragShader = loadShader(fragPath, GL_FRAGMENT_SHADER);

        GLuint program = glCreateProgram();
        glAttachShader(program, vertShader);
        glAttachShader(program, fragShader);
        glLinkProgram(program);

        GLint linkSuccess;
        glGetProgramiv(program, GL_LINK_STATUS, &linkSuccess);
        if (!linkSuccess) {
            GLint logLen;
            glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLen);
            std::vector<char> log(logLen);
            glGetProgramInfoLog(program, logLen, nullptr, log.data());
            std::cerr << "Shader link error:\n" << log.data() << "\n";
            exit(EXIT_FAILURE);
        }

        glDeleteShader(vertShader);
        glDeleteShader(fragShader);

        return program;
    }
    GLuint CreateComputeProgram(const char* path) {
        if (!GLEW_ARB_compute_shader) {
            cout << "[WARNING] Compute shaders not supported (requires OpenGL 4.3+). Using fallback rendering." << endl;
            return 0;
        }

        // 1) read GLSL source
        std::ifstream in(path);
        if(!in.is_open()) {
            std::cerr << "Failed to open compute shader: " << path << "\n";
            exit(EXIT_FAILURE);
        }
        std::stringstream ss;
        ss << in.rdbuf();
        std::string srcStr = ss.str();
        const char* src = srcStr.c_str();

        // 2) compile
        GLuint cs = glCreateShader(GL_COMPUTE_SHADER);
        glShaderSource(cs, 1, &src, nullptr);
        glCompileShader(cs);
        GLint ok; 
        glGetShaderiv(cs, GL_COMPILE_STATUS, &ok);
        if(!ok) {
            GLint logLen;
            glGetShaderiv(cs, GL_INFO_LOG_LENGTH, &logLen);
            std::vector<char> log(logLen);
            glGetShaderInfoLog(cs, logLen, nullptr, log.data());
            cerr << "[ERROR] Compute shader compile error:\n" << log.data() << endl;
            glDeleteShader(cs);
            return 0;
        }

        // 3) link
        GLuint prog = glCreateProgram();
        glAttachShader(prog, cs);
        glLinkProgram(prog);
        glGetProgramiv(prog, GL_LINK_STATUS, &ok);
        if(!ok) {
            GLint logLen;
            glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &logLen);
            std::vector<char> log(logLen);
            glGetProgramInfoLog(prog, logLen, nullptr, log.data());
            cerr << "[ERROR] Compute shader link error:\n" << log.data() << endl;
            glDeleteShader(cs);
            glDeleteProgram(prog);
            return 0;
        }

        glDeleteShader(cs);
        return prog;
    }

    GLuint CreateShaderProgramFromSource(const char* vertSrc, const char* fragSrc) {
        GLuint vs = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vs, 1, &vertSrc, nullptr);
        glCompileShader(vs);
        GLint ok;
        glGetShaderiv(vs, GL_COMPILE_STATUS, &ok);
        if(!ok) {
            GLint logLen;
            glGetShaderiv(vs, GL_INFO_LOG_LENGTH, &logLen);
            vector<char> log(logLen);
            glGetShaderInfoLog(vs, logLen, nullptr, log.data());
            cerr << "[ERROR] Vertex shader compile error:\n" << log.data() << endl;
            return 0;
        }

        GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fs, 1, &fragSrc, nullptr);
        glCompileShader(fs);
        glGetShaderiv(fs, GL_COMPILE_STATUS, &ok);
        if(!ok) {
            GLint logLen;
            glGetShaderiv(fs, GL_INFO_LOG_LENGTH, &logLen);
            vector<char> log(logLen);
            glGetShaderInfoLog(fs, logLen, nullptr, log.data());
            cerr << "[ERROR] Fragment shader compile error:\n" << log.data() << endl;
            return 0;
        }

        GLuint prog = glCreateProgram();
        glAttachShader(prog, vs);
        glAttachShader(prog, fs);
        glLinkProgram(prog);
        glGetProgramiv(prog, GL_LINK_STATUS, &ok);
        if(!ok) {
            GLint logLen;
            glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &logLen);
            vector<char> log(logLen);
            glGetProgramInfoLog(prog, logLen, nullptr, log.data());
            cerr << "[ERROR] Shader link error:\n" << log.data() << endl;
            return 0;
        }

        glDeleteShader(vs);
        glDeleteShader(fs);
        return prog;
    }

    void createSphere() {
        vector<vec3> vertices;
        vector<GLuint> indices;
        int stackCount = 20;
        int sectorCount = 20;
        float radius = 1.0f;

        for (int i = 0; i <= stackCount; ++i) {
            float stackAngle = M_PI / 2 - i * M_PI / stackCount;
            float xy = radius * cosf(stackAngle);
            float z = radius * sinf(stackAngle);
            for (int j = 0; j <= sectorCount; ++j) {
                float sectorAngle = j * 2 * M_PI / sectorCount;
                float x = xy * cosf(sectorAngle);
                float y = xy * sinf(sectorAngle);
                vertices.push_back(vec3(x, y, z));
            }
        }

        for(int i = 0; i < stackCount; ++i) {
            int k1 = i * (sectorCount + 1);
            int k2 = k1 + sectorCount + 1;
            for(int j = 0; j < sectorCount; ++j, ++k1, ++k2) {
                if (i != 0) {
                    indices.push_back(k1);
                    indices.push_back(k2);
                    indices.push_back(k1 + 1);
                }
                if (i != (stackCount-1)) {
                    indices.push_back(k1 + 1);
                    indices.push_back(k2);
                    indices.push_back(k2 + 1);
                }
            }
        }
        sphereIndexCount = indices.size();

        glGenVertexArrays(1, &sphereVAO);
        glBindVertexArray(sphereVAO);
        GLuint vbo, ebo;
        glGenBuffers(1, &vbo);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vec3), vertices.data(), GL_STATIC_DRAW);
        
        glGenBuffers(1, &ebo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), indices.data(), GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (void*)0);
        glEnableVertexAttribArray(0);
    }

    void drawObjects(mat4 viewProj, const vector<ObjectData>& objs) {
        glUseProgram(objectShaderProgram);
        glUniformMatrix4fv(glGetUniformLocation(objectShaderProgram, "viewProj"), 1, GL_FALSE, value_ptr(viewProj));
        glBindVertexArray(sphereVAO);
        for (const auto& obj : objs) {
            mat4 model = mat4(1.0f);
            model = translate(model, vec3(obj.posRadius));
            model = scale(model, vec3(obj.posRadius.w));
            glUniformMatrix4fv(glGetUniformLocation(objectShaderProgram, "model"), 1, GL_FALSE, value_ptr(model));
            glUniform4fv(glGetUniformLocation(objectShaderProgram, "objColor"), 1, value_ptr(obj.color));
            glDrawElements(GL_TRIANGLES, sphereIndexCount, GL_UNSIGNED_INT, 0);
        }
    }

    void dispatchCompute(const Camera& cam) {
        // Check if compute shaders are supported
        if (!computeProgram) {
            // Fallback: create a simple colored texture
            glBindTexture(GL_TEXTURE_2D, texture);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 200, 150, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
            
            // Create a more interesting fallback pattern
            std::vector<unsigned char> fallbackData(200 * 150 * 4, 0);
            for (int i = 0; i < 200 * 150; ++i) {
                int x = i % 200;
                int y = i / 200;
                
                // Black hole in center
                float distFromCenter = sqrt((x - 100) * (x - 100) + (y - 75) * (y - 75));
                if (distFromCenter < 20) {
                    fallbackData[i * 4 + 0] = 255;  // Red (black hole)
                    fallbackData[i * 4 + 1] = 0;
                    fallbackData[i * 4 + 2] = 0;
                } else if (distFromCenter < 40) {
                    fallbackData[i * 4 + 0] = 255;  // Orange (accretion disk)
                    fallbackData[i * 4 + 1] = 165;
                    fallbackData[i * 4 + 2] = 0;
                } else {
                    // Stars and background
                    if (x % 50 == 0 || y % 50 == 0) {
                        fallbackData[i * 4 + 0] = 255;  // White stars
                        fallbackData[i * 4 + 1] = 255;
                        fallbackData[i * 4 + 2] = 255;
                    } else {
                        fallbackData[i * 4 + 0] = 20;   // Dark blue background
                        fallbackData[i * 4 + 1] = 20;
                        fallbackData[i * 4 + 2] = 100;
                    }
                }
                fallbackData[i * 4 + 3] = 255;  // A
            }
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 200, 150, GL_RGBA, GL_UNSIGNED_BYTE, fallbackData.data());
            return;
        }
        
        // determine target compute‐res
        int cw = cam.moving ? COMPUTE_WIDTH  : 200;
        int ch = cam.moving ? COMPUTE_HEIGHT : 150;

        // 1) reallocate the texture if needed
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexImage2D(GL_TEXTURE_2D,
                    0,                // mip
                    GL_RGBA8,         // internal format
                    cw,               // width
                    ch,               // height
                    0, GL_RGBA, 
                    GL_UNSIGNED_BYTE, 
                    nullptr);

        // 2) bind compute program & UBOs
        glUseProgram(computeProgram);
        uploadCameraUBO(cam);
        uploadDiskUBO();
        uploadObjectsUBO(objects);

        // 3) bind it as image unit 0
        glBindImageTexture(0, texture, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA8);

        // 4) dispatch grid
        GLuint groupsX = (GLuint)std::ceil(cw / 16.0f);
        GLuint groupsY = (GLuint)std::ceil(ch / 16.0f);
        glDispatchCompute(groupsX, groupsY, 1);

        // 5) sync
        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    }
    void uploadCameraUBO(const Camera& cam) {
        struct UBOData {
            vec3 pos; float _pad0;
            vec3 right; float _pad1;
            vec3 up; float _pad2;
            vec3 forward; float _pad3;
            float tanHalfFov;
            float aspect;
            bool moving;
            int _pad4;
        } data;
        vec3 fwd = normalize(cam.target - cam.position());
        vec3 up = vec3(0, 1, 0); // y axis is up, so disk is in x-z plane
        vec3 right = normalize(cross(fwd, up));
        up = cross(right, fwd);

        data.pos = cam.position();
        data.right = right;
        data.up = up;
        data.forward = fwd;
        data.tanHalfFov = tan(radians(60.0f * 0.5f));
        data.aspect = float(WIDTH) / float(HEIGHT);
        data.moving = cam.dragging || cam.panning;

        glBindBuffer(GL_UNIFORM_BUFFER, cameraUBO);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(UBOData), &data);
    }
    void uploadObjectsUBO(const vector<ObjectData>& objs) {
        struct UBOData {
            int   numObjects;
            float _pad0, _pad1, _pad2;        // <-- pad out to 16 bytes
            vec4  posRadius[16];
            vec4  color[16];
            float  mass[16]; 
        } data;

        size_t count = std::min(objs.size(), size_t(16));
        data.numObjects = static_cast<int>(count);

        for (size_t i = 0; i < count; ++i) {
            data.posRadius[i] = objs[i].posRadius;
            data.color[i] = objs[i].color;
            data.mass[i] = objs[i].mass;
        }

        // Upload
        glBindBuffer(GL_UNIFORM_BUFFER, objectsUBO);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(data), &data);
    }
    void uploadDiskUBO() {
        // disk
        float r1 = SagA.r_s * 2.2f;    // inner radius just outside the event horizon
        float r2 = SagA.r_s * 5.2f;   // outer radius of the disk
        float num = 2.0;               // number of rays
        float thickness = 1e9f;          // padding for std140 alignment
        float diskData[4] = { r1, r2, num, thickness };

        glBindBuffer(GL_UNIFORM_BUFFER, diskUBO);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(diskData), diskData);
    }
    
    vector<GLuint> QuadVAO(){
        float quadVertices[] = {
            // positions   // texCoords
            -1.0f,  1.0f,  0.0f, 1.0f,  // top left
            -1.0f, -1.0f,  0.0f, 0.0f,  // bottom left
            1.0f, -1.0f,  1.0f, 0.0f,  // bottom right

            -1.0f,  1.0f,  0.0f, 1.0f,  // top left
            1.0f, -1.0f,  1.0f, 0.0f,  // bottom right
            1.0f,  1.0f,  1.0f, 1.0f   // top right
        };
        
        GLuint VAO, VBO;
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);

        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);

        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
        glEnableVertexAttribArray(1);

        GLuint texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexImage2D(GL_TEXTURE_2D,
                    0,             // mip
                    GL_RGBA8,      // internal format
                    COMPUTE_WIDTH,
                    COMPUTE_HEIGHT,
                    0,
                    GL_RGBA,
                    GL_UNSIGNED_BYTE,
                    nullptr);
        vector<GLuint> VAOtexture = {VAO, texture};
        return VAOtexture;
    }
    void renderScene() {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(shaderProgram);
        glBindVertexArray(quadVAO);
        // make sure your fragment shader samples from texture unit 0:
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glfwSwapBuffers(window);
        glfwPollEvents();
    };
};
Engine engine;
void setupCameraCallbacks(GLFWwindow* window) {
    glfwSetWindowUserPointer(window, &camera);

    glfwSetMouseButtonCallback(window, [](GLFWwindow* win, int button, int action, int mods) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processMouseButton(button, action, mods, win);
    });

    glfwSetCursorPosCallback(window, [](GLFWwindow* win, double x, double y) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processMouseMove(x, y);
    });

    glfwSetScrollCallback(window, [](GLFWwindow* win, double xoffset, double yoffset) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processScroll(xoffset, yoffset);
    });

    glfwSetKeyCallback(window, [](GLFWwindow* win, int key, int scancode, int action, int mods) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processKey(key, scancode, action, mods);
    });
}


// -- MAIN -- //
int main() {
    cout << "[INFO] Starting Black Hole 3D Simulation..." << endl;
    cout << "[INFO] Setting up camera callbacks..." << endl;
    setupCameraCallbacks(engine.window);
    cout << "[INFO] Camera callbacks setup complete" << endl;
    
    vector<unsigned char> pixels(engine.WIDTH * engine.HEIGHT * 3);

    auto t0 = Clock::now();
    lastPrintTime = chrono::duration<double>(t0.time_since_epoch()).count();

    double lastTime = glfwGetTime();
    int   renderW  = 800, renderH = 600, numSteps = 80000;
    
    cout << "[INFO] Entering main render loop..." << endl;
    int frameCount = 0;
    while (!glfwWindowShouldClose(engine.window)) {
        frameCount++;
        if (frameCount % 60 == 0) {
            cout << "[INFO] Frame " << frameCount << " - Camera pos: (" 
                 << camera.position().x << ", " << camera.position().y << ", " << camera.position().z << ")" << endl;
        }
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);  // optional, but good practice
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        double now   = glfwGetTime();
        double dt    = now - lastTime;   // seconds since last frame
        lastTime     = now;

        mat4 view = lookAt(camera.position(), camera.target, vec3(0,1,0));
        mat4 proj = perspective(radians(60.0f), float(engine.WIDTH)/float(engine.HEIGHT), 1e9f, 1e14f);
        mat4 viewProj = proj * view;

        if(engine.computeProgram) {
            // ---------- RUN RAYTRACER ------------- //
            glViewport(0, 0, engine.WIDTH, engine.HEIGHT);
            engine.dispatchCompute(camera);
            engine.drawFullScreenQuad();
        } else {
            // Fallback CPU ray tracing
            int texWidth = 400;
            int texHeight = 300;
            static vector<unsigned char> pixels(texWidth * texHeight * 3);
            
            raytraceCPU(pixels, texWidth, texHeight, camera, objects);

            glBindTexture(GL_TEXTURE_2D, engine.texture);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texWidth, texHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
            
            engine.drawFullScreenQuad();
        }
        
        // Gravity
        for (auto& obj : objects) {
            for (auto& obj2 : objects) {
                if (&obj == &obj2) continue; // skip self-interaction
                 float dx  = obj2.posRadius.x - obj.posRadius.x;
                 float dy = obj2.posRadius.y - obj.posRadius.y;
                 float dz = obj2.posRadius.z - obj.posRadius.z;
                 float distance = sqrt(dx * dx + dy * dy + dz * dz);
                 if (distance > 0) {
                        vector<double> direction = {dx / distance, dy / distance, dz / distance};
                        //distance *= 1000;
                        double Gforce = (G * obj.mass * obj2.mass) / (distance * distance);

                        double acc1 = Gforce / obj.mass;
                        std::vector<double> acc = {direction[0] * acc1, direction[1] * acc1, direction[2] * acc1};
                        if (Gravity) {
                            obj.velocity.x += acc[0];
                            obj.velocity.y += acc[1];
                            obj.velocity.z += acc[2];

                            obj.posRadius.x += obj.velocity.x;
                            obj.posRadius.y += obj.velocity.y;
                            obj.posRadius.z += obj.velocity.z;
                            cout << "velocity: " <<obj.velocity.x<<", " <<obj.velocity.y<<", " <<obj.velocity.z<<endl;
                        }
                    }
            }
        }



        // ---------- GRID ------------- //
        // 2) rebuild grid mesh on CPU
        //engine.generateGrid(objects);
        // 5) overlay the bent grid
        //mat4 view = lookAt(camera.position(), camera.target, vec3(0,1,0));
        //mat4 proj = perspective(radians(60.0f), float(engine.COMPUTE_WIDTH)/engine.COMPUTE_HEIGHT, 1e9f, 1e14f);
        //mat4 viewProj = proj * view;
        //engine.drawGrid(viewProj);

        // ---------- RUN RAYTRACER ------------- //
        //glViewport(0, 0, engine.WIDTH, engine.HEIGHT);
        //engine.dispatchCompute(camera);
        //engine.drawFullScreenQuad();

        // 6) present to screen
        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    cout << "[INFO] Render loop ended. Cleaning up..." << endl;
    glfwDestroyWindow(engine.window);
    glfwTerminate();
    cout << "[INFO] Black Hole 3D Simulation terminated successfully" << endl;
    return 0;
}