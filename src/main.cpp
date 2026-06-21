#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>

struct Vec3 {
    float x, y, z;
    Vec3(float x=0, float y=0, float z=0) : x(x), y(y), z(z) {}
    Vec3 operator+(const Vec3& o) const { return Vec3(x+o.x, y+o.y, z+o.z); }
    Vec3 operator-(const Vec3& o) const { return Vec3(x-o.x, y-o.y, z-o.z); }
    Vec3 operator*(float s) const { return Vec3(x*s, y*s, z*s); }
    Vec3 operator/(float s) const { return Vec3(x/s, y/s, z/s); }
    float dot(const Vec3& o) const { return x*o.x + y*o.y + z*o.z; }
    float magnitude() const { return sqrt(dot(*this)); }
    Vec3 normalize() const { float m = magnitude(); return m>0 ? *this/m : Vec3(); }
    Vec3 cross(const Vec3& o) const {
        return Vec3(y*o.z - z*o.y, z*o.x - x*o.z, x*o.y - y*o.x);
    }
};

struct Color {
    float r, g, b;
    Color(float r=0, float g=0, float b=0) : r(r), g(g), b(b) {}
};

struct Electron {
    int n, l, ml, ms;
    Vec3 pos;
    Vec3 vel;
    float energy;
};

struct Nucleon {
    Vec3 pos;
    Color color;
    float radius;
    bool isProton;
};

struct Atom {
    int Z;
    int N;
    std::vector<Nucleon> nucleons;
    std::vector<Electron> electrons;
    float nuclearRadius;
};

struct OrbitalPoint {
    Vec3 pos;
    Color color;
    float density;
};

struct Camera {
    Vec3 pos;
    float theta, phi, distance;
    float fov;
};

static float lerp(float a, float b, float t) { return a + (b-a)*t; }

static Color colorLerp(const Color& a, const Color& b, float t) {
    return Color(lerp(a.r,b.r,t), lerp(a.g,b.g,t), lerp(a.b,b.b,t));
}

static Vec3 sphericalToCart(float r, float theta, float phi) {
    return Vec3(r*sin(theta)*cos(phi), r*cos(theta), r*sin(theta)*sin(phi));
}

static std::vector<int> getElectronConfig(int Z) {
    std::vector<int> config;
    const char* order[] = {"1s","2s","2p","3s","3p","4s","3d","4p","5s","4d","5p","6s","4f","5d","6p","7s","5f","6d","7p"};
    int maxCap[] = {2,2,6,2,6,2,10,6,2,10,6,2,14,10,6,2,14,10,6};
    int total = 0;
    for(int i=0; i<19 && total<Z; i++) {
        int cap = maxCap[i];
        if(total + cap > Z) cap = Z - total;
        config.push_back(cap);
        total += cap;
    }
    return config;
}

static int getSubshellL(int index) {
    int lMap[] = {0,0,1,0,1,0,2,1,0,2,1,0,3,2,1,0,3,2,1};
    if(index < 19) return lMap[index];
    return 0;
}

static int getSubshellN(int index) {
    int nMap[] = {1,2,2,3,3,4,3,4,5,4,5,6,4,5,6,7,5,6,7};
    if(index < 19) return nMap[index];
    return 1;
}

static float getElectronEnergy(int n, int l) {
    return -13.6f / (n*n);
}

static Vec3 getOrbitalPosition(int n, int l, int ml, int ms, float time) {
    std::random_device rd;
    std::mt19937 gen(rd() + (int)(time * 1000) + n*100 + l*10 + ml);
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    std::normal_distribution<float> noise(0.0f, 0.15f);
    std::mt19937 noiseGen(rd() + (int)(time * 2000) + n*200 + l*20 + ml*5);
    
    float r, theta, phi;
    
    if(l == 0) {
        float rho = sqrt(-log(dist(gen)));
        r = rho * 0.5f * n * n;
        theta = acos(2*dist(gen)-1);
        phi = 2 * 3.14159f * dist(gen);
    }
    else if(l == 1) {
        float rho = sqrt(-log(dist(gen)));
        r = rho * 0.4f * n * n;
        theta = acos(2*dist(gen)-1);
        phi = 2 * 3.14159f * dist(gen);
        
        if(ml == -1) {
            float tmp = theta;
            theta = acos(2*dist(gen)-1);
            phi = 2 * 3.14159f * dist(gen);
            Vec3 p = sphericalToCart(r, theta, phi);
            p.x = fabs(p.x);
            if(dist(gen) > 0.5f) p.x = -p.x;
            p.x += noise(noiseGen) * 0.3f;
            p.y += noise(noiseGen) * 0.3f;
            p.z += noise(noiseGen) * 0.3f;
            return p;
        }
        else if(ml == 0) {
            Vec3 p = sphericalToCart(r, theta, phi);
            p.y = fabs(p.y);
            if(dist(gen) > 0.5f) p.y = -p.y;
            p.x += noise(noiseGen) * 0.3f;
            p.y += noise(noiseGen) * 0.3f;
            p.z += noise(noiseGen) * 0.3f;
            return p;
        }
        else if(ml == 1) {
            Vec3 p = sphericalToCart(r, theta, phi);
            p.z = fabs(p.z);
            if(dist(gen) > 0.5f) p.z = -p.z;
            p.x += noise(noiseGen) * 0.3f;
            p.y += noise(noiseGen) * 0.3f;
            p.z += noise(noiseGen) * 0.3f;
            return p;
        }
    }
    else if(l == 2) {
        float rho = sqrt(-log(dist(gen)));
        r = rho * 0.35f * n * n;
        theta = acos(2*dist(gen)-1);
        phi = 2 * 3.14159f * dist(gen);
        
        if(ml == -2) {
            Vec3 p = sphericalToCart(r, theta, phi);
            p.x = fabs(p.x); p.y = fabs(p.y);
            if(dist(gen) > 0.5f) p.x = -p.x;
            if(dist(gen) > 0.5f) p.y = -p.y;
            p.x += noise(noiseGen) * 0.3f;
            p.y += noise(noiseGen) * 0.3f;
            p.z += noise(noiseGen) * 0.3f;
            return p;
        }
        else if(ml == -1) {
            Vec3 p = sphericalToCart(r, theta, phi);
            p.x = fabs(p.x); p.z = fabs(p.z);
            if(dist(gen) > 0.5f) p.x = -p.x;
            if(dist(gen) > 0.5f) p.z = -p.z;
            p.x += noise(noiseGen) * 0.3f;
            p.y += noise(noiseGen) * 0.3f;
            p.z += noise(noiseGen) * 0.3f;
            return p;
        }
        else if(ml == 0) {
            Vec3 p = sphericalToCart(r, theta, phi);
            float factor = 3*cos(theta)*cos(theta) - 1;
            p = p * fabs(factor);
            p.x += noise(noiseGen) * 0.3f;
            p.y += noise(noiseGen) * 0.3f;
            p.z += noise(noiseGen) * 0.3f;
            return p;
        }
        else if(ml == 1) {
            Vec3 p = sphericalToCart(r, theta, phi);
            p.y = fabs(p.y); p.z = fabs(p.z);
            if(dist(gen) > 0.5f) p.y = -p.y;
            if(dist(gen) > 0.5f) p.z = -p.z;
            p.x += noise(noiseGen) * 0.3f;
            p.y += noise(noiseGen) * 0.3f;
            p.z += noise(noiseGen) * 0.3f;
            return p;
        }
        else if(ml == 2) {
            Vec3 p = sphericalToCart(r, theta, phi);
            p.x = fabs(p.x); p.y = fabs(p.y);
            if(dist(gen) > 0.5f) p.x = -p.x;
            if(dist(gen) > 0.5f) p.y = -p.y;
            p.x += noise(noiseGen) * 0.3f;
            p.y += noise(noiseGen) * 0.3f;
            p.z += noise(noiseGen) * 0.3f;
            return p;
        }
    }
    else if(l == 3) {
        float rho = sqrt(-log(dist(gen)));
        r = rho * 0.3f * n * n;
        theta = acos(2*dist(gen)-1);
        phi = 2 * 3.14159f * dist(gen);
        Vec3 p = sphericalToCart(r, theta, phi);
        p.x += noise(noiseGen) * 0.3f;
        p.y += noise(noiseGen) * 0.3f;
        p.z += noise(noiseGen) * 0.3f;
        return p;
    }
    
    Vec3 p = sphericalToCart(r, theta, phi);
    p.x += noise(noiseGen) * 0.3f;
    p.y += noise(noiseGen) * 0.3f;
    p.z += noise(noiseGen) * 0.3f;
    return p;
}

static Color getOrbitalColor(int l, int ml) {
    if(l == 0) return Color(0.2f, 0.5f, 1.0f);
    if(l == 1) {
        if(ml == -1) return Color(1.0f, 0.2f, 0.2f);
        if(ml == 0) return Color(0.2f, 1.0f, 0.2f);
        if(ml == 1) return Color(0.2f, 0.2f, 1.0f);
    }
    if(l == 2) {
        if(ml == -2) return Color(1.0f, 0.5f, 0.0f);
        if(ml == -1) return Color(1.0f, 0.0f, 0.5f);
        if(ml == 0) return Color(0.5f, 0.5f, 0.0f);
        if(ml == 1) return Color(0.0f, 1.0f, 0.5f);
        if(ml == 2) return Color(0.5f, 0.0f, 1.0f);
    }
    return Color(0.5f, 0.5f, 0.5f);
}

static Atom createAtom(int Z, int N) {
    Atom atom;
    atom.Z = Z;
    atom.N = N;
    atom.nuclearRadius = 0.3f * pow(Z + N, 1.0f/3.0f);
    
    atom.nucleons.clear();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    
    for(int i=0; i<Z+N; i++) {
        Nucleon nuc;
        float r = atom.nuclearRadius * pow(fabs(dist(gen)), 2.0f) * 0.8f;
        float theta = acos(2*dist(gen)-1);
        float phi = 2 * 3.14159f * dist(gen);
        nuc.pos = sphericalToCart(r, theta, phi);
        nuc.radius = 0.08f * pow(Z+N, 1.0f/3.0f);
        nuc.isProton = (i < Z);
        nuc.color = nuc.isProton ? Color(0.9f, 0.2f, 0.2f) : Color(0.5f, 0.5f, 0.5f);
        atom.nucleons.push_back(nuc);
    }
    
    atom.electrons.clear();
    std::vector<int> config = getElectronConfig(Z);
    int electronIndex = 0;
    
    for(int shellIdx=0; shellIdx<(int)config.size(); shellIdx++) {
        int electronsInShell = config[shellIdx];
        int n = getSubshellN(shellIdx);
        int l = getSubshellL(shellIdx);
        
        int mlValues = 2*l + 1;
        int elecInSubshell = 0;
        
        for(int ml = -l; ml <= l && elecInSubshell < electronsInShell; ml++) {
            for(int ms = -1; ms <= 1 && elecInSubshell < electronsInShell; ms += 2) {
                if(elecInSubshell < electronsInShell) {
                    Electron e;
                    e.n = n;
                    e.l = l;
                    e.ml = ml;
                    e.ms = ms;
                    e.energy = getElectronEnergy(n, l);
                    e.pos = getOrbitalPosition(n, l, ml, ms, electronIndex * 1.5f);
                    e.vel = Vec3(0,0,0);
                    atom.electrons.push_back(e);
                    elecInSubshell++;
                    electronIndex++;
                }
            }
        }
    }
    
    return atom;
}

static std::vector<OrbitalPoint> generateOrbitalPoints(const Atom& atom, int pointsPerOrbital) {
    std::vector<OrbitalPoint> points;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    
    float nucleusExclusionRadius = atom.nuclearRadius * 1.8f;
    
    for(const Electron& e : atom.electrons) {
        Color col = getOrbitalColor(e.l, e.ml);
        float rFactor = 0.3f + e.n * 0.15f;
        
        int numPoints = pointsPerOrbital / atom.electrons.size();
        if(numPoints < 10) numPoints = 10;
        
        for(int i=0; i<numPoints; i++) {
            Vec3 pos = getOrbitalPosition(e.n, e.l, e.ml, e.ms, i * 0.7f + e.n * 1.3f + e.ml * 0.5f);
            pos = pos * rFactor;
            
            float distFromCenter = pos.magnitude();
            if(distFromCenter < nucleusExclusionRadius) {
                pos = pos.normalize() * nucleusExclusionRadius;
            }
            
            float density = 1.0f - pos.magnitude() / (e.n * e.n * 1.5f);
            density = std::max(0.0f, std::min(1.0f, density));
            
            Color pc = col;
            pc.r *= 0.3f + 0.7f * density;
            pc.g *= 0.3f + 0.7f * density;
            pc.b *= 0.3f + 0.7f * density;
            
            OrbitalPoint op;
            op.pos = pos;
            op.color = pc;
            op.density = density;
            points.push_back(op);
        }
    }
    
    return points;
}

static GLuint createShader(const char* src, GLenum type) {
    GLuint sh = glCreateShader(type);
    glShaderSource(sh, 1, &src, nullptr);
    glCompileShader(sh);
    GLint ok;
    glGetShaderiv(sh, GL_COMPILE_STATUS, &ok);
    if(!ok) {
        char log[512];
        glGetShaderInfoLog(sh, sizeof(log), nullptr, log);
        std::cerr << "Shader error: " << log << std::endl;
    }
    return sh;
}

static GLuint createGridVAO(float size, int divisions) {
    std::vector<float> vertices;
    float step = size / divisions;
    float halfSize = size * 0.5f;
    float col = 0.15f;
    float colAxis = 0.4f;
    
    for(int i=0; i<=divisions; i++) {
        float pos = -halfSize + i * step;
        vertices.push_back(pos);
        vertices.push_back(0.0f);
        vertices.push_back(-halfSize);
        vertices.push_back(col);
        vertices.push_back(col);
        vertices.push_back(col * 1.5f);
        
        vertices.push_back(pos);
        vertices.push_back(0.0f);
        vertices.push_back(halfSize);
        vertices.push_back(col);
        vertices.push_back(col);
        vertices.push_back(col * 1.5f);
        
        vertices.push_back(-halfSize);
        vertices.push_back(0.0f);
        vertices.push_back(pos);
        vertices.push_back(col);
        vertices.push_back(col);
        vertices.push_back(col * 1.5f);
        
        vertices.push_back(halfSize);
        vertices.push_back(0.0f);
        vertices.push_back(pos);
        vertices.push_back(col);
        vertices.push_back(col);
        vertices.push_back(col * 1.5f);
    }
    
    float axisLen = halfSize * 1.2f;
    vertices.push_back(-axisLen);
    vertices.push_back(0.0f);
    vertices.push_back(0.0f);
    vertices.push_back(0.8f);
    vertices.push_back(0.1f);
    vertices.push_back(0.1f);
    
    vertices.push_back(axisLen);
    vertices.push_back(0.0f);
    vertices.push_back(0.0f);
    vertices.push_back(0.8f);
    vertices.push_back(0.1f);
    vertices.push_back(0.1f);
    
    vertices.push_back(0.0f);
    vertices.push_back(-axisLen);
    vertices.push_back(0.0f);
    vertices.push_back(0.1f);
    vertices.push_back(0.8f);
    vertices.push_back(0.1f);
    
    vertices.push_back(0.0f);
    vertices.push_back(axisLen);
    vertices.push_back(0.0f);
    vertices.push_back(0.1f);
    vertices.push_back(0.8f);
    vertices.push_back(0.1f);
    
    vertices.push_back(0.0f);
    vertices.push_back(0.0f);
    vertices.push_back(-axisLen);
    vertices.push_back(0.1f);
    vertices.push_back(0.1f);
    vertices.push_back(0.8f);
    
    vertices.push_back(0.0f);
    vertices.push_back(0.0f);
    vertices.push_back(axisLen);
    vertices.push_back(0.1f);
    vertices.push_back(0.1f);
    vertices.push_back(0.8f);
    
    GLuint VAO, VBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(float), vertices.data(), GL_STATIC_DRAW);
    
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)(3*sizeof(float)));
    glEnableVertexAttribArray(1);
    
    glBindVertexArray(0);
    return VAO;
}

static GLuint createSphereVAO(float radius, int segments) {
    std::vector<float> vertices;
    std::vector<unsigned int> indices;
    
    for(int i=0; i<=segments; i++) {
        float theta = i * 3.14159f / segments;
        for(int j=0; j<=segments; j++) {
            float phi = j * 2 * 3.14159f / segments;
            float x = radius * sin(theta) * cos(phi);
            float y = radius * cos(theta);
            float z = radius * sin(theta) * sin(phi);
            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);
            
            float nx = sin(theta) * cos(phi);
            float ny = cos(theta);
            float nz = sin(theta) * sin(phi);
            vertices.push_back(nx);
            vertices.push_back(ny);
            vertices.push_back(nz);
        }
    }
    
    for(int i=0; i<segments; i++) {
        for(int j=0; j<segments; j++) {
            int a = i*(segments+1) + j;
            int b = i*(segments+1) + j+1;
            int c = (i+1)*(segments+1) + j;
            int d = (i+1)*(segments+1) + j+1;
            indices.push_back(a);
            indices.push_back(b);
            indices.push_back(c);
            indices.push_back(b);
            indices.push_back(d);
            indices.push_back(c);
        }
    }
    
    GLuint VAO, VBO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(float), vertices.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size()*sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);
    
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)(3*sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    
    return VAO;
}

static GLuint createNucleusVAO(float radius, int segments) {
    std::vector<float> vertices;
    std::vector<unsigned int> indices;
    
    for(int i=0; i<=segments; i++) {
        float theta = i * 3.14159f / segments;
        for(int j=0; j<=segments; j++) {
            float phi = j * 2 * 3.14159f / segments;
            float x = radius * sin(theta) * cos(phi);
            float y = radius * cos(theta);
            float z = radius * sin(theta) * sin(phi);
            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);
            
            float nx = sin(theta) * cos(phi);
            float ny = cos(theta);
            float nz = sin(theta) * sin(phi);
            vertices.push_back(nx);
            vertices.push_back(ny);
            vertices.push_back(nz);
        }
    }
    
    for(int i=0; i<segments; i++) {
        for(int j=0; j<segments; j++) {
            int a = i*(segments+1) + j;
            int b = i*(segments+1) + j+1;
            int c = (i+1)*(segments+1) + j;
            int d = (i+1)*(segments+1) + j+1;
            indices.push_back(a);
            indices.push_back(b);
            indices.push_back(c);
            indices.push_back(b);
            indices.push_back(d);
            indices.push_back(c);
        }
    }
    
    GLuint VAO, VBO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(float), vertices.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size()*sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);
    
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)(3*sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    
    return VAO;
}

static GLuint createPointVAO() {
    GLuint VAO, VBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)(3*sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    return VAO;
}

static void framebufferSizeCallback(GLFWwindow*, int width, int height) {
    glViewport(0, 0, width, height);
}

static void errorCallback(int error, const char* description) {
    std::cerr << "GLFW Error [" << error << "]: " << description << std::endl;
}

struct UserData {
    Camera cam;
    bool dragging;
    int currentZ;
};

static void mouseButtonCallback(GLFWwindow* window, int button, int action, int) {
    UserData* ud = (UserData*)glfwGetWindowUserPointer(window);
    if(button == GLFW_MOUSE_BUTTON_LEFT) {
        if(action == GLFW_PRESS) ud->dragging = true;
        else if(action == GLFW_RELEASE) ud->dragging = false;
    }
}

static void cursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
    static double lastX = 0, lastY = 0;
    UserData* ud = (UserData*)glfwGetWindowUserPointer(window);
    if(!ud->dragging) { lastX = xpos; lastY = ypos; return; }
    double dx = xpos - lastX;
    double dy = ypos - lastY;
    lastX = xpos;
    lastY = ypos;
    ud->cam.theta += dy * 0.005f;
    ud->cam.phi += dx * 0.005f;
    if(ud->cam.theta < 0.1f) ud->cam.theta = 0.1f;
    if(ud->cam.theta > 3.0f) ud->cam.theta = 3.0f;
}

static void scrollCallback(GLFWwindow* window, double, double yoffset) {
    UserData* ud = (UserData*)glfwGetWindowUserPointer(window);
    ud->cam.distance -= yoffset * 0.5f;
    if(ud->cam.distance < 1.0f) ud->cam.distance = 1.0f;
    if(ud->cam.distance > 20.0f) ud->cam.distance = 20.0f;
}

static void keyCallback(GLFWwindow* window, int key, int, int action, int) {
    if(action != GLFW_PRESS) return;
    UserData* ud = (UserData*)glfwGetWindowUserPointer(window);
    if(key == GLFW_KEY_RIGHT || key == GLFW_KEY_UP) {
        if(ud->currentZ < 36) { ud->currentZ++; }
    }
    else if(key == GLFW_KEY_LEFT || key == GLFW_KEY_DOWN) {
        if(ud->currentZ > 1) { ud->currentZ--; }
    }
}

static void multiplyMat4(float* a, float* b, float* result) {
    for(int i=0; i<4; i++) {
        for(int j=0; j<4; j++) {
            result[i*4+j] = 0;
            for(int k=0; k<4; k++) {
                result[i*4+j] += a[i*4+k] * b[k*4+j];
            }
        }
    }
}

static void setTranslationMatrix(float* mat, float x, float y, float z) {
    mat[0] = 1; mat[1] = 0; mat[2] = 0; mat[3] = 0;
    mat[4] = 0; mat[5] = 1; mat[6] = 0; mat[7] = 0;
    mat[8] = 0; mat[9] = 0; mat[10] = 1; mat[11] = 0;
    mat[12] = x; mat[13] = y; mat[14] = z; mat[15] = 1;
}

static void setScaleMatrix(float* mat, float sx, float sy, float sz) {
    mat[0] = sx; mat[1] = 0; mat[2] = 0; mat[3] = 0;
    mat[4] = 0; mat[5] = sy; mat[6] = 0; mat[7] = 0;
    mat[8] = 0; mat[9] = 0; mat[10] = sz; mat[11] = 0;
    mat[12] = 0; mat[13] = 0; mat[14] = 0; mat[15] = 1;
}

int main() {
    glfwSetErrorCallback(errorCallback);
    if(!glfwInit()) {
        std::cerr << "Error initializing GLFW" << std::endl;
        return -1;
    }
    
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_SAMPLES, 4);
    
    GLFWwindow* window = glfwCreateWindow(1024, 768, "Atom Simulator 3D", nullptr, nullptr);
    if(!window) {
        std::cerr << "Error creating window" << std::endl;
        glfwTerminate();
        return -1;
    }
    
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);
    glfwSwapInterval(1);
    
    if(!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Error initializing GLAD" << std::endl;
        return -1;
    }
    
    std::cout << "OpenGL: " << glGetString(GL_VERSION) << std::endl;
    std::cout << "Renderer: " << glGetString(GL_RENDERER) << std::endl;
    std::cout << "Controls: Arrow Up/Down = change element, Mouse drag = rotate, Scroll = zoom" << std::endl;
    
    const char* vsSource = "#version 330 core\n"
        "layout(location=0) in vec3 aPos;\n"
        "layout(location=1) in vec3 aCol;\n"
        "uniform mat4 uModelView;\n"
        "uniform mat4 uProjection;\n"
        "out vec3 vCol;\n"
        "void main() {\n"
        "    gl_Position = uProjection * uModelView * vec4(aPos, 1.0);\n"
        "    vCol = aCol;\n"
        "}\n";
    
    const char* fsSource = "#version 330 core\n"
        "in vec3 vCol;\n"
        "out vec4 fCol;\n"
        "void main() { fCol = vec4(vCol, 1.0); }\n";
    
    GLuint vs = createShader(vsSource, GL_VERTEX_SHADER);
    GLuint fs = createShader(fsSource, GL_FRAGMENT_SHADER);
    GLuint program = glCreateProgram();
    glAttachShader(program, vs);
    glAttachShader(program, fs);
    glLinkProgram(program);
    glDeleteShader(vs);
    glDeleteShader(fs);
    glUseProgram(program);
    
    GLuint modelViewLoc = glGetUniformLocation(program, "uModelView");
    GLuint projLoc = glGetUniformLocation(program, "uProjection");
    
    UserData userData;
    userData.currentZ = 6;
    userData.cam.theta = 1.2f;
    userData.cam.phi = 0.8f;
    userData.cam.distance = 8.0f;
    userData.cam.fov = 60.0f;
    userData.dragging = false;
    
    glfwSetWindowUserPointer(window, &userData);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetCursorPosCallback(window, cursorPosCallback);
    glfwSetScrollCallback(window, scrollCallback);
    glfwSetKeyCallback(window, keyCallback);
    
    Atom atom = createAtom(userData.currentZ, userData.currentZ);
    std::vector<OrbitalPoint> orbitalPoints = generateOrbitalPoints(atom, 20000);
    
    GLuint gridVAO = createGridVAO(12.0f, 12);
    GLuint nucleusVAO = createNucleusVAO(0.3f, 16);
    GLuint pointVAO = createPointVAO();
    
    std::vector<Vec3> electronPos;
    std::vector<Color> electronCol;
    for(const Electron& e : atom.electrons) {
        Color c = getOrbitalColor(e.l, e.ml);
        float bright = 0.7f + 0.3f * (1.0f / e.n);
        c.r *= bright; c.g *= bright; c.b *= bright;
        electronPos.push_back(e.pos);
        electronCol.push_back(c);
    }
    
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor(0.02f, 0.02f, 0.05f, 1.0f);
    glPointSize(2.5f);
    
    auto lastTime = std::chrono::steady_clock::now();
    float totalTime = 0.0f;
    
    while(!glfwWindowShouldClose(window)) {
        auto currentTime = std::chrono::steady_clock::now();
        float dt = std::chrono::duration<float>(currentTime - lastTime).count();
        lastTime = currentTime;
        totalTime += dt;
        if(dt > 0.05f) dt = 0.05f;
        
        glfwPollEvents();
        
        if(userData.currentZ != atom.Z) {
            atom = createAtom(userData.currentZ, userData.currentZ);
            orbitalPoints = generateOrbitalPoints(atom, 20000);
            
            electronPos.clear(); electronCol.clear();
            for(const Electron& e : atom.electrons) {
                Color c = getOrbitalColor(e.l, e.ml);
                float bright = 0.7f + 0.3f * (1.0f / e.n);
                c.r *= bright; c.g *= bright; c.b *= bright;
                electronPos.push_back(e.pos);
                electronCol.push_back(c);
            }
            
            pointVAO = createPointVAO();
        }
        
        float nucleusExclusionRadius = atom.nuclearRadius * 1.8f;
        
        for(Electron& e : atom.electrons) {
            Vec3 force(0,0,0);
            float k = 0.5f;
            
            Vec3 toNucleus = e.pos * -1.0f;
            float distN = toNucleus.magnitude() + 0.1f;
            
            if(distN < nucleusExclusionRadius) {
                Vec3 pushOut = e.pos.normalize() * (nucleusExclusionRadius - distN) * 2.0f;
                e.pos = e.pos + pushOut;
            }
            
            float coulomb = k * atom.Z / (distN * distN);
            force = force + toNucleus.normalize() * (-coulomb);
            
            for(const Electron& other : atom.electrons) {
                if(&e == &other) continue;
                Vec3 diff = e.pos - other.pos;
                float d = diff.magnitude() + 0.1f;
                float rep = k / (d * d);
                force = force + diff.normalize() * rep * 0.3f;
            }
            
            e.vel = e.vel + force * dt * 2.0f;
            e.vel = e.vel * 0.99f;
            e.pos = e.pos + e.vel * dt;
            
            float distAfter = e.pos.magnitude();
            if(distAfter < nucleusExclusionRadius) {
                e.pos = e.pos.normalize() * nucleusExclusionRadius;
                Vec3 velNormal = e.vel.normalize();
                Vec3 posNormal = e.pos.normalize();
                float dotProduct = velNormal.dot(posNormal);
                if(dotProduct < 0) {
                    e.vel = e.vel - posNormal * (velNormal.dot(posNormal)) * 2.0f;
                }
            }
            
            float maxR = e.n * e.n * 1.2f;
            if(e.pos.magnitude() > maxR) {
                e.vel = e.vel - e.pos.normalize() * e.vel.dot(e.pos.normalize()) * 0.5f;
                e.pos = e.pos.normalize() * maxR * 0.95f;
            }
        }
        
        for(size_t i=0; i<electronPos.size() && i<atom.electrons.size(); i++) {
            electronPos[i] = atom.electrons[i].pos;
        }
        
        int pointsPerOrbital = 20000 / atom.electrons.size();
        if(pointsPerOrbital < 100) pointsPerOrbital = 100;
        orbitalPoints = generateOrbitalPoints(atom, pointsPerOrbital * atom.electrons.size());
        
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        float aspect = (float)width / (float)height;
        
        Vec3 eye = sphericalToCart(userData.cam.distance, userData.cam.theta, userData.cam.phi);
        Vec3 center(0,0,0);
        Vec3 up(0,1,0);
        
        float modelView[16] = {0};
        
        Vec3 f = (center - eye).normalize();
        Vec3 s = f.cross(up).normalize();
        Vec3 u = s.cross(f);
        
        modelView[0] = s.x; modelView[1] = u.x; modelView[2] = -f.x; modelView[3] = 0;
        modelView[4] = s.y; modelView[5] = u.y; modelView[6] = -f.y; modelView[7] = 0;
        modelView[8] = s.z; modelView[9] = u.z; modelView[10] = -f.z; modelView[11] = 0;
        modelView[12] = -s.dot(eye); modelView[13] = -u.dot(eye); modelView[14] = f.dot(eye); modelView[15] = 1;
        
        float fovRad = userData.cam.fov * 3.14159f / 180.0f;
        float tanHalf = tan(fovRad / 2.0f);
        float near = 0.1f;
        float far = 50.0f;
        float proj[16] = {0};
        proj[0] = 1.0f / (aspect * tanHalf);
        proj[5] = 1.0f / tanHalf;
        proj[10] = -(far + near) / (far - near);
        proj[11] = -1.0f;
        proj[14] = -(2.0f * far * near) / (far - near);
        proj[15] = 0.0f;
        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        glUniformMatrix4fv(modelViewLoc, 1, GL_FALSE, modelView);
        glUniformMatrix4fv(projLoc, 1, GL_FALSE, proj);
        
        glBindVertexArray(gridVAO);
        int gridLines = (12 * 2 + 2) * 4;
        glDrawArrays(GL_LINES, 0, gridLines);
        
        float nucleusScale = atom.nuclearRadius * 1.5f;
        float scaleMat[16];
        setScaleMatrix(scaleMat, nucleusScale, nucleusScale, nucleusScale);
        
        float nucleusMV[16];
        multiplyMat4(modelView, scaleMat, nucleusMV);
        
        glUniformMatrix4fv(modelViewLoc, 1, GL_FALSE, nucleusMV);
        
        glBindVertexArray(nucleusVAO);
        glDrawElements(GL_TRIANGLES, 16*16*6, GL_UNSIGNED_INT, 0);
        
        std::vector<float> pointData;
        pointData.reserve(orbitalPoints.size() * 6);
        for(size_t i=0; i<orbitalPoints.size(); i++) {
            pointData.push_back(orbitalPoints[i].pos.x);
            pointData.push_back(orbitalPoints[i].pos.y);
            pointData.push_back(orbitalPoints[i].pos.z);
            pointData.push_back(orbitalPoints[i].color.r * 0.8f);
            pointData.push_back(orbitalPoints[i].color.g * 0.8f);
            pointData.push_back(orbitalPoints[i].color.b * 0.8f);
        }
        
        glUniformMatrix4fv(modelViewLoc, 1, GL_FALSE, modelView);
        glBindVertexArray(pointVAO);
        GLuint buf;
        glGenBuffers(1, &buf);
        glBindBuffer(GL_ARRAY_BUFFER, buf);
        glBufferData(GL_ARRAY_BUFFER, pointData.size()*sizeof(float), pointData.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)(3*sizeof(float)));
        glEnableVertexAttribArray(1);
        glDrawArrays(GL_POINTS, 0, orbitalPoints.size());
        glDeleteBuffers(1, &buf);
        
        glfwSwapBuffers(window);
    }
    
    glDeleteProgram(program);
    glDeleteVertexArrays(1, &gridVAO);
    glDeleteVertexArrays(1, &nucleusVAO);
    glDeleteVertexArrays(1, &pointVAO);
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}