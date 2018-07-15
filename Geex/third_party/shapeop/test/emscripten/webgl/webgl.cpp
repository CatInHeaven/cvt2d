
#include "emscripten.h"
#include <cmath>
#include <cstdlib>
#include <GLES2/gl2.h>
#include <EGL/egl.h>

#define PI 3.14159265359


/*
	OpenGL ES helpers
*/

// Create an OpenGL ES context
bool createContext() {
	// Initialize OpenGL ES and get display
	EGLDisplay display = eglGetDisplay(EGL_DEFAULT_DISPLAY);
	if (display == EGL_NO_DISPLAY || !eglInitialize(display, NULL, NULL)) {
		emscripten_log(EM_LOG_CONSOLE | EM_LOG_ERROR, "Failed to initialize OpenGL ES");
		return false;
	}
	// Get config
	EGLint numConfigs;
	EGLint attribs[] = {
		EGL_RED_SIZE, 5,
		EGL_GREEN_SIZE, 6,
		EGL_BLUE_SIZE, 5,
		EGL_NONE
	};
	EGLConfig config;
	if (!eglGetConfigs(display, NULL, 0, &numConfigs) || !eglChooseConfig(display, attribs, &config, 1, &numConfigs)) {
		emscripten_log(EM_LOG_CONSOLE | EM_LOG_ERROR, "Failed to get config");
		return false;
	}
	// Create surface
	EGLNativeWindowType window;
	EGLSurface surface = eglCreateWindowSurface(display, config, window, NULL);
    if (surface == EGL_NO_SURFACE) {
		emscripten_log(EM_LOG_CONSOLE | EM_LOG_ERROR, "Failed to create surface");
		return false;
	}
	// Create context and make current
	EGLint contextAttribs[] =
	{
		EGL_CONTEXT_CLIENT_VERSION, 2,
		EGL_NONE
	};
	EGLContext context = eglCreateContext(display, config, EGL_NO_CONTEXT, contextAttribs);
    if (context == EGL_NO_CONTEXT || !eglMakeCurrent(display, surface, surface, context)) {
		emscripten_log(EM_LOG_CONSOLE | EM_LOG_ERROR, "Failed to create context");
		return false;
	}
	return true;
}

// Create a shader object of given type (GL_VERTEX_SHADER or GL_FRAGMENT_SHADER) with specified source code
GLuint createShader(GLenum type, GLchar const * src) {
	GLuint id = glCreateShader(type);
	glShaderSource(id, 1, &src, NULL);
	glCompileShader(id);
	GLint compiled;
	glGetShaderiv(id, GL_COMPILE_STATUS, &compiled);
	if (!compiled) {
		GLint length;
		glGetShaderiv(id, GL_INFO_LOG_LENGTH, &length);
		GLchar * buffer = new GLchar[length];
		glGetShaderInfoLog(id, length, NULL, buffer);
		emscripten_log(EM_LOG_CONSOLE | EM_LOG_ERROR, buffer);
		delete[] buffer;
		glDeleteShader(id);
		return 0;
	}
	return id;
}

// Create a program object using source codes and specified attributes.
GLuint createProgram(GLchar const * vs, GLchar const * fs, GLchar const * const * attribs, int attribsCount) {
	// Create shader objects
	GLuint v = createShader(GL_VERTEX_SHADER, vs);
	GLuint f = createShader(GL_FRAGMENT_SHADER, fs);
	if (!v || !f) {
		emscripten_log(EM_LOG_CONSOLE | EM_LOG_ERROR, "Failed to compile shader!\n");
		glDeleteShader(v);
		glDeleteShader(f);
		return 0;
	}
	// Create program object
	GLuint id = glCreateProgram();
	glAttachShader(id, v);
	glAttachShader(id, f);
	glDeleteShader(v);
	glDeleteShader(f);
	for (unsigned i = 0; i < attribsCount; ++i)
		glBindAttribLocation(id, i, attribs[i]);
	glLinkProgram(id);
	GLint linked;
	glGetProgramiv(id, GL_LINK_STATUS, &linked);
	if (!linked) {
		GLint length;
		glGetProgramiv(id, GL_INFO_LOG_LENGTH, &length);
		GLchar* buffer = new GLchar[length];
		glGetProgramInfoLog(id, length, NULL, buffer);
		emscripten_log(EM_LOG_CONSOLE | EM_LOG_ERROR, buffer);
		delete[] buffer;
		glDeleteProgram(id);
		return 0;
	}
	return id;
}

// Load identity in given array
void matrixIdentity(GLfloat * m) {
	m[ 0] = 1; m[ 4] = 0; m[ 8] = 0; m[12] = 0;
	m[ 1] = 0; m[ 5] = 1; m[ 9] = 0; m[13] = 0;
	m[ 2] = 0; m[ 6] = 0; m[10] = 1; m[14] = 0;
	m[ 3] = 0; m[ 7] = 0; m[11] = 0; m[15] = 1;
}

// Load conical projection matrix in given array
void matrixConical(GLfloat * m, float width, float height, float near, float far) {
	float r = width * 0.5f, t = height * 0.5f;
	m[ 0] = near / r; m[ 4] =        0; m[ 8] =                            0; m[12] =                              0;
	m[ 1] =        0; m[ 5] = near / t; m[ 9] =                            0; m[13] =                              0;
	m[ 2] =        0; m[ 6] =        0; m[10] = -(far + near) / (far - near); m[14] = -2 * far * near / (far - near);
	m[ 3] =        0; m[ 7] =        0; m[11] =                           -1; m[15] =                              0;
}

// Add translation effect to given array (applied after current transform)
void matrixTranslate(GLfloat * m, float x, float y, float z) {
	m[12] += x;
	m[13] += y;
	m[14] += z;
}

// Multiply matrix by given matrix (m = n * m)
void matrixMultiply(GLfloat * m, GLfloat const * n) {
	GLfloat b[16];
	for (unsigned i = 0; i < 16; ++i) {
		b[i] = m[i];
		m[i] = 0;
	}
	for (unsigned i = 0; i < 4; ++i)
		for (unsigned j = 0; j < 4; ++j)
			for (unsigned k = 0; k < 4; ++k)
				m[i + j * 4] += n[i + k * 4] * b[k + j * 4];
}

// Add rotation effect around base axis to given array (applied after current transform)
void matrixRotateX(GLfloat * m, float a) {
	float c = cos(a), s = sin(a);
	GLfloat n[16] = {1, 0, 0, 0, 0, c, s, 0, 0, -s, c, 0, 0, 0, 0, 1};
	matrixMultiply(m, n);
}
void matrixRotateY(GLfloat * m, float a) {
	float c = cos(a), s = sin(a);
	GLfloat n[16] = {c, 0, -s, 0, 0, 1, 0, 0, s, 0, c, 0, 0, 0, 0, 1};
	matrixMultiply(m, n);
}
void matrixRotateZ(GLfloat * m, float a) {
	float c = cos(a), s = sin(a);
	GLfloat n[16] = {c, s, 0, 0, -s, c, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
	matrixMultiply(m, n);
}


/*
	Main application
*/

// Time
double now;
double delta;

// Shader
GLint renderProgram;
GLuint renderProjection;
GLuint renderWorld;

// Mesh
GLuint meshVerticesObject;
int meshVertices;

// Input
float viewX;
float viewY;

// Upload indices and vertices to GPU
void updateVertices(GLfloat const * array, int count) {
	glBindBuffer(GL_ARRAY_BUFFER, meshVerticesObject);
	glBufferData(GL_ARRAY_BUFFER, count * sizeof(GLfloat), array, GL_DYNAMIC_DRAW);
	meshVertices = count / 6;
}

// JS interface to load model
int loadCount, loadIndex;
GLfloat * loadBuffer;
extern "C" void loadBegin(int count) {
	loadCount = count;
	loadIndex = 0;
	loadBuffer = new GLfloat[count * 6];
}
extern "C" void loadVertex(float x, float y, float z, float xn, float yn, float zn) {
	if (isnan(x) || isnan(y) || isnan(z)) {
		x = y = z = 0;
		emscripten_log(EM_LOG_CONSOLE | EM_LOG_WARN, "Vertex is NaN");
	}
	loadBuffer[loadIndex++] = x;
	loadBuffer[loadIndex++] = y;
	loadBuffer[loadIndex++] = z;
	float n = sqrt(xn * xn + yn * yn + zn * zn);
	if (n > 0.001f) {
		xn /= n;
		yn /= n;
		zn /= n;
	} else {
		xn = 1;
		yn = zn = 0;
	}
	loadBuffer[loadIndex++] = xn;
	loadBuffer[loadIndex++] = yn;
	loadBuffer[loadIndex++] = zn;
	// emscripten_log(EM_LOG_CONSOLE, "%f %f %f (%f %f %f)", x, y, z, xn, yn, zn);
}
extern "C" void loadEnd() {
	updateVertices(loadBuffer, loadIndex);
	delete loadBuffer;
}

// Mouse interaction interface
extern "C" void mouseMoved(float dx, float dy) {
	viewX += dx;
	viewY += dy;
}

// Called each frame
void tick() {
	// Update time
	static int frame = 0;
	static double last = 0;
	++frame;
	double time = emscripten_get_now() * 0.001;
	delta = time - now;
	now = time;
	if (now - last >= 5) {
		double fps = frame / (now - last);
		emscripten_log(EM_LOG_CONSOLE, "%f fps", fps);
		frame = 0;
		last = now;
	}
	// Render
	glEnable(GL_DEPTH_TEST);
	glClearColor(0.6f, 0.6f, 0.6f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glBindBuffer(GL_ARRAY_BUFFER, meshVerticesObject);
	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void *)0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void *)(3 * sizeof(GLfloat)));
	glUseProgram(renderProgram);
	float matrix[16];
	matrixIdentity(matrix);
	matrixRotateX(matrix, 0.2f);
	matrixRotateY(matrix, fmod(viewX, 2 * PI));
	matrixTranslate(matrix, 0, 0, viewY);
	glUniformMatrix4fv(renderWorld, 1, GL_FALSE, matrix);
	glDrawArrays(GL_TRIANGLES, 0, meshVertices);
	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
}

// Entry point called by hand-written JS
extern "C" bool init(int width, int height) {
	// Initialize OpenGL ES
	if (!createContext())
		return false;
	glGenBuffers(1, &meshVerticesObject);
	glBindBuffer(GL_ARRAY_BUFFER, meshVerticesObject);
	GLchar const * renderVs =
		"precision mediump float;\n"
		"attribute vec3 pos, nor;\n"
		"varying vec3 v_nor;\n"
		"uniform mat4 projection;\n"
		"uniform mat4 world;\n"
		"void main() {\n"
		"	gl_Position = projection * world * vec4(pos, 1.0);\n"
		"	v_nor = (world * vec4(nor, 0.0)).xyz;\n"
		"}\n";
	GLchar const * renderFs =
		"precision mediump float;\n"
		"varying vec3 v_nor;\n"
		"void main() {\n"
		"	gl_FragColor = vec4((v_nor + vec3(1.0, 1.0, 1.0)) * 0.5, 1.0);\n"
		"}\n";
	GLchar const * renderAttribs[] = {"pos", "nor"};
	renderProgram = createProgram(renderVs, renderFs, renderAttribs, 2);
	renderProjection = glGetUniformLocation(renderProgram, "projection");
	renderWorld = glGetUniformLocation(renderProgram, "world");
	// Set default values
	now = emscripten_get_now() * 0.001;
	glUseProgram(renderProgram);
	float matrix[16];
	matrixConical(matrix, 1.0f, 3.0f / 4.0f, 0.1f, 100.0f);
	glUniformMatrix4fv(renderProjection, 1, GL_FALSE, matrix);
	meshVertices = 0;
	viewX = 0;
	viewY = -4;
	// Ready, start main loop
	emscripten_log(EM_LOG_CONSOLE, "Loaded");
	emscripten_set_main_loop(tick, 0, 0);
	return true;
}
