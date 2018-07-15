
/*
	About mapping a C array into an Eigen vector: http://eigen.tuxfamily.org/dox/classEigen_1_1Map.html
*/

#include "emscripten.h"
#include <cmath>
#include <cstdlib>
#include <GLES2/gl2.h>
#include <EGL/egl.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>


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


/*
	Main application
*/

// Properties
double now, delta;
GLuint program, vbo;
int size;
float * data;
Eigen::Map<Eigen::VectorXf> vector(NULL, 0);
Eigen::SparseMatrix<float> matrix;
Eigen::SimplicialLLT<Eigen::SparseMatrix<float> > cholesky;

// Set info text
void setInfo(float fps, float error) {
	EM_ASM_ARGS({
		info.innerHTML = $0 + " values, " + $1 + " fps, " + $2 + " average error";
	}, size, fps, error);
}

// Called each frame
void tick() {
	// Create random vector
	for (int i = 0; i < size; ++i)
		vector[i] = emscripten_random();
	// Solve system with new vector
	Eigen::VectorXf result = cholesky.solve(vector);
	float dist = (matrix * result - vector).norm();
	// Update time
	static int frame = 0;
	static double last = 0;
	static float distsum = 0;
	++frame;
	distsum += dist;
	double time = emscripten_get_now() * 0.001;
	delta = time - now;
	now = time;
	if (now - last >= 5) {
		double fps = frame / (now - last);
		float distavg = distsum / frame;
		emscripten_log(EM_LOG_CONSOLE, "%f fps, %f average error", fps, distavg);
		setInfo(fps, distavg);
		frame = 0;
		last = now;
		distsum = 0;
	}
	// Render vector as dummy triangles
	glClear(GL_COLOR_BUFFER_BIT);
	glBufferData(GL_ARRAY_BUFFER, size * sizeof(float), data, GL_STREAM_DRAW);
	glDrawArrays(GL_TRIANGLES, 0, size / 3);
}

// Entry point called by hand-written JS
extern "C" bool init(int n) {
	// Initialize OpenGL ES
	if (!createContext())
		return false;
	GLchar const * vs =
		"precision mediump float;\n"
		"attribute vec3 pos;\n"
		"varying float v_hue;\n"
		"void main() {\n"
		"	gl_Position = vec4(pos.x * 2.0 - 1.0, pos.y * 2.0 - 1.0, 0.0, 1.0);\n"
		"	v_hue = pos.z;\n"
		"}\n";
	GLchar const * fs =
		"precision mediump float;\n"
		"varying float v_hue;\n"
		"void main() {\n"
		"	if (v_hue < 1.0 / 6.0)"
		"		gl_FragColor = vec4(1.0, v_hue * 6.0, 0.0, 1.0);\n"
		"	else if (v_hue < 2.0 / 6.0)"
		"		gl_FragColor = vec4(2.0 - v_hue * 6.0, 1.0, 0.0, 1.0);\n"
		"	else if (v_hue < 3.0 / 6.0)"
		"		gl_FragColor = vec4(0.0, 1.0, v_hue * 6.0 - 2.0, 1.0);\n"
		"	else if (v_hue < 4.0 / 6.0)"
		"		gl_FragColor = vec4(0.0, 4.0 - v_hue * 6.0, 1.0, 1.0);\n"
		"	else if (v_hue < 5.0 / 6.0)"
		"		gl_FragColor = vec4(v_hue * 6.0 - 4.0, 0.0, 1.0, 1.0);\n"
		"	else"
		"		gl_FragColor = vec4(1.0, 0.0, 6.0 - v_hue * 6.0, 1.0);\n"
		"}\n";
	GLchar const * attribs[] = {"pos"};
	program = createProgram(vs, fs, attribs, 1);
	glUseProgram(program);
	size = n;
	data = new float[size];
	new (&vector) Eigen::Map<Eigen::VectorXf>(data, size); // Well, this seems to be the legit way...
	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void *)0);
	glClearColor(0.6f, 0.6f, 0.6f, 1.0f);
	setInfo(0, 0);
	// Build matrix and factorize
	Eigen::Triplet<float> * triplets = new Eigen::Triplet<float>[8 * size];
	for (int k = 0; k < 8 * size; ++k) {
		int i = rand() % size;
		int j = rand() % size;
		if (j > i) {
			int t = i;
			i = j;
			j = t;
		}
		float v = emscripten_random() * 10;
		triplets[k] = Eigen::Triplet<float>(i, j, v);
	}
	Eigen::SparseMatrix<float> sub(size, size);
	sub.setFromTriplets(triplets, triplets + 8 * size);
	delete triplets;
	matrix.resize(size, size);
	matrix.setIdentity();
	matrix += sub * sub.transpose();
	int count = 0;
	for (int k = 0; k < matrix.outerSize(); ++k)
		for (Eigen::SparseMatrix<float>::InnerIterator it(matrix, k); it; ++it)
			if (it.value() >= 0.000001f || it.value() <= -0.000001f)
				++count;
	emscripten_log(EM_LOG_CONSOLE, "%i non-zero values", count);
	double before = emscripten_get_now() * 0.001;
	cholesky.compute(matrix);
	now = emscripten_get_now() * 0.001;
	emscripten_log(EM_LOG_CONSOLE, "Factorized in %f seconds", (float)(now - before));
	// Ready, start main loop
	emscripten_log(EM_LOG_CONSOLE, "Loaded");
	emscripten_set_main_loop(tick, 0, 0);
	return true;
}
