#version 330 core
layout (location = 0) out vec4 normalDepth_c;
layout (location = 1) out vec3 diffuse_c;
in float depth_mv;
in vec3 normal_mv;
uniform vec3 diffuse;
void main() {
    normalDepth_c = vec4(normalize(normal_mv), depth_mv);
    diffuse_c = diffuse;
}
