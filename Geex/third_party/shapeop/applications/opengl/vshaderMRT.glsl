#version 330 core 
layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
out float depth_mv;
out vec3 normal_mv;
void main() {
    mat4 model_view = view * model;
    normal_mv = mat3(model_view) * normal;
    vec3 position_mv = (model_view * vec4(position , 1.0)).xyz;
    depth_mv = -position_mv.z;
    gl_Position = projection * vec4(position_mv, 1.0);
}


