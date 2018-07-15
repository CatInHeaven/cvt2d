#version 330 core 
layout(location = 0) in vec3 position;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
void main() {
    mat4 model_view = view * model;
    vec3 position_mv = (model_view * vec4(position , 1.0)).xyz;
    gl_Position = projection * vec4(position_mv, 1.0);
}


