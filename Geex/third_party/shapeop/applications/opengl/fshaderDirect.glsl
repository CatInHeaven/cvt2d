#version 330 core
layout (location = 0) out vec3 color;
in vec2 uv;
uniform vec3 light;
uniform mat4 inverseProjection;
uniform mat4 lightMatrix;
uniform sampler2D normalDepthTex;
uniform sampler2D diffuseTex;
uniform sampler2DShadow shadowTex;
uniform sampler2D randomTex;

vec3 positionFromDepth(in vec2 uv,
                       in float depth) {
     vec3 ray = mat3(inverseProjection)*vec3((uv-0.5)*2.0, 1.0);
     return vec3(ray.xy*depth, -depth);
}
vec3 directLight(in vec3 position, in vec3 normal, in vec3 diffuse) {
  float a = 0.4;
  float b = 0.0;
  float lambert = max(dot(normal, normalize(light - position)), 0.0);
  return (1.0-b)*(a*diffuse + (1.0-a)*diffuse*lambert) + b*vec3(1.0);
}
void main() {
    vec3 diffuse1 = texture(diffuseTex, uv).xyz;
    vec4 normalDepth1 = texture(normalDepthTex, uv);
    vec3 normal1 = normalDepth1.xyz;
    vec3 position1 = positionFromDepth(uv, normalDepth1.w);
    float s = 1.0;
    vec2 inc = s*vec2(1.0/1024.0, 1.0/768.0);
    float mu = 5.0/128.0;
    vec2 fres = mu*vec2(1024.0, 768.0);
    vec2 random = texture(randomTex, uv*fres).xy;
    random = inc*(random*2.0-vec2(1.0));
    vec2 p = inc;
    float bias = 0.0004;
    vec4 shadowCoord = lightMatrix * vec4(position1, 1);
    float cd = shadowCoord.z/shadowCoord.w;
    vec3 uv2 = vec3(shadowCoord.xy/shadowCoord.w, cd) - vec3(0.0, 0.0, bias);
    float visibility = 0.0;
    if(uv2.x >= 0.001 && uv2.x <= 0.999 && uv2.y >= 0.001 && uv2.y <= 0.999) {
      visibility += texture(shadowTex, uv2);
      for(float i = 0.0; i < 6.0; ++i) {
        vec2 np =  (p+random.xy)/cd;
        visibility += texture(shadowTex, uv2 + vec3(np.x, np.y, 0.0));
        visibility += texture(shadowTex, uv2 + vec3(np.x, -np.y, 0.0));
        visibility += texture(shadowTex, uv2 + vec3(-np.x, np.y, 0.0));
        visibility += texture(shadowTex, uv2 + vec3(-np.x, -np.y, 0.0));
        visibility += texture(shadowTex, uv2 + vec3(0.0, np.y, 0.0));
        visibility += texture(shadowTex, uv2 + vec3(0.0, -np.y, 0.0));
        visibility += texture(shadowTex, uv2 + vec3(np.x, 0.0, 0.0));
        visibility += texture(shadowTex, uv2 + vec3(-np.x, 0.0, 0.0));
        p += inc;
      }
      visibility = 1.0 - 0.2*(1.0-visibility/49.0);
    } else {
      visibility = 1.0;
    }

    color = visibility * directLight(position1, normal1, diffuse1);
}

