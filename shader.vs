#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;

out vec3 ourColor;

struct UniformBufferObject {
    mat4 model;
    mat4 view;
    mat4 proj;
};

uniform UniformBufferObject ubo;

void main()
{
    gl_Position = ubo.proj *ubo.view *  ubo.model * vec4(aPos, 1.0);
    gl_PointSize = 10;
    // gl_Position = ubo.model * vec4(aPos, 1.0);
    // gl_Position = vec4(aPos, 1.0);
    // vec4 res = ubo.model * vec4(aPos, 1.0)
    // ourColor = vec3(aColor.x, aColor.y, aColor.z) ;
    // ourColor = vec3(ubo.model[0][0], 1, 0) ;
    ourColor = aColor;
}