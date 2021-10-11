
#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;

out vec3 ourColor;

struct UniformBufferObject
{
    mat4 model;
    mat4 view;
    mat4 proj;
};

uniform UniformBufferObject ubo;

uniform vec3 ball_color;
void main()
{
    vec3 aColor = ball_color;
    if (length(aNormal) > 1e-3)
    {
        vec3 lightPos = vec3(0, 10, 0);
        vec4 vertPos4 = ubo.view * ubo.model * vec4(aPos, 1.0);
        vec3 vertPos = vec3(vertPos4) / vertPos4.w;
        vec3 normalInterp = aNormal;
        vec3 N = normalize(normalInterp);
        vec3 L = normalize(lightPos - vertPos);
        float lambertian = max(dot(N, L), 0.0);
        float specular = 0.0f;
        float shininessVal = 80;
        if (lambertian > 0.0)
        {
            vec3 R = reflect(-L, N);      // Reflected light vector
            vec3 V = normalize(-vertPos); // Vector to viewer
            // Compute the specular term
            float specAngle = max(dot(R, V), 0.0);

            specular = pow(specAngle, shininessVal);
        }
        vec3 diffuseColor = vec3(1.0, 0.5, 0.5);
        vec3 specularColor = vec3(1, 1, 1);
        vec3 ambientColor = vec3(aColor.x, aColor.y, aColor.z);

        float Ka = 1;
        float Kd = 1;
        float Ks = 1;
        vec3 color = vec3(Ka * ambientColor +
                            Kd * lambertian * diffuseColor +
                            Ks * specular * specularColor);

        ourColor = color;
    }
    else
    {
        ourColor = aColor;
    }
    gl_Position = ubo.proj * ubo.view * ubo.model * vec4(aPos, 1.0);
    // ourNormal = aNormal;

    // gl_Position = ubo.model * vec4(aPos, 1.0);
    // gl_Position = vec4(aPos, 1.0);
    // vec4 res = ubo.model * vec4(aPos, 1.0)
    // ourColor = vec3(aColor.x, aColor.y, aColor.z) ;
    // ourColor = vec3(ubo.model[0][0], 1, 0) ;
}