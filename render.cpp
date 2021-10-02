#include "render.h"
#include "cameras/ArcBallCamera.h"
#include "RenderCallback.h"

cRender::cRender()
{
}
cRender::~cRender()
{

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    // glDeleteBuffers(1, &EBO);
    glDeleteProgram(shaderProgram);
}

void cRender::InitCam()
{
    tVector3f pos = tVector3f(1, 1, 1);
    tVector3f center = tVector3f(0, 0, 0);
    tVector3f up = tVector3f(0, 1, 0);
    float near_plane_dist = 1e-3;
    float far_plane_dist = 1e2;
    mCam = std::make_shared<cArcBallCamera>(pos,
                                            center,
                                            up,
                                            60,
                                            near_plane_dist,
                                            far_plane_dist);
}
void cRender::Init()
{
    InitGL();

    InitCam();
}

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

glm::mat4 E2GLM(const tMatrix4f &em)
{
    glm::mat4 mat;
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            mat[j][i] = em(i, j);
        }
    }
    return mat;
}
void cRender::Update()
{
    ourShader->use();

    // ourShader->setMat4("model", glm::mat4(1.0f));
    ourShader->setMat4("ubo.model", glm::mat4(1.0f));
    tMatrix4f eigen_view = this->mCam->ViewMatrix();
    // eigen_view.setIdentity();
    tMatrix4f eigen_proj = this->mCam->ProjMatrix(mWidth, mHeight, false);
    // eigen_proj.setIdentity();
    glm::mat4 glm_view = E2GLM(eigen_view);
    glm::mat4 glm_proj = E2GLM(eigen_proj);
    ourShader->setMat4("ubo.view", glm_view);
    ourShader->setMat4("ubo.proj", glm_proj);
    // std::cout << "cur proj = \n" << eigen_proj << std::endl;
    mCam->MouseRotate();
    glBindVertexArray(VAO); // seeing as we only have a single VAO there's no need to bind it every time, but we'll do so to keep things a bit more organized

    glDrawArrays(GL_TRIANGLES, 0, 3);
    //glDrawArrays(GL_TRIANGLES, 0, 6);
}

void cRender::InitGL()
{
    // init glfw
    if (!glfwInit())
    {
        std::cout << "[error] InitGLFW:: glfw inti failed" << std::endl;
        glfwTerminate();
    }
    glfwSetErrorCallback(ErrorCallback);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE); // fixed size
    mWindow = glfwCreateWindow(mWidth, mHeight, mWindowName.c_str(), NULL, NULL);
    if (mWindow == NULL)
    {
        std::cout << "[error] Failed to create GLFW window" << std::endl;
        glfwTerminate();
    }
    glfwSetWindowPos(mWindow, mStartX, mStartY);
    glfwMakeContextCurrent(mWindow);

    glfwSwapInterval(1); // enable vsync from ImGUI

    glfwSetKeyCallback(mWindow, KeyEventCallback);
    glfwSetCursorPosCallback(mWindow, MouseMoveEventCallback);
    glfwSetMouseButtonCallback(mWindow, MouseButtonEventCallback);
    glfwSetInputMode(mWindow, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    glfwSetFramebufferSizeCallback(mWindow, ResizeCallback);
    glfwSetScrollCallback(mWindow, ScrollCallback);
    if (GLEW_OK != glewInit())
    {
        std::cout << "[errpr] glew init failed " << std::endl;
        exit(1);
    }

    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_DEPTH_TEST);
    // glClearColor(0.2, 0.3, 0.4, 1);
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // link shaders
    ourShader = new Shader("../shader.vs", "../shader.fs");

    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------
    float vertices[] = {

        // positions         // colors
        0.5f, -0.5f, 0.0f, 1.0f, 0.0f, 0.0f,  // bottom right
        -0.5f, -0.5f, 0.0f, 0.0f, 1.0f, 0.0f, // bottom left
        0.0f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f    // top
    };
    unsigned int indices[] = {
        // note that we start from 0!
        0, 1, 2 // first Triangle
    };

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);

    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // note that this is allowed, the call to glVertexAttribPointer registered VBO as the vertex attribute's bound vertex buffer object so afterwards we can safely unbind
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // remember: do NOT unbind the EBO while a VAO is active as the bound element buffer object IS stored in the VAO; keep the EBO bound.
    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    // You can unbind the VAO afterwards so other VAO calls won't accidentally modify this VAO, but this rarely happens. Modifying other
    // VAOs requires a call to glBindVertexArray anyways so we generally don't unbind VAOs (nor VBOs) when it's not directly necessary.
    glBindVertexArray(0);
}