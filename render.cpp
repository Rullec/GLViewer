#include "render.h"
#include "cameras/ArcBallCamera.h"
#include "RenderCallback.h"

cRender::cRender()
{
    mNumOfPts = 1;
}
cRender::~cRender()
{

    glDeleteVertexArrays(1, &triangle_VAO);
    glDeleteBuffers(1, &triangle_VBO);
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
    InitAxesGL();
    InitPtsGL();
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
    glBindVertexArray(triangle_VAO); // seeing as we only have a single VAO there's no need to bind it every time, but we'll do so to keep things a bit more organized

    glDrawArrays(GL_TRIANGLES, 0, 3);
    //glDrawArrays(GL_TRIANGLES, 0, 6);

    glBindVertexArray(axes_VAO);
    glDrawArrays(GL_LINES, 0, 6);

    glBindVertexArray(pts_VAO);
    UpdatePts();
    // glDrawArrays(GL_LINES, 0, 100);
    glDrawArrays(GL_POINTS, 0, mNumOfPts);
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
    glfwSetCursorPosCallback(mWindow, MouseMoveEventCallback);
    glfwSetInputMode(mWindow, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    glfwSetFramebufferSizeCallback(mWindow, ResizeEventCallback);
    glfwSetScrollCallback(mWindow, ScrollEventCallback);
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
    float tri_vertices[] = {
        // positions         // colors
        0.5f, -0.5f, 0.0f, 1.0f, 0.0f, 0.0f,  // bottom right
        -0.5f, -0.5f, 0.0f, 0.0f, 1.0f, 0.0f, // bottom left
        0.0f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f    // top
    };
    unsigned int indices[] = {
        // note that we start from 0!
        0, 1, 2 // first Triangle
    };

    // 1. generate VAO and VBO
    glGenVertexArrays(1, &triangle_VAO);
    glGenBuffers(1, &triangle_VBO);

    glBindVertexArray(triangle_VAO);
    glBindBuffer(GL_ARRAY_BUFFER, triangle_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(tri_vertices), tri_vertices, GL_STATIC_DRAW);

    // pos attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);

    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
}

void cRender::InitAxesGL()
{
    float axes_vertices[] = {
        // positions         // colors
        // X axis
        0.0f,
        0.0f,
        0.0f,
        1.0f,
        0.0f,
        0.0f,

        100.0f,
        0.0f,
        0.0f,
        1.0f,
        0.0f,
        0.0f,

        // Y axis
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        1.0f,
        0.0f,

        0.0f,
        100.0f,
        0.0f,
        0.0f,
        1.0f,
        0.0f,

        // Z axis
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        1.0f,

        0.0f,
        0.0f,
        100.0f,
        0.0f,
        0.0f,
        1.0f,

    };

    // 1. generate VAO and VBO
    glGenVertexArrays(1, &axes_VAO);
    glGenBuffers(1, &axes_VBO);

    glBindVertexArray(axes_VAO);
    glBindBuffer(GL_ARRAY_BUFFER, axes_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(axes_vertices), axes_vertices, GL_STATIC_DRAW);

    // pos attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);

    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
}

void cRender::MouseMoveCallback(double xpos, double ypos)
{
    std::cout << "render move " << xpos << " " << ypos << std::endl;
}
void cRender::MouseButtonCallback(int but, int action, int mods)
{
}
void cRender::KeyCallback(int key, int scancode, int action, int mods)
{
}
void cRender::ResizeCallback(int w, int h)
{
}
void cRender::ScrollCallback(double xoff, double yoff)
{
}

void cRender::InitPtsGL()
{
    int num_of_pts = mNumOfPts;
    mPtVec = tVectorXf::Ones(num_of_pts * 3 * 2) * 0.1;
    for (int i = 0; i < num_of_pts; i++)
    {
        // mPtVec[6 * i + 0] = float(i) * 0.01;
        // mPtVec[6 * i + 1] = float(i) * 0.01;
        // mPtVec[6 * i + 2] = float(i) * 0.01;
        mPtVec[6 * i + 3] = 1.0f;
        mPtVec[6 * i + 4] = 1.0f;
        mPtVec[6 * i + 5] = 1.0f;
    }
    std::cout << mPtVec.transpose() << std::endl;
    glGenVertexArrays(1, &pts_VAO);
    glGenBuffers(1, &pts_VBO);

    glBindVertexArray(pts_VAO);
    glBindBuffer(GL_ARRAY_BUFFER, pts_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * mPtVec.size(), mPtVec.data(), GL_STATIC_DRAW);

    // pos attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);

    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
}

void cRender::UpdatePts()
{
    mNumOfPts += 1;
    int num_of_pts = mNumOfPts;
    tVectorXf new_pts_vertices = tVectorXf::Random(num_of_pts * 3 * 2);
    new_pts_vertices.segment(0, 6 * (mNumOfPts - 1)) = mPtVec;
    for (int i = 0; i < num_of_pts; i++)
    {
        // pts_vertices[6 * i + 0] = float(i) * 0.01;
        // pts_vertices[6 * i + 1] = float(i) * 0.01;
        // pts_vertices[6 * i + 2] = float(i) * 0.01;
        new_pts_vertices[6 * i + 3] = 1.0f;
        new_pts_vertices[6 * i + 4] = 1.0f;
        new_pts_vertices[6 * i + 5] = 1.0f;
    }

    glBindBuffer(GL_ARRAY_BUFFER, pts_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * new_pts_vertices.size(), new_pts_vertices.data(), GL_STATIC_DRAW);
    mPtVec = new_pts_vertices;
}