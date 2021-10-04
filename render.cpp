#include "render.h"
#include "cameras/ArcBallCamera.h"
#include "RenderCallback.h"
#include "utils/ObjUtil.h"

#include "geometries/Primitives.h"
void cRender::InitGLFormat()
{
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);

    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void *)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // normal attribute
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(float), (void *)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);
}
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
    InitResource();
    InitPtsGL();
    InitCam();
    InitBallGL();
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
    normal_shader->use();

    // ourShader->setMat4("model", glm::mat4(1.0f));
    normal_shader->setMat4("ubo.model", glm::mat4(1.0f));
    tMatrix4f eigen_view = this->mCam->ViewMatrix();
    // eigen_view.setIdentity();
    tMatrix4f eigen_proj = this->mCam->ProjMatrix(mWidth, mHeight, false);
    // eigen_proj.setIdentity();
    glm::mat4 glm_view = E2GLM(eigen_view);
    glm::mat4 glm_proj = E2GLM(eigen_proj);
    normal_shader->setMat4("ubo.view", glm_view);
    normal_shader->setMat4("ubo.proj", glm_proj);
    // std::cout << "cur proj = \n" << eigen_proj << std::endl;
    mCam->MouseRotate();
    glBindVertexArray(triangle_VAO); // seeing as we only have a single VAO there's no need to bind it every time, but we'll do so to keep things a bit more organized

    glDrawArrays(GL_TRIANGLES, 0, 3);
    //glDrawArrays(GL_TRIANGLES, 0, 6);

    glBindVertexArray(axes_VAO);
    glDrawArrays(GL_LINES, 0, 6);

    glBindVertexArray(pts_VAO);
    // UpdatePts();
    // glDrawArrays(GL_LINES, 0, 100);
    glDrawArrays(GL_POINTS, 0, mNumOfPts);

    glBindVertexArray(ball_VAO);
    // glDrawArrays(GL_TRIANGLES, 0, mNumOfVertexBall * 6);
    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ball_EBO);
    glDrawElements(GL_TRIANGLES, 1500, GL_UNSIGNED_INT, 0);
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
    normal_shader = new Shader("assets/shader.vs", "assets/shader.fs");

    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------
    float tri_vertices[] = {
        // positions         // colors          // normals
        0.5f, -0.5f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,  // bottom right
        -0.5f, -0.5f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, // bottom left
        0.0f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f    // top
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
    InitGLFormat();
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
        0.0f,
        0.0f,
        0.0f,
        100.0f,
        0.0f,
        0.0f,
        1.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,

        100.0f,
        0.0f,
        0.0f,
        1.0f,
        0.0f,
        0.0f,
        0.0f,
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
        0.0f,
        0.0f,

        0.0f,
        100.0f,
        0.0f,
        0.0f,
        1.0f,
        0.0f,
        0.0f,
        0.0f,
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
        0.0f,
        0.0f,
        0.0f,
        100.0f,
        0.0f,
        0.0f,
        1.0f,
        0.0f,
        0.0f,
        0.0f,

    };

    // 1. generate VAO and VBO
    glGenVertexArrays(1, &axes_VAO);
    glGenBuffers(1, &axes_VBO);

    glBindVertexArray(axes_VAO);
    glBindBuffer(GL_ARRAY_BUFFER, axes_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(axes_vertices), axes_vertices, GL_STATIC_DRAW);

    // pos attribute
    InitGLFormat();
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
    mPtVec = tVectorXf::Ones(num_of_pts * 3 * 3) * 0.1;
    for (int i = 0; i < num_of_pts; i++)
    {
        const tVector3f &cur_pt = point_coords[i];
        mPtVec[9 * i + 0] = cur_pt[0];
        mPtVec[9 * i + 1] = cur_pt[1];
        mPtVec[9 * i + 2] = cur_pt[2];
        mPtVec[9 * i + 3] = 1.0f;
        mPtVec[9 * i + 4] = 1.0f;
        mPtVec[9 * i + 5] = 1.0f;
        mPtVec[9 * i + 6] = 0.0f;
        mPtVec[9 * i + 7] = 0.0f;
        mPtVec[9 * i + 8] = 0.0f;
    }
    // std::cout << mPtVec.transpose() << std::endl;
    glGenVertexArrays(1, &pts_VAO);
    glGenBuffers(1, &pts_VBO);

    glBindVertexArray(pts_VAO);
    glBindBuffer(GL_ARRAY_BUFFER, pts_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * mPtVec.size(), mPtVec.data(), GL_STATIC_DRAW);

    // pos attribute
    InitGLFormat();
}

#include "png2pointcloud.h"
#include "utils/OpenCVUtil.h"

void cRender::InitResource()
{
    tMatrixXf depth_png = cOpencvUtil::LoadGrayscalePngEigen("depth0.png");
    depth_png /= 255;
    // int height = 100;
    // int width = height;
    // depth_png
    // depth_map.setOnes();
    mPng2PointCloud->Resource(depth_png, 59, mNumOfPts, point_coords);
}

void UpdateVertexNormalFromTriangleNormal(
    std::vector<tVertex *> &v_array,
    std::vector<tTriangle *> &t_array)
{
    for (auto &tri : t_array)
    {
        const tVector &v0 = v_array[tri->mId0]->mPos;
        const tVector &v1 = v_array[tri->mId1]->mPos;
        const tVector &v2 = v_array[tri->mId2]->mPos;
        tri->mNormal = (v1 - v0).cross3(v2 - v1).normalized();
        // std::cout << tri->mNormal.transpose() << std::endl;
    }

    // 1. clear all vertex normal
    // cTimeUtil::Begin("update_v_normal");
    for (auto &x : v_array)
        x->mNormal.setZero();
    // 2. iter each edge
    for (auto &x : t_array)
    {
        v_array[x->mId0]->mNormal += x->mNormal;
        v_array[x->mId1]->mNormal += x->mNormal;
        v_array[x->mId2]->mNormal += x->mNormal;
    }

    // 3. averge each vertex
    for (int i = 0; i < v_array.size(); i++)
    {
        auto &v = v_array[i];
        v->mNormal.normalize();
    }
    // cTimeUtil::End("update_v_normal");
}

void cRender::InitBallGL()
{
    std::string ball_path = "./assets/ball.obj";
    cObjUtil::tParams params;
    params.mPath = ball_path;

    std::vector<tVertex *> v_array;
    std::vector<tEdge *> e_array;
    std::vector<tTriangle *> t_array;
    cObjUtil::LoadObj(params, v_array, e_array, t_array);
    UpdateVertexNormalFromTriangleNormal(v_array, t_array);
    std::cout
        << "load ball from " << ball_path << " succ, num of v " << v_array.size() << " num of e " << e_array.size() << " num of t " << t_array.size() << std::endl;

    tVectorXf render_buf = tVectorXf::Zero(9 * v_array.size());

    for (int i = 0; i < v_array.size(); i++)
    {
        const auto &v = v_array[i];
        render_buf[9 * i + 0] = v->mPos[0] * 1e-1;
        render_buf[9 * i + 1] = v->mPos[1] * 1e-1;
        render_buf[9 * i + 2] = v->mPos[2] * 1e-1;
        render_buf[9 * i + 3] = 1;
        render_buf[9 * i + 4] = 0;
        render_buf[9 * i + 5] = 0;
        render_buf[9 * i + 6] = v->mNormal[0];
        render_buf[9 * i + 7] = v->mNormal[1];
        render_buf[9 * i + 8] = v->mNormal[2];
        std::cout << "normal = " << v->mNormal.transpose() << std::endl;
    }

    std::vector<unsigned int> indices(0);
    for (auto &t : t_array)
    {
        indices.push_back(t->mId0);
        indices.push_back(t->mId1);
        indices.push_back(t->mId2);
    }

    glGenVertexArrays(1, &ball_VAO);
    glGenBuffers(1, &ball_VBO);
    glGenBuffers(1, &ball_EBO);

    glBindVertexArray(ball_VAO);

    // VBO
    glBindBuffer(GL_ARRAY_BUFFER, ball_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * render_buf.size(), render_buf.data(), GL_STATIC_DRAW);

    // EBO
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ball_EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * indices.size(), indices.data(), GL_STATIC_DRAW);

    mBallNumIndices = indices.size();
    // pos attribute
    InitGLFormat();
}