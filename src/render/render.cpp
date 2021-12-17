#include "render.h"
#include "utils/FileUtil.h"
#include "cameras/ArcBallCamera.h"
#include "RenderCallback.h"
#include "utils/ObjUtil.h"
#include "geometries/Primitives.h"
double fov = 59;
void cRender::InitGLFormat()
{
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);

    // color or normal attribute attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
}
cRender::cRender()
{
    // mNumOfPts = 1;
    // mNeedToRedrawPointCloud = true;
    mLeftButtonPress = false;
}
cRender::~cRender()
{
    // glDeleteVertexArrays(1, &triangle_VAO);
    // glDeleteBuffers(1, &triangle_VBO);
    // glDeleteBuffers(1, &EBO);
    glDeleteProgram(shaderProgram);
}

void cRender::InitCam(const Json::Value &conf)
{

    tVector3f pos = cJsonUtil::ReadVectorJson(conf["gui_cam_pos"]).segment(0, 3).cast<float>();
    tVector3f center = cJsonUtil::ReadVectorJson(conf["gui_cam_focus"]).segment(0, 3).cast<float>();
    tVector3f up = cJsonUtil::ReadVectorJson(conf["gui_cam_up"]).segment(0, 3).cast<float>();
    float near_plane_dist = 1e-3;
    float far_plane_dist = 1e2;

    mCam = std::make_shared<cArcBallCamera>(pos,
                                            center,
                                            up,
                                            60,
                                            near_plane_dist,
                                            far_plane_dist);
}
#include "utils/JsonUtil.h"
#include "src/sim_kinect/SimKinect.h"
cSimKinectPtr sim_kinect = nullptr;
void cRender::Init(std::string conf_path)
{
    Json::Value root;
    SIM_ASSERT(cJsonUtil::LoadJson(conf_path, root) == true);
    Json::Value resource_lst = cJsonUtil::ParseAsValue("resource_lst", root);
    sim_kinect = std::make_shared<cSimKinect>();
    sim_kinect->Init();
    {
        // tMatrixXf img = cOpencvUtil::LoadGrayscalePngEigen("assets\\cufanggezi_ps160\\0.png");
        // img = sim_kinect->AddContinousNoise(img);
        // exit(1);
    }
    // exit(1);
    // Json::Value mesh_4view_lst = cJsonUtil::ParseAsValue("mesh_4view_lst", root);

    auto format_ptr = std::make_shared<tImageFormat>(cJsonUtil::ParseAsValue("image_format_info", root));
    mEnableRenderResource.resize(resource_lst.size(), true);
    mEnableTransformAdjust.resize(resource_lst.size(), true);
    // std::cout << "[debug] resource single list num = " << resource_lst.size() << std::endl;
    // std::cout << "[debug] mesh_4view list num = " << mesh_4view_lst.size() << std::endl;

    for (auto &x : resource_lst)
    {
        std::string type_str = cJsonUtil::ParseAsString("type", x);
        auto cur_type = BuildRenderResourceTypeFromStr(type_str);
        tRenderResourceBasePtr new_res = nullptr;
        switch (cur_type)
        {
        case eRenderResourceType::SINGLE_VIEW_IMAGE_RESOURCE_TYPE:
        {
            new_res = std::make_shared<tRenderResourceSingleImage>(x, format_ptr);
        }
        break;
        case eRenderResourceType::MULTIVIEW_IMAGE_RESOURCE_TYPE:
        {
            new_res = std::make_shared<tRenderResourceMesh4View>(x, format_ptr);
        }
        break;
        case eRenderResourceType::OBJ_RESOURCE_TYPE:
        {
            new_res = std::make_shared<tRenderResourceObj>(x);
            // std::cout << "need to support for obj resource\n";
            // exit(1);
        }
        break;
        default:
        {
            SIM_ERROR("we have no option for resource type {}", cur_type);
            exit(1);
        }
        }
        mRenderResources.push_back(new_res);
    }

    mPngCamPos = cJsonUtil::ReadVectorJson(cJsonUtil::ParseAsValue("camera_pos", root)).segment(0, 3).cast<float>();
    mPngCamFocus = cJsonUtil::ReadVectorJson(cJsonUtil::ParseAsValue("camera_focus", root)).segment(0, 3).cast<float>();
    mPngCamUp = cJsonUtil::ReadVectorJson(cJsonUtil::ParseAsValue("camera_up", root)).segment(0, 3).cast<float>();
    mPngCamUp = mPngCamUp.normalized();
    // for png resource, coords to world coords

    for (auto &res : mRenderResources)
    {
        auto img_res = std::dynamic_pointer_cast<tRenderResourceImageBase>(res);
        if (img_res != nullptr)
        {
            img_res->ApplyCameraPose(mPngCamPos, mPngCamFocus, mPngCamUp);
        }
        res->InitPointCloudArray();
    }

    // for

    InitGL();
    InitAxesGL();
    InitPtsGL();
    InitCam(root["GUI_configuration"]);
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
    SetCamInShader(normal_shader);
    SetCamInShader(ball_shader);

    normal_shader->use();
    // std::cout << "cur proj = \n" << eigen_proj << std::endl;
    // mCam->MouseRotate();
    // glBindVertexArray(triangle_VAO); // seeing as we only have a single VAO there's no need to bind it every time, but we'll do so to keep things a bit more organized

    // glDrawArrays(GL_TRIANGLES, 0, 3);
    //glDrawArrays(GL_TRIANGLES, 0, 6);

    glBindVertexArray(axes_VAO);
    glDrawArrays(GL_LINES, 0, 6);

    // glBindVertexArray(pts_VAO);
    // UpdatePts();
    // glDrawArrays(GL_LINES, 0, 100);
    // glDrawArrays(GL_POINTS, 0, mNumOfPts);

    // glDrawArrays(GL_TRIANGLES, 0, mNumOfVertexBall * 6);
    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ball_EBO);

    // if (mNeedToRedrawPointCloud)
    {
        ball_shader->use();
        glBindVertexArray(mBallObj.mVAO);
        tMatrix4f model = tMatrix4f::Identity();
        for (int i = 0; i < mRenderResources.size(); i++)
        {
            if (mEnableRenderResource[i] == false)
                continue;
            auto &res = mRenderResources[i];

            if (std::dynamic_pointer_cast<tRenderResourceMesh4View>(res) == nullptr)
            {
                // std::cout << "single resource\n";
                ball_shader->setVec3("ball_color", glm::vec3(res->mColor[0], res->mColor[1], res->mColor[2]));
                for (auto &pt : res->mPointCloudArray)
                {
                    model.block(0, 3, 3, 1) = pt;

                    ball_shader->setMat4("ubo.model", E2GLM(model));
                    glDrawElements(GL_TRIANGLES, mBallObj.mNumIndices, GL_UNSIGNED_INT, 0);
                }
            }
            else
            {
                // 4 views
                // std::cout << "multiview resource\n";
                for (int i = 0; i < 2; i++)
                {
                    tVector3f cur_color = res->mColor * (i + 1) * 0.25;
                    ball_shader->setVec3("ball_color", glm::vec3(cur_color[0], cur_color[1], cur_color[2]));
                    int num_of_pts_gap = res->mPointCloudArray.size() / 4;
                    for (int _idx = num_of_pts_gap * i; _idx < num_of_pts_gap * (i + 1); _idx++)
                    {
                        model.block(0, 3, 3, 1) = res->mPointCloudArray[_idx];

                        ball_shader->setMat4("ubo.model", E2GLM(model));
                        glDrawElements(GL_TRIANGLES, mBallObj.mNumIndices, GL_UNSIGNED_INT, 0);
                    }
                }
            }
        }
        // mNeedToRedrawPointCloud = true;
    }
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
    normal_shader = new Shader("assets/normal_shader.vs", "assets/shader.fs");
    ball_shader = new Shader("assets/ball_shader.vs", "assets/shader.fs");

    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------
    // float tri_vertices[] = {
    //     // positions         // colors          // normals
    //     0.5f, -0.5f, 0.0f, 1.0f, 0.0f, 0.0f,  // bottom right
    //     -0.5f, -0.5f, 0.0f, 0.0f, 1.0f, 0.0f, // bottom left
    //     0.0f, 0.5f, 0.0f, 0.0f, 0.0f, 1.0f    // top
    // };
    // unsigned int indices[] = {
    //     // note that we start from 0!
    //     0, 1, 2 // first Triangle
    // };

    // // 1. generate VAO and VBO
    // glGenVertexArrays(1, &triangle_VAO);
    // glGenBuffers(1, &triangle_VBO);

    // glBindVertexArray(triangle_VAO);
    // glBindBuffer(GL_ARRAY_BUFFER, triangle_VBO);
    // glBufferData(GL_ARRAY_BUFFER, sizeof(tri_vertices), tri_vertices, GL_STATIC_DRAW);

    // pos attribute
    InitGLFormat();
}

void cRender::InitAxesGL()
{
    float axes_vertices[] = {
        // positions         // colors

        0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,   // X axis start
        100.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, // X axis end
        0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,   // Y axis start
        0.0f, 100.0f, 0.0f, 0.0f, 1.0f, 0.0f, // Y axis end
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,   // Z axis sart
        0.0f, 0.0f, 100.0f, 0.0f, 0.0f, 1.0f  // Z axis end
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
    // printf("mouse move to %f %f, left button %d\n", xpos, ypos, mLeftButtonPress);
    if (mLeftButtonPress)
    {
        // printf("mouse move to %d %d\n", xpos, ypos);
        mCam->MouseMove(xpos, ypos);
    }
}
void cRender::MouseButtonCallback(int button, int action, int mods)
{
    // std::cout << "button = " << button << std::endl;
    // std::cout << "action = " << action << std::endl;
    if (button == GLFW_MOUSE_BUTTON_1)
    {
        if (action == GLFW_RELEASE)
        {
            mLeftButtonPress = false;
            mCam->ResetFlag();
        }

        else if (action == GLFW_PRESS)
        {
            mLeftButtonPress = true;
            // if (mCamera->IsFirstMouse() == true)
            // {
            //     mCamera->MouseMove()
            // }
        }
    }
}
void cRender::KeyCallback(int key, int scancode, int action, int mods)
{
}
void cRender::ResizeCallback(int w, int h)
{
}
void cRender::ScrollCallback(double xoff, double yoff)
{
    if (yoff > 0)
        mCam->MoveForward();
    else
        mCam->MoveBackward();
}

void cRender::InitPtsGL()
{
}

void cRender::AddResource(const Json::Value &conf)
{

    // std::cout << "load resource png from " << png_path << std::endl;
    // cPng2PointCloud::Resource(png_path, cam_vfov_deg, mRenderScene.mNumOfPoint, mRenderScene.mPointCloudArray);
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
        render_buf[6 * i + 0] = v->mPos[0] * 1e-3;
        render_buf[6 * i + 1] = v->mPos[1] * 1e-3;
        render_buf[6 * i + 2] = v->mPos[2] * 1e-3;
        render_buf[6 * i + 3] = v->mNormal[0];
        render_buf[6 * i + 4] = v->mNormal[1];
        render_buf[6 * i + 5] = v->mNormal[2];
        // std::cout << "normal = " << v->mNormal.transpose() << std::endl;
    }

    std::vector<unsigned int> indices(0);
    for (auto &t : t_array)
    {
        indices.push_back(t->mId0);
        indices.push_back(t->mId1);
        indices.push_back(t->mId2);
    }

    glGenVertexArrays(1, &mBallObj.mVAO);
    glGenBuffers(1, &mBallObj.mVBO);
    glGenBuffers(1, &mBallObj.mEBO);

    glBindVertexArray(mBallObj.mVAO);

    // VBO
    glBindBuffer(GL_ARRAY_BUFFER, mBallObj.mVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * render_buf.size(), render_buf.data(), GL_STATIC_DRAW);

    // EBO
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mBallObj.mEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * indices.size(), indices.data(), GL_STATIC_DRAW);

    mBallObj.mNumIndices = indices.size();
    // pos attribute
    InitGLFormat();
}

void cRender::SetCamInShader(Shader *this_shader) const
{
    this_shader->use();

    // ourShader->setMat4("model", glm::mat4(1.0f));
    this_shader->setMat4("ubo.model", glm::mat4(1.0f));
    tMatrix4f eigen_view = this->mCam->ViewMatrix();
    // eigen_view.setIdentity();
    tMatrix4f eigen_proj = this->mCam->ProjMatrix(mWidth, mHeight, false);
    // eigen_proj.setIdentity();
    glm::mat4 glm_view = E2GLM(eigen_view);
    glm::mat4 glm_proj = E2GLM(eigen_proj);
    this_shader->setMat4("ubo.view", glm_view);
    this_shader->setMat4("ubo.proj", glm_proj);
}
