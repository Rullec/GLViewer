#pragma once
#include "utils/DefUtil.h"
#include "utils/GLUtil.h"
#include "shader.h"
#include "utils/MathUtil.h"
SIM_DECLARE_CLASS_AND_PTR(CameraBase);
SIM_DECLARE_CLASS_AND_PTR(cPng2PointCloud);

class cRender
{
public:
    explicit cRender();
    virtual void Init(std::string conf_path);
    virtual ~cRender();
    virtual GLFWwindow *GetWindow()
    {
        return mWindow;
    }
    void InitResource(std::string png_path);
    virtual void Update();
    void MouseMoveCallback(double xpos, double ypos);
    void MouseButtonCallback(int but, int action, int mods);
    void KeyCallback(int key, int scancode, int action, int mods);
    void ResizeCallback(int w, int h);
    void ScrollCallback(double xoff, double yoff);

protected:
    GLFWwindow *mWindow;
    int mWidth = 800, mHeight = 600;
    int mStartX = 100, mStartY = 100;
    int mNumOfPts;
    std::string mWindowName = "gl viewer";
    unsigned int triangle_VBO, triangle_VAO;
    unsigned int axes_VBO, axes_VAO;
    // unsigned int pts_VBO, pts_VAO;
    unsigned int ball_VAO, ball_VBO, ball_EBO;
    unsigned int shaderProgram;
    tVectorXf mPtVec;
    Shader *normal_shader;
    CameraBasePtr mCam;
    cPng2PointCloudPtr mPng2PointCloud;
    tEigenArr<tVector3f> point_coords;
    bool mNeedToRedrawPointCloud;
    unsigned int mBallNumIndices;
    bool mLeftButtonPress;
    virtual void InitCam();
    virtual void InitGL();
    virtual void InitAxesGL();
    virtual void InitPtsGL();
    virtual void InitBallGL();
    void InitGLFormat();
};

SIM_DECLARE_PTR(cRender);