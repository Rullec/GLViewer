#pragma once
#include "RenderResource.h"
#include <deque>
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
    virtual void AddResource(const Json::Value &conf);
    virtual void Update();
    virtual void MouseMoveCallback(double xpos, double ypos);
    virtual void MouseButtonCallback(int but, int action, int mods);
    virtual void KeyCallback(int key, int scancode, int action, int mods);
    virtual void ResizeCallback(int w, int h);
    virtual void ScrollCallback(double xoff, double yoff);

protected:
    GLFWwindow *mWindow;
    int mWidth = 1200, mHeight = 800;
    int mStartX = 100, mStartY = 100;
    // int mNumOfPts;
    std::string mWindowName = "gl viewer";
    // unsigned int triangle_VBO, triangle_VAO;
    unsigned int axes_VBO, axes_VAO;
    tRenderObj mBallObj;
    // unsigned int pts_VBO, pts_VAO;
    // unsigned int ball_VAO, ball_VBO, ball_EBO;
    unsigned int shaderProgram;
    // tVectorXf mPtVec;
    Shader *normal_shader, *ball_shader;
    CameraBasePtr mCam;
    cPng2PointCloudPtr mPng2PointCloud;
    std::deque<bool> mEnableRenderResource; // enable render resource or not?
    std::deque<bool> mEnableTransformAdjust;
    std::vector<tRenderResourceBasePtr> mRenderResources;
    tVector3f mPngCamPos, mPngCamFocus, mPngCamUp;
    // bool mEnableTransformDepthImageToWorldCoords; // enable transofrm depth image to wrold coords or not
    // tEigenArr<tVector3f> point_coords;
    // bool mNeedToRedrawPointCloud;
    // unsigned int mBallNumIndices;
    bool mLeftButtonPress;
    virtual void InitCam(const Json::Value &conf);
    virtual void InitGL();
    virtual void InitAxesGL();
    virtual void InitPtsGL();
    virtual void InitBallGL();
    virtual void InitGLFormat();
    virtual void SetCamInShader(Shader *this_shader) const;
};

SIM_DECLARE_PTR(cRender);