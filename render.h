#pragma once
#include "utils/DefUtil.h"
#include "utils/GLUtil.h"
#include "shader.h"

SIM_DECLARE_CLASS_AND_PTR(CameraBase);

class cRender
{
public:
    explicit cRender();
    virtual void Init();
    virtual ~cRender();
    virtual GLFWwindow *GetWindow()
    {
        return mWindow;
    }
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
    std::string mWindowName = "gl viewer";
    unsigned int triangle_VBO, triangle_VAO;
    unsigned int axes_VBO, axes_VAO;
    unsigned int pts_VBO, pts_VAO;
    unsigned int shaderProgram;
    Shader *ourShader;
    CameraBasePtr mCam;

    virtual void InitCam();
    virtual void InitGL();
    virtual void InitAxesGL();
    virtual void InitPtsGL();
};

SIM_DECLARE_PTR(cRender);