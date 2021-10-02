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

protected:
    GLFWwindow *mWindow;
    int mWidth = 300, mHeight = 300;
    int mStartX = 300, mStartY = 300;
    std::string mWindowName = "gl viewer";
    unsigned int VBO, VAO;
    unsigned int shaderProgram;
    Shader *ourShader;
    CameraBasePtr mCam;


    virtual void InitCam();
    virtual void InitGL();
};

SIM_DECLARE_PTR(cRender);