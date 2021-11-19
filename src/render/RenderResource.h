#pragma once
#include "utils/DefUtil.h"
#include "utils/GLUtil.h"
#include "shader.h"
#include "utils/MathUtil.h"
#include "utils/JsonUtil.h"

struct tRenderObj
{
    unsigned int mVAO, mVBO, mEBO;
    unsigned int mNumIndices;
};

struct tRenderResourceBase
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    tRenderResourceBase(const Json::Value &conf);
    virtual void ApplyCameraPose(const tVector3f &cam_pos, const tVector3f &cam_focus, const tVector3f &cam_up);
    std::string mName;
    tEigenArr<tVector3f> mPointCloudArray;
    int mNumOfPoint;
    tVector3f mColor;

    bool mEnableWindow;    // enable window
    tVector2i mRawImgSize; // raw image size
    tVector2i mWindowSt;   //  window st
    tVector2i mWindowSize; // window size
};
struct tRenderResourceSingleImage : public tRenderResourceBase
{
    tRenderResourceSingleImage(const Json::Value &conf);
    std::string mPngPath;
};

struct tRenderResourceMesh4View : public tRenderResourceBase
{
    tRenderResourceMesh4View(const Json::Value &conf);
    virtual void ApplyCameraPose(const tVector3f &cam_pos, const tVector3f &cam_focus, const tVector3f &cam_up);
    std::vector<std::string> mPngPathList;
};

// SIM_DECLARE_PTR(tRenderResourceSingleImage);
SIM_DECLARE_PTR(tRenderResourceBase);
