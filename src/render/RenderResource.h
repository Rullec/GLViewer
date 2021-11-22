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

/**
 * \brief           the image can be represented by two types: png and txt. txt is more accurate
*/
enum eRenderResourceImageType
{
    PNG_RENDER_RESOURCE_IMAGE_TYPE = 0,
    TXT_RENDER_RESOURCE_IMAGE_TYPE,
    NUM_RENDER_RESOURCE_IMAGE_TYPE
};
enum eRenderResourceType
{
    SINGLE_VIEW_IMAGE_RESOURCE_TYPE = 0,
    MULTIVIEW_IMAGE_RESOURCE_TYPE = 1,
    OBJ_RESOURCE_TYPE,
    NUM_RESOURCE_TYPE
};

eRenderResourceType BuildRenderResourceTypeFromStr(std::string);
std::string BuildStrFromRenderResourceType(eRenderResourceType);
struct tImageFormat
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    tImageFormat(const Json::Value &conf);

    bool mEnableWindow; // enable window

    tVector2i mRawImgSize; // raw image size
    tVector2i mWindowSt;   //  window st
    tVector2i mWindowSize; // window size
};
SIM_DECLARE_PTR(tImageFormat);

struct tRenderResourceBase
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    tRenderResourceBase(const Json::Value &conf);
    std::string mName;
    tEigenArr<tVector3f> mPointCloudArray;

    int mNumOfPoint;
    tVector3f mColor;
    tMatrix mTransform;
    virtual void CalcAABB(tVector3f &aabb_min, tVector3f &aabb_max);
    virtual void ApplyTransform();
    virtual void SetPos(const tVector3f &new_pos);
    virtual tVector3f GetPos() const;
    virtual void SetRotAxisAngle(const tVector3f & new_aa);
    virtual tVector3f GetRotAxisAngle() const;
    virtual void InitPointCloudArray();

protected:
    tEigenArr<tVector3f> mInitPointCloudArray;
};
SIM_DECLARE_PTR(tRenderResourceBase);

struct tRenderResourceObj : public tRenderResourceBase
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    tRenderResourceObj(const Json::Value &conf);
    std::string mObjPath;
    float mScale, mRandom;
};

struct tRenderResourceImageBase : tRenderResourceBase
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    tRenderResourceImageBase(const Json::Value &conf, tImageFormatPtr ptr);
    static eRenderResourceImageType GetImageTypeFromPath(std::string name);
    virtual void ApplyCameraPose(const tVector3f &cam_pos, const tVector3f &cam_focus, const tVector3f &cam_up);

    tImageFormatPtr mFormat;
};

struct tRenderResourceSingleImage : public tRenderResourceImageBase
{
    tRenderResourceSingleImage(const Json::Value &conf, tImageFormatPtr ptr);
    std::string mResourcePath;
    eRenderResourceImageType mRenderResourceType;
};

struct tRenderResourceMesh4View : public tRenderResourceImageBase
{
    tRenderResourceMesh4View(const Json::Value &conf, tImageFormatPtr ptr);
    virtual void ApplyCameraPose(const tVector3f &cam_pos, const tVector3f &cam_focus, const tVector3f &cam_up);
    std::vector<std::string> mResourcePathList;
    std::vector<eRenderResourceImageType> mResourceTypeList;
};

// SIM_DECLARE_PTR(tRenderResourceSingleImage);
SIM_DECLARE_PTR(tRenderResourceImageBase);
SIM_DECLARE_PTR(tRenderResourceMesh4View);
