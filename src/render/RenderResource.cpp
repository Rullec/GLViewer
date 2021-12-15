#include "RenderResource.h"
#include "restore/png2pointcloud.h"
#include "utils/OpenCVUtil.h"
#include "utils/FileUtil.h"
#include "utils/StringUtil.h"
#include "sim_kinect/SimKinect.h"
extern double fov;
extern cSimKinectPtr sim_kinect;
const std::string gRenderResourceTypStr[eRenderResourceType::NUM_RESOURCE_TYPE] = {
    "single_view_image",
    "multi_view_image",
    "obj"};
tImageFormat::tImageFormat(const Json::Value &conf)
{
    mRawImgSize = cJsonUtil::ReadVectorJson(
                      cJsonUtil::ParseAsValue("raw_img_size", conf))
                      .segment(0, 2)
                      .cast<int>();
    mEnableWindow = cJsonUtil::ParseAsBool("enable_window", conf);
    if (mEnableWindow)
    {
        int window_height_st = cJsonUtil::ParseAsInt("window_height_st", conf);
        int window_width_st = cJsonUtil::ParseAsInt("window_width_st", conf);
        int window_height_size = cJsonUtil::ParseAsInt("window_height_size", conf);
        int window_width_size = cJsonUtil::ParseAsInt("window_width_size", conf);
        mWindowSt[0] = window_height_st;
        mWindowSt[1] = window_width_st;
        mWindowSize[0] = window_height_size;
        mWindowSize[1] = window_width_size;
    }
}

/**
 * \brief            load depth from txt, return in meter unit
*/
tMatrixXf LoadTxtDepthMeter(const std::string &path)
{
    std::string content = cFileUtil::ReadTextFile(path);
    // std::cout << "content = " << content << std::endl;
    std::vector<std::string> lines = cStringUtil::SplitString(content, "\n");
    cStringUtil::RemoveEmptyLine(lines);
    lines.erase(lines.begin());

    std::vector<std::vector<float>> imgs = {};
    for (auto &line : lines)
    {
        auto nums = cStringUtil::SplitString(line, " ");
        cStringUtil::RemoveEmptyLine(nums);
        std::vector<float> floats = {};
        for (auto &x : nums)
        {
            floats.push_back(std::stof(x));
        }
        if (imgs.size() != 0)
        {
            int cur_size = floats.size();
            int pre_size = imgs[imgs.size() - 1].size();
            if (cur_size != pre_size)
            {
                SIM_ERROR("cur size {} != pre size {}, parse depth image failed", cur_size, pre_size);
            }
        }
        imgs.push_back(floats);
    }

    tMatrixXf img = tMatrixXf::Zero(imgs.size(), imgs[0].size());
    for (int i = 0; i < imgs.size(); i++)
    {
        for (int j = 0; j < imgs[0].size(); j++)
        {
            img(i, j) = imgs[i][j];
        }
    }
    std::cout << "parse text file img = " << img.rows() << " " << img.cols() << std::endl;
    return img;
    // for (int i = 0; i < lines.size(); i++)
    // {
    //     auto &line = lines[i];
    //     std::cout << "line " << i << " = " << line << std::endl;
    // }
    return tMatrixXf::Zero(0, 0);
}

#include "utils/ObjUtil.h"
#include "geometries/Primitives.h"
tRenderResourceObj::tRenderResourceObj(const Json::Value &conf)
    : tRenderResourceBase(conf)
{
    mObjPath = cJsonUtil::ParseAsString("obj_path", conf);
    // load the obj, into the vertices
    mScale = cJsonUtil::ParseAsFloat("scale", conf);
    mRandom = cJsonUtil::ParseAsFloat("random_sample", conf);
    {
        cObjUtil::tParams params;
        params.mPath = mObjPath;
        std::vector<tVertex *> v_array;
        std::vector<tEdge *> e_array;
        std::vector<tTriangle *> t_array;
        cObjUtil::LoadObj(params, v_array, e_array, t_array);
        mNumOfPoint = v_array.size();
        mPointCloudArray.clear();
        for (auto &v : v_array)
        {
            if (cMathUtil::RandDouble(0, 1) < mRandom)
            {

                mPointCloudArray.push_back(v->mPos.segment(0, 3).cast<float>() * mScale);
            }
        }
        mNumOfPoint = mPointCloudArray.size();
    }
    std::cout << "[debug] load obj from " << mObjPath << " num of points " << mNumOfPoint << std::endl;
}

tRenderResourceSingleImage::tRenderResourceSingleImage(const Json::Value &conf, tImageFormatPtr ptr) : tRenderResourceImageBase(conf, ptr)
{
    mResourcePath = cJsonUtil::ParseAsString("resource_path", conf);
    // 1. load the png, check the shape
    SIM_ASSERT(cFileUtil::ExistsFile(mResourcePath) == true);
    mRenderResourceType = GetImageTypeFromPath(mResourcePath);
    LoadResource();
}

void tRenderResourceSingleImage::LoadResource()
{
    tMatrixXf img = tMatrixXf::Zero(0, 0);
    switch (mRenderResourceType)
    {
    case eRenderResourceImageType::PNG_RENDER_RESOURCE_IMAGE_TYPE:
    {
        img = cOpencvUtil::LoadGrayscalePngEigen(mResourcePath);
        img /= 255;
    }
    break;
    case eRenderResourceImageType::TXT_RENDER_RESOURCE_IMAGE_TYPE:
    {
        img = LoadTxtDepthMeter(mResourcePath);
        // SIM_ERROR("need to work on txt resource");
        // exit(1);
    }
    break;
    default:
        SIM_ERROR("no policy can work for type {}", this->mRenderResourceType);
        exit(1);
    }

    if (mFormat->mEnableWindow == true)
    {
        // the image must be windowed
        SIM_ASSERT(img.rows() == mFormat->mWindowSize[0]);
        SIM_ASSERT(img.cols() == mFormat->mWindowSize[1]);
    }
    else
    {
        // the image must be full, raw window
        std::cout << "raw img size = " << mFormat->mRawImgSize.transpose() << std::endl;
        SIM_ASSERT(img.rows() == mFormat->mRawImgSize[0]);
        SIM_ASSERT(img.cols() == mFormat->mRawImgSize[1]);
    }

    // 2. convert it to the resource
    if (mFormat->mEnableWindow == false)
    {
        cPng2PointCloud::ResourceWhole(img, fov, mNumOfPoint, mPointCloudArray);
    }
    else
    {
        cPng2PointCloud::ResourceWindow(img, fov,
                                        mFormat->mRawImgSize, mFormat->mWindowSt,
                                        mNumOfPoint,
                                        mPointCloudArray);
    }
    SIM_INFO("load {} points from {}", mNumOfPoint, this->mResourcePath);
}

tRenderResourceMesh4View::tRenderResourceMesh4View(const Json::Value &conf, tImageFormatPtr ptr)
    : tRenderResourceImageBase(conf, ptr)
{
    mResourcePathList.clear();

    // 1. load the png, check the shape
    Json::Value resource_path_lst_json = cJsonUtil::ParseAsValue("resource_path_lst", conf);
    mNumOfPoint = 0;
    mPointCloudArray.clear();

    for (int i = 0; i < resource_path_lst_json.size(); i++)
    {
        std::string cur_str = resource_path_lst_json[i].asString();
        mResourcePathList.push_back(cur_str);
        if (cFileUtil::ExistsFile(cur_str) == false)
        {
            SIM_ERROR("{} doesn't exist", cur_str);
            exit(1);
        }
    }

    LoadResource();
    // SIM_INFO("load {} points from {}, {}, {}, {}", mNumOfPoint, mResourcePathList[0], mResourcePathList[1], mResourcePathList[2], mResourcePathList[3]);
}

void tRenderResourceMesh4View::LoadResource()
{
    mNumOfPoint = 0;
    mPointCloudArray.clear();
    for (int i = 0; i < mResourcePathList.size(); i++)
    {
        auto cur_str = mResourcePathList[i];
        // std::cout << "load str " << cur_str << std::endl;
        tMatrixXf img;
        if (this->mEnableKinectNoise == true)
        {
            // std::cout << "enable kinect noise! load from " << cur_str << std::endl;
            img = sim_kinect->Calculate(cur_str);
            img /= 255;
        }
        else
        {
            img.noalias() = cOpencvUtil::LoadGrayscalePngEigen(cur_str);
            img /= 255;
        }
        if (mFormat->mEnableWindow == true)

        {
            // the image must be windowed
            SIM_ASSERT(img.rows() == mFormat->mWindowSize[0]);
            SIM_ASSERT(img.cols() == mFormat->mWindowSize[1]);
        }
        else
        {
            // the image must be full, raw window
            std::cout << "raw img size = " << mFormat->mRawImgSize.transpose() << std::endl;
            SIM_ASSERT(img.rows() == mFormat->mRawImgSize[0]);
            SIM_ASSERT(img.cols() == mFormat->mRawImgSize[1]);
        }
        int num_of_point_cur = 0;
        // 2. convert it to the resource
        tEigenArr<tVector3f> cur_point_cloud_array(0);
        if (mFormat->mEnableWindow == false)
        {
            cPng2PointCloud::ResourceWhole(img, fov, num_of_point_cur, cur_point_cloud_array);
        }
        else
        {
            cPng2PointCloud::ResourceWindow(img, fov,
                                            mFormat->mRawImgSize, mFormat->mWindowSt,
                                            num_of_point_cur,
                                            cur_point_cloud_array);
        }
        this->mNumOfPoint += num_of_point_cur;
        for (int j = 0; j < cur_point_cloud_array.size(); j++)
        {
            mPointCloudArray.push_back(cur_point_cloud_array[j]);
        }
    }
    InitPointCloudArray();
}
/**
 * \brief           Get type from path
*/
eRenderResourceImageType tRenderResourceImageBase::GetImageTypeFromPath(std::string name)
{
    std::string ext = cFileUtil::GetExtension(name);
    if ("png" == ext)
    {
        return eRenderResourceImageType::PNG_RENDER_RESOURCE_IMAGE_TYPE;
    }
    else if ("txt" == ext)
    {
        return eRenderResourceImageType::TXT_RENDER_RESOURCE_IMAGE_TYPE;
    }
    else
    {
        SIM_ERROR("fail to get render resource type from path {}", name);
        exit(1);
    }
    return eRenderResourceImageType::NUM_RENDER_RESOURCE_IMAGE_TYPE;
}

tRenderResourceBase::tRenderResourceBase(const Json::Value &conf)
{
    mName = cJsonUtil::ParseAsString("name", conf);
    mColor = cJsonUtil::ReadVectorJson(cJsonUtil::ParseAsValue("color", conf)).segment(0, 3).cast<float>();
    mTransform.setIdentity();
}

tRenderResourceImageBase::tRenderResourceImageBase(const Json::Value &conf, tImageFormatPtr ptr) : tRenderResourceBase(conf)
{
    mEnableKinectNoise = false;
    mFormat = ptr;
}

void tRenderResourceImageBase::SetEnableKinectNoise(bool val)
{
    if (mEnableKinectNoise == val)
    {
        return;
    }
    else
    {
        mEnableKinectNoise = val;
        std::cout << "begin to relaod resource\n";
        LoadResource();
    }
}
/**
 * \brief           Apply camera pos, to change the point cloud array
*/
void tRenderResourceImageBase::ApplyCameraPose(const tVector3f &cam_pos, const tVector3f &cam_focus, const tVector3f &cam_up_)
{
    // up : Y axis
    // (focus - pos) : -Z axis
    // X: Y.cross3(Z)
    tVector3f Z = (cam_pos - cam_focus).normalized();
    // std::cout << "Z = " << Z.transpose() << std::endl;
    tVector3f cam_up = cam_up_.dot(Z) * (-Z) + cam_up_;
    // std::cout << "new up = " << mPngCamUp.transpose() << std::endl;
    cam_up.normalize();
    tVector3f Y = cam_up;
    tVector3f X = (Y.cross(Z)).normalized();
    tMatrix3f R = tMatrix3f::Zero();
    R.col(0) = X;
    R.col(1) = Y;
    R.col(2) = Z;
    // std::cout << "[inner] R = \n" << R << std::endl;
    // std::cout << "[inner] cam_pos = " << cam_pos.transpose() << std::endl;
    for (auto &x : mPointCloudArray)
    {
        x = R * x + cam_pos;
    }
}
void tRenderResourceBase::InitPointCloudArray()
{
    mInitPointCloudArray = mPointCloudArray;
}

void tRenderResourceBase::ApplyTransform()
{
    for (int i = 0; i < mNumOfPoint; i++)
    {
        mPointCloudArray[i] = (mTransform * cMathUtil::Expand(mInitPointCloudArray[i], 1)).segment(0, 3).cast<float>();
    }
}

void tRenderResourceBase::SetPos(const tVector3f &new_pos)
{
    tVector3f now_res = this->mTransform.block(0, 3, 3, 1).cast<float>();
    float diff = (new_pos - now_res).norm();
    if (diff > 1e-6)
    {
        mTransform.block(0, 3, 3, 1) = new_pos.cast<double>();
        ApplyTransform();
    }
}
tVector3f tRenderResourceBase::GetPos() const
{
    return tVector3f(mTransform(0, 3), mTransform(1, 3), mTransform(2, 3));
}

/**
 * \brief           Apply camera pos, to change the point cloud array
*/
void tRenderResourceMesh4View::ApplyCameraPose(const tVector3f &the_first_cam_pos, const tVector3f &cam_focus, const tVector3f &the_first_cam_up)
{
    SIM_ASSERT(mResourcePathList.size() == 4);

    // int num_of_images = 1;
    // for (int i = 0; i < num_of_images; i++)
    for (int i = 0; i < mResourcePathList.size(); i++)
    {
        // tMatrix3f rotmat = cMathUtil::AxisAngleToRotmat(-tVector(0, 1, 0, 0) * i * M_PI / 2).topLeftCorner<3, 3>().cast<float>();
        auto png_path = mResourcePathList[i];
        tVector3f cam_pos = the_first_cam_pos;
        tVector3f cam_up_ = the_first_cam_up;
        tVector3f Z = (cam_pos - cam_focus).normalized();
        // std::cout << "Z = " << Z.transpose() << std::endl;
        tVector3f cam_up = cam_up_.dot(Z) * (-Z) + cam_up_;
        // std::cout << "new up = " << mPngCamUp.transpose() << std::endl;
        cam_up.normalize();
        tVector3f Y = cam_up;
        tVector3f X = (Y.cross(Z)).normalized();
        tMatrix3f R = tMatrix3f::Zero();
        R.col(0) = X;
        R.col(1) = Y;
        R.col(2) = Z;

        int gap = int(this->mPointCloudArray.size() / 4);
        // tMatrix3f rotmat_again = cMathUtil::AxisAngleToRotmat(-tVector(0, 1, 0, 0) * i * M_PI / 2).topLeftCorner<3, 3>().cast<float>();
        tMatrix3f rotmat_again = cMathUtil::AxisAngleToRotmat(tVector(0, 1, 0, 0) * i * M_PI / 2).topLeftCorner<3, 3>().cast<float>();
        for (int j = gap * i; j < gap * (i + 1); j++)
        {
            mPointCloudArray[j] = R * mPointCloudArray[j] + cam_pos;
            mPointCloudArray[j] = rotmat_again * mPointCloudArray[j];
        }
        // std::cout << "png path = " << png_path << std::endl;
    }
}

/**
 * \brief               calculate AABB for rendering resource
*/
void tRenderResourceBase::CalcAABB(tVector3f &aabb_min, tVector3f &aabb_max)
{
    aabb_min = tVector3f::Ones() * std::numeric_limits<float>::max();
    aabb_max = -tVector3f::Ones() * std::numeric_limits<float>::max();
    for (auto &x : mPointCloudArray)
    {
        for (int i = 0; i < 3; i++)
        {
            if (x[i] > aabb_max[i])
            {
                aabb_max[i] = x[i];
            }
            if (x[i] < aabb_min[i])
            {
                aabb_min[i] = x[i];
            }
        }
    }
}

/**
 * \brief           build the rendering resource from string
*/
eRenderResourceType BuildRenderResourceTypeFromStr(std::string cur_type_str)
{
    for (int i = 0; i < eRenderResourceType::NUM_RESOURCE_TYPE; i++)
    {

        if (gRenderResourceTypStr[i] == cur_type_str)
        {
            return static_cast<eRenderResourceType>(i);
        }
    }
    SIM_ERROR("failed to judge the render resource type from string {}", cur_type_str);
    exit(1);
    return eRenderResourceType::NUM_RESOURCE_TYPE;
}

/**
 * \brief           build the rendering string  from resource
*/
std::string BuildStrFromRenderResourceType(eRenderResourceType cur_type)
{
    return gRenderResourceTypStr[cur_type];
}

/**
 * \brief           set the rotate part
*/
void tRenderResourceBase::SetRotAxisAngle(const tVector3f &new_aa)
{
    tVector3f cur_aa = GetRotAxisAngle();
    double diff_norm = (cur_aa - new_aa).norm();
    if (diff_norm > 1e-3)
    {
        tMatrix new_trans = cMathUtil::AxisAngleToRotmat(cMathUtil::Expand(new_aa, 0));
        mTransform.topLeftCorner<3, 3>() = new_trans.topLeftCorner<3, 3>();
        ApplyTransform();
    }
}
/**
 * \brief           get the rotate part
*/
tVector3f tRenderResourceBase::GetRotAxisAngle() const
{
    tVector aa;
    double theta;
    cMathUtil::RotMatToAxisAngle(mTransform, aa, theta);
    aa *= theta;

    return aa.segment(0, 3).cast<float>();
}