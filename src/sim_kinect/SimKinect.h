#pragma once
#include "utils/JsonUtil.h"
#include "render/RenderResource.h"
#include "utils/OpenCVUtil.h"
class cSimKinect
{
public:
    cSimKinect();
    // virtual void Init(const Json::Value &conf);
    virtual void Init();
    virtual tMatrixXf LoadAndCalculate(std::string img_path);
    virtual tMatrixXf Calculate(const tMatrixXf &raw_img);
    int &GetFocalLengthRef();
    float &GetBaselineRef();
    int GetFocalLength() const;
    float GetBaseline() const;
    void SetFocalLength(int focal_length);
    void SetBaseline(float base_line);

protected:
    const float invalid_disp = 99999999.9;
    const float pixel_to_m = 256.0;
    int focal_length = 200;
    float baseline_m = 0.1;
    // std::string mDepthImagePath;
    int rows, cols;
    tVectorXf xx, yy;
    tMatrixXf xp, yp;
    tMatrixXf noised0, noised1;
    tMatrixXf xp_interp, yp_interp;
    cv::Mat raw_depth_m_cv;
    cv::Mat xp_interp_cv;
    cv::Mat yp_interp_cv;
    cv::Mat depth_interp;
    tMatrixXf xf, yf;
    tMatrixXf kinect_pattern_mat;
    tMatrixXf interpolation_map;

    void add_gaussian_shifts(const tMatrixXf &raw_depth_m, tMatrixXf &depth_interp_eigen, float std_var = 1 / 2.0);
    virtual void GenerateDataForImage(int rows, int cols);
    void filter_disp(tMatrixXf disp,
                     tMatrixXf &output_disp,
                     tMatrixXf kinect_pattern,
                     const int h_d,
                     const int w_d,
                     const int h_kp,
                     const int w_kp,
                     const float invalid_disp);
};
SIM_DECLARE_PTR(cSimKinect);