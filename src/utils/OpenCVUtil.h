#pragma once
#include <opencv2/opencv.hpp>
#include "utils/MathUtil.h"
class cOpencvUtil
{
public:
    static cv::Mat LoadGrayscalePng(std::string path);
    static tMatrixXf LoadGrayscalePngEigen(std::string path);
    static cv::Mat ConvertFloatArrayToRGBMat(int height, int width, float *array);
    static cv::Mat ConvertFloatArrayToGrayscaleMat(int height, int width, int input_channels, float *array);
    static void ConvertRGBMatToFloatArray(const cv::Mat &mat, int &height, int &width, std::vector<float> &array);
    static cv::Mat Undistort(const cv::Mat &Mat, const cv::Mat &cameraMatrix, const cv::Mat &distCoeffs);

    static cv::Mat ConvertEigenMatToOpencv(const tMatrixXd &mat);
    static cv::Mat ConvertEigenVecToOpencv(const tVectorXd &mat);
    static cv::Mat ConvertEigenMatFloatToOpencv(const tMatrixXf &mat);
    static cv::Mat ConvertEigenVecFloatToOpencv(const tVectorXf &mat);

    static tVectorXd ConvertOpencvToEigenVec(const cv::Mat &mat);
    static tVectorXf ConvertOpencvToEigenVecFloat(const cv::Mat &mat);
    static tMatrixXd ConvertOpencvToEigenMat(const cv::Mat &mat);
    static tMatrixXf ConvertOpencvToEigenMatFloat(const cv::Mat &mat);
    static std::string type2str(int type);
    static tVector2i GetOpencvMatSize(const cv::Mat &mat);
};