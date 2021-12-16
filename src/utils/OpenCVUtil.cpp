#include "OpenCVUtil.h"
#include "utils/LogUtil.h"
using namespace cv;

uint8_t Clamp(int val)
{
    uint8_t new_val = 0;

    if (val < 0)
        new_val = 0;

    else if (val > 255)
        new_val = 255;
    else
        new_val = val;
    return new_val;
}
#include "utils/FileUtil.h"
cv::Mat cOpencvUtil::LoadGrayscalePng(std::string path)
{
    if(cFileUtil::ExistsFile(path) == false)
    {
        SIM_ERROR("{} doesn't exist ", path);
    }
    cv::Mat img = cv::imread(path, cv::IMREAD_GRAYSCALE);
    return img;
}

tMatrixXf cOpencvUtil::LoadGrayscalePngEigen(std::string path)
{
    cv::Mat res = LoadGrayscalePng(path);
    tMatrixXf png_eigen = cOpencvUtil::ConvertOpencvToEigenMatFloat(res);
    return png_eigen;
}

cv::Mat cOpencvUtil::ConvertFloatArrayToRGBMat(int height, int width, float *array)
{
    Mat new_mat(height, width, CV_8UC3);
    for (int row = 0; row < height; row++)
        for (int col = 0; col < width; col++)
        {
            auto &pixel = new_mat.at<Vec3b>(row, col);
            int buffer_idx = ((height - 1 - row) * width + col) * 3;
            int val_B = array[buffer_idx + 2] * 255.99;
            int val_G = array[buffer_idx + 1] * 255.99;
            int val_R = array[buffer_idx + 0] * 255.99;

            pixel[0] = Clamp(val_B); // B
            pixel[1] = Clamp(val_G); // G
            pixel[2] = Clamp(val_R); // R
        }
    // imshow("hello", new_mat);
    // waitKey(0);
    // exit(1);
    return new_mat;
}

cv::Mat cOpencvUtil::ConvertFloatArrayToGrayscaleMat(int height, int width, int input_channels, float *array)
{
    Mat new_mat(height, width, CV_8UC1);
    for (int row = 0; row < height; row++)
        for (int col = 0; col < width; col++)
        {
            int buf_idx = ((height - 1 - row) * width + col) * input_channels;
            int val = int(array[buf_idx] * 255.99);
            if (val > 255)
                val = 255;
            new_mat.at<uint8_t>(row, col) = val;
        }
    return new_mat;
}

void cOpencvUtil::ConvertRGBMatToFloatArray(const cv::Mat &mat, int &height, int &width, std::vector<float> &array)
{
    height = mat.rows;
    width = mat.cols;
    array.resize(3 * height * width);

    for (int row = 0; row < height; row++)
        for (int col = 0; col < width; col++)
        {
            const auto &pixel = mat.at<Vec3b>(row, col);
            int buffer_idx = ((height - 1 - row) * width + col) * 3;
            array[buffer_idx + 0] = float(pixel[2]) / 255; // R
            array[buffer_idx + 1] = float(pixel[1]) / 255; // G
            array[buffer_idx + 2] = float(pixel[0]) / 255; // B
        }
}

cv::Mat cOpencvUtil::Undistort(const cv::Mat &Mat, const cv::Mat &cameraMatrix, const cv::Mat &distCoeffs)
{
    cv::Mat undistorted;
    cv::undistort(Mat, undistorted, cameraMatrix, distCoeffs);

    return undistorted;
}
#include <opencv2/core/eigen.hpp>
cv::Mat cOpencvUtil::ConvertEigenMatToOpencv(const tMatrixXd &mat)
{
    cv::Mat new_mat;
    cv::eigen2cv(mat, new_mat);
    return new_mat;
}
cv::Mat cOpencvUtil::ConvertEigenVecToOpencv(const tVectorXd &mat)
{
    cv::Mat new_mat;
    cv::eigen2cv(mat, new_mat);
    return new_mat;
}
cv::Mat cOpencvUtil::ConvertEigenMatFloatToOpencv(const tMatrixXf &mat)
{
    cv::Mat new_mat;
    cv::eigen2cv(mat, new_mat);
    return new_mat;
}
cv::Mat cOpencvUtil::ConvertEigenVecFloatToOpencv(const tVectorXf &mat)
{
    cv::Mat new_mat;
    cv::eigen2cv(mat, new_mat);
    return new_mat;
}

tVectorXd cOpencvUtil::ConvertOpencvToEigenVec(const cv::Mat &mat)
{
    tVectorXd new_mat;
    cv::cv2eigen(mat, new_mat);
    return new_mat;
}
tVectorXf cOpencvUtil::ConvertOpencvToEigenVecFloat(const cv::Mat &mat)
{
    tVectorXf new_mat;
    cv::cv2eigen(mat, new_mat);
    return new_mat;
}
tMatrixXd cOpencvUtil::ConvertOpencvToEigenMat(const cv::Mat &mat)
{
    tMatrixXd new_mat;
    cv::cv2eigen(mat, new_mat);
    return new_mat;
}
tMatrixXf cOpencvUtil::ConvertOpencvToEigenMatFloat(const cv::Mat &mat)
{

    tMatrixXf new_mat;
    cv::cv2eigen(mat, new_mat);
    return new_mat;
}
std::string cOpencvUtil::type2str(int type)
{
    std::string r;

    uchar depth = type & CV_MAT_DEPTH_MASK;
    uchar chans = 1 + (type >> CV_CN_SHIFT);

    switch (depth)
    {
    case CV_8U:
        r = "8U";
        break;
    case CV_8S:
        r = "8S";
        break;
    case CV_16U:
        r = "16U";
        break;
    case CV_16S:
        r = "16S";
        break;
    case CV_32S:
        r = "32S";
        break;
    case CV_32F:
        r = "32F";
        break;
    case CV_64F:
        r = "64F";
        break;
    default:
        r = "User";
        break;
    }

    r += "C";
    r += (chans + '0');

    return r;
}

tVector2i cOpencvUtil::GetOpencvMatSize(const cv::Mat &mat)
{
    return tVector2i(mat.rows, mat.cols);
}