#include "GLUtil.h"
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include "utils/MathUtil.h"
void ConvertImageToTexutre(GLuint &texture)
{
    texture = 0;
    cv::Mat image = cv::imread("now.png");
    if (image.empty())
    {
        std::cout << "image empty" << std::endl;
    }
    else
    {
        cv::flip(image, image, 0);
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        // Set texture clamping method
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

        glTexImage2D(GL_TEXTURE_2D,    // Type of texture
                     0,                // Pyramid level (for mip-mapping) - 0 is the top level
                     GL_RGBA,          // Internal colour format to convert to
                     image.cols,       // Image width  i.e. 640 for Kinect in standard mode
                     image.rows,       // Image height i.e. 480 for Kinect in standard mode
                     0,                // Border width in pixels (can either be 1 or 0)
                     GL_BGR,           // Input image format (i.e. GL_RGB, GL_RGBA, GL_BGR etc.)
                     GL_UNSIGNED_BYTE, // Image data type
                     image.ptr());     // The actual image data itself

        glGenerateMipmap(GL_TEXTURE_2D);
    }
}

void UpdateValue(float *data, int height, int width, int comp = 3)
{
    for (int row = 0; row < height; row++)
    {
        for (int col = 0; col < width; col++)
        {
            int bias = (row * width + col) * 3;
            // data[bias + 0] = 0;
            data[bias + 1] += (cMathUtil::RandDouble() - 0.5) / 10;
            data[bias + 2] += (cMathUtil::RandDouble() - 0.5) / 10;
        }
    }
}

void ConvertDepthImageToRGB(const tMatrixXi &depth_image,
                            tMatrixXf &R, tMatrixXf &G, tMatrixXf &B)
{
    R = depth_image.cast<float>();
    G.noalias() = G;
    B.noalias() = B;
    float max_amp = 1000;
    R /= max_amp;
    G /= max_amp;
    B /= max_amp;
}
