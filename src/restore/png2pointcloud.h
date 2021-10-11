#pragma once
#include "utils/MathUtil.h"
#include <memory>

class cPng2PointCloud : std::enable_shared_from_this<cPng2PointCloud>
{
public:
    cPng2PointCloud();

    static void ResourceWhole(const tMatrixXf &img, double cam_vfov_deg, int &num_of_pts, tEigenArr<tVector3f> &pt_coords);
    static void ResourceWindow(const tMatrixXf &img, double cam_vfov_deg, const tVector2i &whole_image_size, const tVector2i &window_start, int &num_of_pts, tEigenArr<tVector3f> &pt_coords);
    
protected:
};