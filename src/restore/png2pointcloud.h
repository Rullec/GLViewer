#pragma once
#include "utils/MathUtil.h"
#include <memory>

class cPng2PointCloud : std::enable_shared_from_this<cPng2PointCloud>
{
public:
    cPng2PointCloud();
    void Resource(const tMatrixXf &depth_map, double cam_vfov_deg, int &num_of_pts, tEigenArr<tVector3f> & pt_coords);
};