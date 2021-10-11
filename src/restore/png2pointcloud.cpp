#include "png2pointcloud.h"
#include <cmath>

float calculate_hfov(int width, int height, double vfov_deg)
{
    double vfov_rad = gDegreesToRadians * vfov_deg;
    double hfov_rad = 2 * std::atan(width * std::tan(vfov_rad / 2) / float(height));
    double hfov_deg = hfov_rad * gRadiansToDegrees;
    return hfov_deg;
}

tVector3f CalcPointPosition(int height, int width, int row_id, int col_id, double hfov_rad, double vfov_rad, float depth)
{
    float angle0 = hfov_rad * float(col_id) / width - hfov_rad / 2;
    float angle1 = -vfov_rad * float(row_id) / height + vfov_rad / 2;
    double x = std::tan(angle0);
    double y = std::tan(angle1);
    tVector3f dir = tVector3f(x, y, -1).normalized();
    double length = -depth / (dir.z());
    tVector3f point_pos = dir * length;
    return point_pos;
}

void cPng2PointCloud::Resource(const tMatrixXf &depth_map, double cam_vfov_deg, int &num_of_pts, tEigenArr<tVector3f> &pt_coords)
{
    int height = depth_map.rows();
    int width = depth_map.cols();
    double cam_hfov_deg = calculate_hfov(width, height, cam_vfov_deg);
    num_of_pts = 0;

    double cam_hfov_rad = cam_hfov_deg * gDegreesToRadians;
    double cam_vfov_rad = cam_vfov_deg * gDegreesToRadians;
    /*
    hfov = 2 atan(w * tan (vfov/2) / h)
    */
    // tEigenArr<tVector3f> pos_lst(0);
    pt_coords.reserve(height * width);
    for (int row = 0; row < height; row++)
    {
        for (int col = 0; col < width; col++)
        {
            float depth = depth_map(row, col);
            if (std::fabs(depth) < 1e-4)
                continue;
            // get point position
            tVector3f point_pos = CalcPointPosition(
                height, width, row, col, cam_hfov_rad, cam_vfov_rad, depth);
            pt_coords.push_back(point_pos);
        }
    }

    // cover the result
    num_of_pts = pt_coords.size();
    // point_coords.resize(3 * num_of_pts);
    // for (int i = 0; i < num_of_pts; i++)
    // {
    //     point_coords.segment(3 * i, 3 * (i + 1)).noalias() = pos_lst[i];
    // }
}
