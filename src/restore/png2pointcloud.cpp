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

#include "utils/OpenCVUtil.h"
void cPng2PointCloud::ResourceWhole(const tMatrixXf &depth_png, double cam_vfov_deg, int &num_of_pts, tEigenArr<tVector3f> &pt_coords)
{

    int height = depth_png.rows();
    int width = depth_png.cols();
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
            float depth = depth_png(row, col);
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
#include "utils/LogUtil.h"
void cPng2PointCloud::ResourceWindow(const tMatrixXf &img, double cam_vfov_deg, const tVector2i &whole_image_size, const tVector2i &window_start, int &num_of_pts, tEigenArr<tVector3f> &pt_coords)
{
    int raw_height = whole_image_size[0];
    int raw_width = whole_image_size[1];
    int window_st_height = window_start[0];
    int window_st_width = window_start[1];
    int window_size_height = img.rows();
    int window_size_width = img.cols();

    double cam_hfov_deg = calculate_hfov(raw_width, raw_height, cam_vfov_deg);
    num_of_pts = 0;

    double cam_hfov_rad = cam_hfov_deg * gDegreesToRadians;
    double cam_vfov_rad = cam_vfov_deg * gDegreesToRadians;

    num_of_pts = 0;
    pt_coords.clear();

    for (int window_row_id = 0; window_row_id < window_size_height; window_row_id++)
    {
        for (int window_col_id = 0; window_col_id < window_size_width; window_col_id++)
        {
            float depth = img(window_row_id, window_col_id);
            if (std::fabs(depth) > 1e-6)
            {

                tVector3f pos = CalcPointPosition(
                    raw_height, raw_width, window_row_id + window_st_height, window_col_id + window_st_width, cam_hfov_rad, cam_vfov_rad, depth);
                pt_coords.push_back(pos);
            }
        }
    }
    num_of_pts = pt_coords.size();
}