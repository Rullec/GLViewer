#include "SimKinect.h"
#include "utils/OpenCVUtil.h"
#include "utils/TimeUtil.hpp"
cSimKinect::cSimKinect()
{
}

// void MeshGrid(const tVectorXi &xx, const tVectorXi &yy,
//               tMatrixXf &xp, tMatrixXf &yp)
// {
//     int x_size = xx.size(),
//         y_size = yy.size();
//     xp.resize(y_size, x_size);
//     yp.resize(y_size, x_size);
//     xp.rowwise() = xx.cast<float>().transpose();
//     yp.colwise() = yy.cast<float>();
// }
void MeshGrid(const tVectorXf &xx, const tVectorXf &yy,
              tMatrixXf &xp, tMatrixXf &yp)
{
    int x_size = xx.size(),
        y_size = yy.size();
    xp.resize(y_size, x_size);
    yp.resize(y_size, x_size);
    xp.rowwise() = xx.cast<float>().transpose();
    yp.colwise() = yy.cast<float>();
}

void cSimKinect::Init(const Json::Value &conf)
{
    // mDepthImagePath = cJsonUtil::ParseAsString("depth_path", conf);
    // std::cout << "[debug] depth image path = " << mDepthImagePath << std::endl;

    // tMatrixXf raw_depth = cOpencvUtil::LoadGrayscalePngEigen(mDepthImagePath);
    // rows = raw_depth.rows();
    // cols = raw_depth.cols();

    kinect_pattern_mat = cOpencvUtil::LoadGrayscalePngEigen("assets/kinect-pattern_3x3.png");
    // Calculate();
    // Calculate();
    // Calculate();
    // cv::waitKey(0);
}

tMatrixXf GenerateNoisedGausssian(int rows, int cols, float mean, float std)
{
    tMatrixXf res = tMatrixXf::Zero(rows, cols);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
        {
            *(res.data() + j * rows + i) = cMathUtil::RandDoubleNorm(mean, std);
        }
    return res;
}
void cSimKinect::GenerateDataForImage(int rows_, int cols_)
{
    if (rows_ != rows || cols_ != cols)
    {
        rows = rows_;
        cols = cols_;
        xx = tVectorXf::LinSpaced(cols, 0, cols - 1);
        yy = tVectorXf::LinSpaced(rows, 0, rows - 1);
        MeshGrid(xx, yy, xp, yp);
    }
}
void cSimKinect::add_gaussian_shifts(const tMatrixXf &raw_depth_m, tMatrixXf &depth_interp_eigen, float std_var /* = 1 / 2.0*/)
{
    // cTimeUtil::Begin("add_gaussian_shifts_inner");
    // cTimeUtil::Begin("add_gaussian_shifts0");
    int rows = raw_depth_m.rows();
    int cols = raw_depth_m.cols();
    noised0.noalias() = GenerateNoisedGausssian(rows, cols, 0, std_var);
    noised1.noalias() = GenerateNoisedGausssian(rows, cols, 0, std_var);
    // cTimeUtil::End("add_gaussian_shifts0");
    // // cTimeUtil::Begin("add_gaussian_shifts1");
    // std::cout << "xx = " << xx.transpose() << std::endl;
    // std::cout << "yy = " << yy.transpose() << std::endl;

    // // cTimeUtil::End("add_gaussian_shifts1");
    // cTimeUtil::Begin("add_gaussian_shifts2");
    xp_interp.noalias() = (xp + noised0).cwiseMax(0).cwiseMin(cols);
    yp_interp.noalias() = (yp + noised1).cwiseMax(0).cwiseMin(rows);

    raw_depth_m_cv = cOpencvUtil::ConvertEigenMatFloatToOpencv(raw_depth_m);
    xp_interp_cv = cOpencvUtil::ConvertEigenMatFloatToOpencv(xp_interp);
    yp_interp_cv = cOpencvUtil::ConvertEigenMatFloatToOpencv(yp_interp);

    cv::remap(raw_depth_m_cv, depth_interp, xp_interp_cv, yp_interp_cv, cv::INTER_LINEAR);
    depth_interp_eigen.noalias() = cOpencvUtil::ConvertOpencvToEigenMatFloat(depth_interp);
    // cTimeUtil::End("add_gaussian_shifts2");
    // cTimeUtil::End("add_gaussian_shifts_inner");
}

void cSimKinect::filter_disp(tMatrixXf disp,
                             tMatrixXf &output_disp,
                             tMatrixXf kinect_pattern,
                             const int h_d,
                             const int w_d,
                             const int h_kp,
                             const int w_kp,
                             const float invalid_disp)
{
    const int num_threads = omp_get_num_procs();
    const int size_filter = 9;
    int center = int(size_filter / 2.f);
    tMatrixXf xf, yf;
    auto x = tVectorXf::LinSpaced(size_filter, 0, size_filter - 1).transpose(),
         y = tVectorXf::LinSpaced(size_filter, 0, size_filter - 1).transpose();

    MeshGrid(x, y, xf, yf);

    xf.array() -= center;
    yf.array() -= center;

    tMatrixXf sqr_radius = (xf.array().square() + yf.array().square()); // (x**2 + y **2)
    // std::cout << "sqr_radius = " << sqr_radius << std::endl;
    tMatrixXf vals = sqr_radius * std::pow(1.2, 2); // sqr * 1.2 ** 2;
    // std::cout << "vals first = " << vals << std::endl;
    float *vals_ptr = vals.data();

#pragma omp parallel for num_threads(2 * num_threads - 1)
    for (int i = 0; i < vals.size(); i++)
    { // vals[vals == 0] = 1;
        if (vals_ptr[i] == 0)
            vals_ptr[i] = 1;
    }
    // std::cout << "removed vals = " << vals << std::endl;

    tMatrixXf weights = vals.array().cwiseInverse(); // weigths = 1 / vals;
    // std::cout << "weights_ = " << weights << std::endl;
    tMatrixXf fill_weights = (sqr_radius.array() + float(1)).array().cwiseInverse(); // weights = 1 / (1 + sqr_radius);
    float *fill_weights_ptr = fill_weights.data();

    // std::cout << "fill_weights first = " << fill_weights << std::endl;
#pragma omp parallel for num_threads(2 * num_threads - 1)
    for (int i = 0; i < fill_weights.size(); i++)
    { // vals[vals == 0] = 1;
        if (sqr_radius.data()[i] > size_filter)
            fill_weights_ptr[i] = -1;
    }

    // std::cout << "fill_weights = " << fill_weights << std::endl;
    // exit(1);
    int lim_rows = h_d < h_kp ? h_d - size_filter : h_kp - size_filter;
    int lim_cols = w_d < w_kp ? w_d - size_filter : w_kp - size_filter;

    float window_inlier_distance_ = 0.1f;

    output_disp.resize(h_d, w_d);
    output_disp.fill(invalid_disp);
    if (interpolation_map.rows() != h_d || interpolation_map.cols() != w_d)
    {
        interpolation_map.noalias() = tMatrixXf::Zero(h_d, w_d);
    }
    else
    {
        interpolation_map.setZero();
    }

#pragma omp parallel for collapse(2)
    for (int r = 0; r < lim_rows; r++)
    {
        for (int c = 0; c < lim_cols; c++)
        {
            if (kinect_pattern(r + center, c + center) > 0)
            {
                tMatrixXf window = disp.block(r, c, size_filter, size_filter);
                tMatrixXf dot_window = kinect_pattern.block(r, c, size_filter, size_filter);
                float *window_ptr = window.data();
                float *dot_window_ptr = dot_window.data();

                std::vector<float> valid_dots, valid_dots_;
                assert(window.size() == dot_window.size());
                for (int i = 0; i < window.size(); i++)
                {
                    if (window_ptr[i] < invalid_disp)
                    {
                        valid_dots.emplace_back(dot_window_ptr[i]);
                        valid_dots_.emplace_back(window_ptr[i]);
                    }
                }

                tVectorXf valid_dots_v = Eigen::Map<tVectorXf, Eigen::Unaligned>(valid_dots.data(), valid_dots.size());
                tVectorXf valid_dots_v_ = Eigen::Map<tVectorXf, Eigen::Unaligned>(valid_dots_.data(), valid_dots_.size());
                float n_valids = valid_dots_v.sum() / 255.f;
                float n_thresh = dot_window.sum() / 255.f;

                // std::cout << "[first] n_valids = " << n_valids << " n_thresh = " << n_thresh << " at " << r << " " << c << std::endl;
                if (n_valids > n_thresh / 1.2f)
                {
                    float mean = valid_dots_v_.mean();
                    tMatrixXf diffs = (window.array() - mean).cwiseAbs();
                    tMatrixXf diffs_ = diffs.cwiseProduct(weights);
                    float *diffs_ptr = diffs_.data();

                    tMatrixXf tmp_mat0, tmp_mat1;
                    tmp_mat0.resize(size_filter, size_filter);
                    tmp_mat1.resize(size_filter, size_filter);
                    float *tmp_mat0_ptr = tmp_mat0.data();
                    float *tmp_mat1_ptr = tmp_mat1.data();

                    for (int i = 0; i < window.size(); i++)
                    {
                        if (window_ptr[i] < invalid_disp)
                            tmp_mat0_ptr[i] = dot_window_ptr[i];
                        else
                            tmp_mat0_ptr[i] = 0;

                        if (diffs_ptr[i] < window_inlier_distance_)
                            tmp_mat1_ptr[i] = 1;
                        else
                            tmp_mat1_ptr[i] = 0;
                    }

                    tMatrixXf cur_valid_dots = tmp_mat0.cwiseProduct(tmp_mat1);
                    float n_valids = cur_valid_dots.sum() / 255.f;
                    // if (r == 6 && c == 39)
                    // {
                    //     std::cout << "[stop] n valids = " << n_valids << " at " << r << " " << c << std::endl;
                    //     exit(1);
                    // }
                    if (n_valids > n_thresh / 1.2f)
                    {
                        // std::cout << "[inner] n valids " << n_valids << " at " << r << " " << c << std::endl;
                        float accu = window(center, center);
                        assert(accu <= invalid_disp);
                        // std::cout << "[inner] accu " << accu << std::endl;
                        output_disp(r + center, c + center) = round((accu)*8.0) / 8.0;
                        // std::cout << "[inner] set output_disp[r + center, c + center] = " << output_disp(r + center, c + center) << std::endl;

                        // Eigen::Ref<tMatrixXf> interpolation_map.block(r, c, size_filter, size_filter) = interpolation_map.block(r, c, size_filter, size_filter);
                        // Eigen::Ref<tMatrixXf> output_disp.block(r, c, size_filter, size_filter) = output_disp.block(r, c, size_filter, size_filter);

                        // for (int i = 0; i < size_filter; i++)
                        // {
                        //     for (int j = 0; j < size_filter; j++)
                        //     {
                        //         if (interpolation_map.block(r, c, size_filter, size_filter)(i, j) < fill_weights(i, j) + 1e-6)
                        //         {
                        //             // printf("[inner] [inner] judge condition %.6f < %.6f\n", interpolation_map.block(r, c, size_filter, size_filter)(i, j), fill_weights(i, j));
                        //             interpolation_map.block(r, c, size_filter, size_filter)(i, j) = fill_weights(i, j);
                        //             output_disp.block(r, c, size_filter, size_filter)(i, j) = output_disp.block(r, c, size_filter, size_filter)(center, center);
                        //             // printf("[inner] [inner] set interpolation_map (%d, %d) = %.1f\n", i, j, interpolation_map.block(r, c, size_filter, size_filter)(i, j));
                        //             // printf("[inner] [inner] set output_disp (%d, %d) = %.1f\n", i, j, output_disp.block(r, c, size_filter, size_filter)(i, j));
                        //         }
                        //     }
                        // }
                        // std::cout << "[inner] after interpolation_map sum = " << interpolation_map.block(r, c, size_filter, size_filter).sum() << std::endl;
                        // std::cout << "[inner] after output_disp sum = " << output_disp.block(r, c, size_filter, size_filter).sum() << std::endl;

                        float *interpolation_window_ptr = interpolation_map.data();
                        float *output_disp_tep_ptr = output_disp.data();

                        for (int i = 0; i < size_filter; i++)
                        {
                            for (int j = 0; j < size_filter; j++)
                            {
                                int idx = (c + i) * h_d + r + j;
                                int idx_ = i * size_filter + j;
                                if (interpolation_window_ptr[idx] < fill_weights_ptr[idx_])
                                {
                                    interpolation_window_ptr[idx] = fill_weights_ptr[idx_];
                                    output_disp_tep_ptr[idx] = output_disp(r + center, c + center);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

tMatrixXf cSimKinect::Calculate(const tMatrixXf &raw_depth)
{
    int h = raw_depth.rows();
    int w = raw_depth.cols();
    GenerateDataForImage(h, w);

    float scale_factor = 256;
    tMatrixXf raw_depth_m = raw_depth / scale_factor;
    // 2. add gaussian shift
    tMatrixXf depth_interp, output_disp;
    depth_interp = raw_depth_m;
    // cTimeUtil::Begin("add_gaussian_shifts");
    add_gaussian_shifts(raw_depth_m, depth_interp);
    // cTimeUtil::End("add_gaussian_shifts");

    // 3. filter disp
    // 1. ---------------------------------
    cTimeUtil::Begin("filter_disp");
    tMatrixXf disp_ = (depth_interp.array() + float(1e-8)).cwiseAbs().cwiseInverse() * (focal_length * baseline_m);
    // introduce the quantization error
    tMatrixXf depth_f = (disp_ * 8.0f).array().round() / 8.0f; // round(disp * 8) / 8

    // 2. ---------------------------------
    filter_disp(depth_f, output_disp, kinect_pattern_mat, h, w, kinect_pattern_mat.rows(), kinect_pattern_mat.cols(), this->invalid_disp);
    cTimeUtil::End("filter_disp");

    // 3. -----------------------------------
    Eigen::MatrixXf depth = (output_disp.array() + float(1e-8)).cwiseInverse() * (focal_length * baseline_m);

    // filter the valid depth value.
    const int num_threads = omp_get_num_procs();
    int zero_count = 0;
#pragma omp parallel for num_threads(2 * num_threads - 1)
    for (int i = 0; i < h * w; i++)
    {
        if (std::fabs(output_disp.data()[i] - invalid_disp) < 1e-6)
        {
            depth.data()[i] = 0;
            // std::cout << "set pixel " << i << " to 0\n";
            zero_count += 1;
        }
    }
    // std::cout << "zero count = " << zero_count << std::endl;
    tMatrixXf noisy_depth = ((depth.array() * scale_factor).round().cwiseInverse() * 35130).round().cwiseInverse() * 35130;
    return noisy_depth;
    // noisy_depth *= pixel_to_m;

    // depth_aug *= scale_factor;
    // # The depth here needs to converted to cms so scale factor is introduced
    // # though often this can be tuned from [100, 200] to get the desired banding / quantisation effects
    // noisy_depth = (35130 /
    //                np.round((35130 / np.round(depth * scale_factor)) +
    //                         np.random.normal(size = (h, w)) *
    //                             (1.0 / 6.0) +
    //                         0.5)) /
    //               scale_factor

    //                   noisy_depth = noisy_depth * 5000.0 noisy_depth = noisy_depth.astype('uint16')

    // cTimeUtil::End("filter_disp");

    // cv::imshow("result" + std::to_string(cMathUtil::RandDouble()), cOpencvUtil::ConvertEigenMatFloatToOpencv(noisy_depth));
    // cv::imshow("raw" + std::to_string(cMathUtil::RandDouble()), cOpencvUtil::ConvertEigenMatFloatToOpencv(raw_depth / pixel_to_m));
    // cv::waitKey(0);
    return noisy_depth;
}

tMatrixXf cSimKinect::Calculate(std::string img_path)
{
    // cTimeUtil::Begin("calc");
    // 1. read depth image
    tMatrixXf raw_depth = cOpencvUtil::LoadGrayscalePngEigen(img_path);
    // std::cout << "raw depth mean = " << raw_depth.cwiseAbs().mean() << std::endl;
    raw_depth = Calculate(raw_depth);
    std::cout << "[sim_kinect] calc: focal length = " << this->focal_length << std::endl;
    std::cout << "[sim_kinect] calc: baseline = " << this->baseline_m << std::endl;
    // cTimeUtil::End("calc");
    // std::cout << "new depth mean = " << raw_depth.cwiseAbs().mean() << std::endl;
    return raw_depth;
}

int &cSimKinect::GetFocalLength()
{
    return focal_length;
}

float &cSimKinect::GetBaseline()
{
    return baseline_m;
}