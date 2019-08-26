/* ----------------------------------------------------------------------------
 * Copyright 2019, Tzu-yuan Lin <tzuyuan@umich.edu>, Maani Ghaffari <maanigj@umich.edu>
 * All Rights Reserved
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

/**
 *  @file   data_type.h
 *  @author Tzu-yuan Lin, Maani Ghaffari 
 *  @brief  Data type definition
 *  @date   August 4, 2019
 **/
#ifndef DATA_TYPE_H
#define DATA_TYPE_H

#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <tbb/concurrent_vector.h>

#define PYR_LEVELS 3

namespace cvo{

struct frame{

    int frame_id;

    int h;        // height of the image without downsampling
    int w;        // width of the image without downsampling

    cv::Mat image;
    cv::Mat intensity;
    cv::Mat depth;
    

    Eigen::Vector3f* dI;    // flattened image gradient, (w*h,3). 0: magnitude, 1: dx, 2: dy
    Eigen::Vector3f* dI_pyr[PYR_LEVELS];  // pyramid for dI. dI_pyr[0] = dI
    float* abs_squared_grad[PYR_LEVELS];  // pyramid for absolute squared gradient (dx^2+dy^2)

    // destructor
    ~frame(){
        // release memories
        for(int i=0;i<PYR_LEVELS;i++){
			delete[] dI_pyr[i];
			delete[] abs_squared_grad[i];
		}
    }
};

struct point_cloud{

    int num_points;

    std::vector<Eigen::Vector3f> positions;  // points position. x,y,z
    Eigen::Matrix<float, Eigen::Dynamic, 5> features;   // features are rgb dx dy
};



}

#endif // DATA_TYPE_H