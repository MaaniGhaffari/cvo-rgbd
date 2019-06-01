/* ----------------------------------------------------------------------------
 * Copyright 2019, Tzu-yuan Lin <tzuyuan@umich.edu>, Maani Ghaffari <maanigj@umich.edu>
 * All Rights Reserved
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

/**
 *  @file   rkhs_se3.hpp
 *  @author Tzu-yuan Lin, Maani Ghaffari 
 *  @brief  Header file for contineuous visual odometry rkhs_se3 registration
 *  @date   May 31, 2019
 **/

#ifndef RKHS_SE3_H
#define RKHS_SE3_H

// #define EIGEN_USE_MKL_ALL

#include <LieGroup.h>
#include <iostream>
#include <pcl/filters/filter.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
#include <string.h>
#include <pcl/common/transforms.h>
#include <pcl/visualization/cloud_viewer.h>
#include <boost/timer/timer.hpp>

using namespace std;

class rkhs_se3{

    private:
        // private variables
        pcl::PointCloud< pcl::PointXYZRGB > fixed;  // target point cloud, remain unchanged
        pcl::PointCloud< pcl::PointXYZRGB > moving; // source point cloud, remain unchanged
        pcl::PointCloud< pcl::PointXYZRGB > moved; // source point cloud that changes in each iteration

        int num_fixed;              // target point cloud counts
        int num_moving;             // source point cloud counts
        Eigen::Matrix<float, Eigen::Dynamic, 3> cloud_x;    // target points represented as a matrix (num_fixed,3)
        Eigen::Matrix<float, Eigen::Dynamic, 3> cloud_y;    // source points represented as a matrix (num_moving,3)

        float ell;          // kernel characteristic length-scale
        float sigma;        // kernel signal variance (set as std)      
        float sp_thres;     // kernel sparsification threshold       
        float c;            // so(3) inner product scale     
        float d;            // R^3 inner product scale
        float color_scale;  // color space inner product scale
        int MAX_ITER;       // maximum number of iteration
        float eps;          // the program stops if norm(omega)+norm(v) < eps
        float eps_2;        // threshold for se3 distance
        float min_step;     // minimum step size for integration
        float step;         // integration step size

        Eigen::Matrix3f R;   // orientation
        Eigen::Vector3f T;   // translation
        Eigen::MatrixXf K;   // isotropic SE kernel
        Eigen::SparseMatrix<float,Eigen::RowMajor> A;      // coefficient matrix, represented in sparse
        vector< vector<int> > group_ids;    // nonzero indicies of coefficient matrix A
        Eigen::Vector3f omega;  // so(3) part of twist
        Eigen::Vector3f v;      // R^3 part of twist

    public:
        // public variables
        int iter;           // final iteration for display
        Eigen::Affine3f transform;  // transformation matrix

    private:
        // private functions
        
        /**
         * @brief remove nan points in pcd
         *        this function is currently unused
         */
        void remove_nan_points();

        /**
         * @brief a polynomial root finder
         * @param coef: coefficeints of the polynomial in descending order
         *              ex. [2,1,3] for 2x^2+x+3
         * @return roots: roots of the polynomial in complex number. ex. a+bi
         */
        Eigen::VectorXcf poly_solver(Eigen::VectorXf& coef);

        /**
         * @brief calculate the se3 distance for giving R and T
         * @return d: se3 distance of given R and T
         */
        float dist_se3(Eigen::Matrix3f& R, Eigen::Vector3f& T);

        /**
         * @brief update transformation matrix
         */
        void update_tf();

        /**
         * @brief compute color inner product of ith row in fixed and jth row in moving
         * @param i: index of desired row in fixed
         * @param j: indxe of desired row in moving
         * @return CI: the inner product
         */
        float color_inner_product(const int i, const int j);
        
        /**
         * @brief isotropic (same length-scale for all dimensions) squared-exponential kernel
         * @param l: kernel characteristic length-scale, aka rkhs_se3.ell
         * @prarm s2: signal variance, square of rkhs_se3.sigma
         * @return k: n-by-m kernel matrix 
         */
        void se_kernel(const float l, const float s2);

        /**
         * @brief computes the Lie algebra transformation elements
         *        twist = [omega; v] will be updated in this function
         */
        void compute_flow();
        
        /**
         * @brief compute the integration step size
         *        step will be updated in this function
         */
        void compute_step_size();

    public:
        // public funcitons

        // constructor and destructor
        rkhs_se3();
        ~rkhs_se3();

        /**
         * @brief initialize new point cloud and extract pcd as matrices
         */
        void set_pcd(string fixed_pth, string moving_pth);

        /**
         * @brief align two rgbd pointcloud
         *        the function will iterate MAX_ITER times unless break conditions are met
         */
        void align();
};
#endif  // RKHS_SE3_H
