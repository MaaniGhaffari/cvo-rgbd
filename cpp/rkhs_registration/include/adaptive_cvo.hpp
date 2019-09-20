/* ----------------------------------------------------------------------------
 * Copyright 2019, Tzu-yuan Lin <tzuyuan@umich.edu>, Maani Ghaffari <maanigj@umich.edu>
 * All Rights Reserved
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

/**
 *  @file   adaptive_cvo.hpp
 *  @author Tzu-yuan Lin, Maani Ghaffari 
 *  @brief  Header file for adaptive contineuous visual odometry
 *  @date   September 20, 2019
 **/

#ifndef ADAPTIVE_CVO_H
#define ADAPTIVE_CVO_H


#include "data_type.h"
#include "LieGroup.h"
#include "pcd_generator.hpp"
#include "../thirdparty/nanoflann.hpp"
#include "../thirdparty/KDTreeVectorOfVectorsAdaptor.h"

#include <vector>
#include <string.h>
#include <iostream>
#include <memory>
#include <utility>
#include <future>
#include <thread>

// #include <pcl/filters/filter.h>
// #include <pcl/point_types.h>
// #include <pcl/point_cloud.h>
// #include <pcl/io/pcd_io.h>
// #include <pcl/common/transforms.h>
// #include <pcl/visualization/cloud_viewer.h>
#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky> 
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/StdVector>
#include <opencv2/core/mat.hpp>
#include <boost/timer/timer.hpp>
// #include <omp.h>
#include <tbb/tbb.h>



using namespace std;
using namespace nanoflann;
using namespace cvo;

namespace acvo{
class acvo{

    private:
        // private variables
        unique_ptr<frame> ptr_fixed_fr;
        unique_ptr<frame> ptr_moving_fr;

        unique_ptr<point_cloud> ptr_fixed_pcd;
        unique_ptr<point_cloud> ptr_moving_pcd;

        int num_fixed;              // target point cloud counts
        int num_moving;             // source point cloud counts
        cloud_t *cloud_x;    // target points represented as a matrix (num_fixed,3)
        cloud_t *cloud_y;    // source points represented as a matrix (num_moving,3)

        float ell;          // kernel characteristic length-scale
        const float ell_init;
        const float ell_min;
        float ell_max;
        double dl;           // changes for ell in each iteration
        double dl_step;
        const float min_dl_step;
        const float max_dl_step;
        const float sigma;        // kernel signal variance (set as std)      
        const float sp_thres;     // kernel sparsification threshold       
        const float c;            // so(3) inner product scale     
        const float d;            // R^3 inner product scale
        const float c_ell;        // kernel characteristic length-scale for color kernel
        const float c_sigma;      // kernel signal variance for color kernel
        const float c_sp_thres;
        const int MAX_ITER;       // maximum number of iteration
        const float eps;          // the program stops if norm(omega)+norm(v) < eps
        const float eps_2;        // threshold for se3 distance
        const float min_step;     // minimum step size for integration
        float step;         // integration step size

        Eigen::Matrix3f R;   // orientation
        Eigen::Vector3f T;   // translation
        Eigen::SparseMatrix<float,Eigen::RowMajor> A;      // coefficient matrix, represented in sparse
        Eigen::SparseMatrix<float,Eigen::RowMajor> Axx;      // coefficient matrix, represented in sparse
        Eigen::SparseMatrix<float,Eigen::RowMajor> Ayy;      // coefficient matrix, represented in sparse
        Eigen::Vector3f omega;  // so(3) part of twist
        Eigen::Vector3f v;      // R^3 part of twist

        // variables for cloud manipulations
        typedef Eigen::Triplet<float> Trip_t;
        tbb::concurrent_vector<Trip_t> A_trip_concur;

        int frame_id = 0;

    public:
        // public variables
        bool init;          // initialization indicator
        int iter;           // final iteration for display
        Eigen::Affine3f transform;  // transformation matrix
        Eigen::Affine3f prev_transform;
        Eigen::Affine3f accum_transform;
        
    private:
        // private functions
        
        /**
         * @brief a polynomial root finder
         * @param coef: coefficeints of the polynomial in descending order
         *              ex. [2,1,3] for 2x^2+x+3
         * @return roots: roots of the polynomial in complex number. ex. a+bi
         */
        inline Eigen::VectorXcf poly_solver(const Eigen::VectorXf& coef);

        /**
         * @brief calculate the se3 distance for giving R and T
         * @return d: se3 distance of given R and T
         */
        inline float dist_se3(const Eigen::Matrix3f& R, const Eigen::Vector3f& T);

        /**
         * @brief update transformation matrix
         */
        inline void update_tf();


        /**
         * @brief isotropic (same length-scale for all dimensions) squared-exponential kernel
         * @param l: kernel characteristic length-scale, aka acvo.ell
         * @prarm s2: signal variance, square of acvo.sigma
         * @return k: n-by-m kernel matrix 
         */
        void se_kernel(point_cloud* cloud_a, point_cloud* cloud_b, \
                       cloud_t* cloud_a_pos, cloud_t* cloud_b_pos,\
                       Eigen::SparseMatrix<float,Eigen::RowMajor>& A_temp);

        /**
         * @brief computes the Lie algebra transformation elements
         *        twist = [omega; v] will be updated in this function
         */
        void compute_flow(); 
        

        /**
         * @brief compute the integration step size
         *        step will be updated in  this function
         */
        void compute_step_size();

        /**
         * @brief transform cloud_y for current update
         */
        void transform_pcd();

        void write_pcl_point_cloud_to_disk(point_cloud* ptr_pcd, string pcd_pth, string pcd_dso_pth);

    public:
        // public funcitons

        // constructor and destructor
        acvo();
        ~acvo();

        /**
         * @brief compute function inner product and return a scalar
         */
        float function_inner_product(point_cloud* cloud_a, point_cloud* cloud_b);

        /**
         * @brief initialize new point cloud and extract pcd as matrices
         */
        void set_pcd(const int dataset_seq,const cv::Mat& RGB_img,const cv::Mat& dep_img, string pcd_pth, string pcd_dso_pth);

        /**
         * @brief align two rgbd pointcloud
         *        the function will iterate MAX_ITER times unless break conditions are met
         */
        void align();

        /**
         * @brief run acvo
         */ 
        void run_cvo(const int dataset_seq,const cv::Mat& RGB_img,const cv::Mat& dep_img, string pcd_pth, string pcd_dso_pth);
};
}
#endif  // ADAPTIVE_CVO
