/* ----------------------------------------------------------------------------
 * Copyright 2019, Tzu-yuan Lin <tzuyuan@umich.edu>, Maani Ghaffari <maanigj@umich.edu>
 * All Rights Reserved
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

/**
 *  @file   rkhs_se3.cpp
 *  @author Tzu-yuan Lin, Maani Ghaffari 
 *  @brief  Source file for contineuous visual odometry rkhs_se3 registration
 *  @date   May 31, 2019
 **/

#include "rkhs_se3.hpp"

rkhs_se3::rkhs_se3():
    // initialize parameters
    ell(0.15),             // kernel characteristic length-scale
    sigma(0.1),            // kernel signal variance (set as std)      
    sp_thres(1.0e-3),      // kernel sparsification threshold       
    c(7.0),                  // so(3) inner product scale     
    d(7.0),                  // R^3 inner product scale
    color_scale(1.0e-5),   // color space inner product scale
    MAX_ITER(2000),        // maximum number of iteration
    min_step(2*1.0e-1),    // minimum integration step
    eps(5*1.0e-5),         // threshold for stopping the function
    eps_2(1.0e-5),         // threshold for se3 distance
    // A(num_fixed,num_moving),           // initialize sparse coefficient matrix
    R(Eigen::Matrix3f::Identity(3,3)), // initialize rotation matrix to I
    T(Eigen::Vector3f::Zero()),        // initialize translation matrix to zeros
    transform(Eigen::Affine3f::Identity())    // initialize transformation to I
{
    
}

rkhs_se3::~rkhs_se3(){
}

void rkhs_se3::remove_nan_points(){
    // remove nan points
    vector<int> nan_ids;
    pcl::removeNaNFromPointCloud(fixed, fixed, nan_ids);
    nan_ids.clear();
    pcl::removeNaNFromPointCloud(moving,moving,nan_ids);
}

Eigen::VectorXcf rkhs_se3::poly_solver(Eigen::VectorXf& coef){
    // extract order
    int order = coef.size()-1;
    Eigen::VectorXcf roots;
    
    // create M = diag(ones(n-1,1),-1)
    Eigen::MatrixXf M = Eigen::MatrixXf::Zero(order,order);
    M.bottomLeftCorner(order-1,order-1) = Eigen::MatrixXf::Identity(order-1,order-1);
    
    // M(1,:) = -p(2:n+1)./p(1)
    Eigen::VectorXf coefdiv = coef/coef(0);
    M.row(0) = -(coef/coef(0)).segment(1,order).transpose();

    // eigen(M) and get the answer
    roots = M.eigenvalues();

    return roots;
}

float rkhs_se3::dist_se3(Eigen::Matrix3f& R, Eigen::Vector3f& T){
    // create transformation matrix
    Eigen::Matrix4f temp_transform = Eigen::Matrix4f::Identity();
    temp_transform.block<3,3>(0,0)=R;
    temp_transform.block<3,1>(0,3)=T;
    
    // distance = frobenius_norm(logm(trans))
    float d = temp_transform.log().norm();
    
    return d;
}

void rkhs_se3::update_tf(){
    // transform = [R', -R'*T; 0,0,0,1]
    transform.matrix().block<3,3>(0,0) = R.transpose();
    transform.matrix().block<3,1>(0,3) = -R.transpose()*T;
}


float rkhs_se3::color_inner_product(const int i, const int j){
    // compute color inner product for fixed[i] and moving[j]
    float CI = color_scale*(float(fixed.points[i].r)*float(moving.points[j].r)+
                    float(fixed.points[i].g)*float(moving.points[j].g)+
                    float(fixed.points[i].b)*float(moving.points[j].b));
    return CI;
}

void rkhs_se3::se_kernel(const float l, const float s2){
    /**
     * @brief isotropic (same length-scale for all dimensions) squared-exponential kernel
     * @param cloud_x: n-by-d target point cloud as matrix [x0,y0,z0; x1,y1,x1; ...]
     * @param cloud_y: m-by-d source point cloud as matrix
     *                 n, m are number of observation, d is the dimension.
     * @return A: sparsified coefficient matrix 
     */
    
    // convert k threshold to d2 threshold (so that we only need to calculate k when needed)
    float d2_thres = -2.0*l*l*log(sp_thres/s2);
    typedef Eigen::Triplet<float> T;
    std::vector<T> tripletList;
    tripletList.reserve(10000);
    // loop through points
    for(int i=0; i<num_fixed; ++i){
        for(int j=0; j<num_moving; ++j){
            // d2 = (x-y)^2
            float d2 = (cloud_x.row(i)-cloud_y.row(j)).squaredNorm();
            if(d2<d2_thres){
                float k = s2*exp(-d2/(2.0*l*l));
                float ci = color_inner_product(i,j);
                tripletList.push_back(T(i,j,ci*k));
            }
        }
    }
    // form A
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    A.makeCompressed();
}

void rkhs_se3::compute_flow(){
    // compute SE kernel
    se_kernel(ell, sigma*sigma);

    // some initialization of the variables
    omega = Eigen::Vector3f::Zero();
    v = Eigen::Vector3f::Zero();
    group_ids.clear();
    group_ids.reserve(10000);

    // loop through points in cloud_x
    for(int i=0; i<num_fixed; ++i){
        // initialize reused varaibles
        vector<int> used_ids;
        int num_used_ids = A.innerVector(i).nonZeros();
        Eigen::MatrixXf Ai = Eigen::MatrixXf::Zero(1,num_used_ids);
        Eigen::MatrixXf cross_xy = Eigen::MatrixXf::Zero(num_used_ids,3);
        Eigen::MatrixXf diff_yx = Eigen::MatrixXf::Zero(num_used_ids,3);
        Eigen::Matrix<float, 1, 3> partial_omega;
        Eigen::Matrix<float, 1, 3> partial_v;

        int j = 0;
        // loop through used ids in ith row
        for(Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(A,i); it; ++it){
            used_ids.push_back(it.col());   // take out current index
            Ai(0,j) = it.value();    // extract current value in A
            cross_xy.row(j) = cloud_x.row(i).cross(cloud_y.row(used_ids[j]));
            diff_yx.row(j) = cloud_y.row(used_ids[j])-cloud_x.row(i);
            ++j;
        }

        partial_omega = 1/c*Ai*cross_xy;
        partial_v = 1/d*Ai*diff_yx;

        // collect this used_ids into group_ids
        group_ids.push_back(used_ids);

        // sum up them to class-wide variable
        omega += partial_omega.transpose();
        v += partial_v.transpose();
    }
}


void rkhs_se3::compute_step_size(){
    // compute skew matrix
    Eigen::Matrix3f omega_hat = skew(omega);
    
    // compute xi*z+v, xi^2*z+xi*v, xi^3*z+xi^2*v, xi^4*z+xi^3*v
    Eigen::MatrixXf xiz(num_moving,cloud_y.cols());
    
    for(int i=0; i<num_moving; i++) 
        xiz.row(i) = omega.transpose().cross(cloud_y.row(i))+v.transpose(); // (xi*z+v)
    
    Eigen::MatrixXf xi2z = ((omega_hat*omega_hat*cloud_y.transpose()).colwise()\
                            +(omega_hat*v)).transpose();    // (xi^2*z+xi*v)
    Eigen::MatrixXf xi3z = ((omega_hat*omega_hat*omega_hat*cloud_y.transpose()).colwise()\
                            +(omega_hat*omega_hat*v)).transpose();  // (xi^3*z+xi^2*v)
    Eigen::MatrixXf xi4z = ((omega_hat*omega_hat*omega_hat*omega_hat*cloud_y.transpose()).colwise()\
                            +(omega_hat*omega_hat*omega_hat*v)).transpose();    // (xi^4*z+xi^3*v)
    
    // perform row-wise norm/dot for coefficients
    Eigen::MatrixXf normxiz2 = xiz.rowwise().norm().array().pow(2.0).matrix();    // norm(xiz).^2
    Eigen::MatrixXf xiz_dot_xi2z = (-xiz*xi2z.transpose()).diagonal();    // dot(-xiz,xi2z) 
    Eigen::MatrixXf epsil_const =  xi2z.rowwise().norm().array().pow(2.0).matrix()\
                                    + 2.0*(xiz*xi3z.transpose()).diagonal(); // norm(xi2z).^2+2*dot(xiz,xi3z)
    
    // initialize coefficients
    float temp_coef = 1/(2.0*ell*ell);   // 1/(2*l^2)
    float B = 0;
    float C = 0;
    float D = 0;
    float E = 0;

    // loops through points in cloud_x to calculate the derivatives
    for(int i=0; i<num_fixed; ++i){
        // initialization
        Eigen::MatrixXf diff_xy(group_ids[i].size(),3);
        Eigen::MatrixXf beta_i(group_ids[i].size(),1);
        Eigen::MatrixXf gamma_i(group_ids[i].size(),1);
        Eigen::MatrixXf delta_i(group_ids[i].size(),1);
        Eigen::MatrixXf epsil_i(group_ids[i].size(),1);
        Eigen::MatrixXf Ai = Eigen::MatrixXf::Zero(1,group_ids[i].size());
        
        int j=0;
        // loop through used index in ith row

        for(Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(A,i); it; ++it){
            int idx = group_ids[i][j];

            // diff_xy = x[i] - y[used_idx[j]]
            diff_xy.row(j) = cloud_x.row(i) - cloud_y.row(idx);    
            // beta_i = -1/l^2 * dot(xiz,diff_xy)
            beta_i.row(j) = -2.0*temp_coef * xiz.row(idx)*diff_xy.row(j).transpose();
            // gamma_i = -1/(2*l^2) * (norm(xiz).^2 + 2*dot(xi2z,diff_xy))
            gamma_i.row(j) = -temp_coef * (normxiz2.row(idx)\
                            + 2.0*xi2z.row(idx)*diff_xy.row(j).transpose());
            // delta_i = 1/l^2 * (dot(-xiz,xi2z) + dot(-xi3z,diff_xy))
            delta_i.row(j) = 2.0*temp_coef * (xiz_dot_xi2z.row(idx)\
                            + (-xi3z.row(idx)*diff_xy.row(j).transpose()));
            // epsil_i = -1/(2*l^2) * (norm(xi2z).^2 + 2*dot(xiz,xi3z) + 2*dot(xi4z,diff_xy))
            epsil_i.row(j) = -temp_coef * (epsil_const.row(idx)\
                            + 2.0*xi4z.row(idx)*diff_xy.row(j).transpose());

            // form Ai
            Ai(0,j) = it.value();

            ++j;
        }
        
        Eigen::MatrixXf beta2_i = beta_i.array().pow(2.0).matrix();   // beta_i.^2
        Eigen::MatrixXf beta3_i = beta2_i.array().cwiseProduct(beta_i.array()).matrix();   // beta_i.^3
        Eigen::MatrixXf beta4_i = beta3_i.array().cwiseProduct(beta_i.array()).matrix();   // beta_i.^4
        Eigen::MatrixXf beta_gamma_i = beta_i.array().cwiseProduct(gamma_i.array()).matrix();   // beta.*gamma
        Eigen::MatrixXf beta_delta_i = beta_i.array().cwiseProduct(delta_i.array()).matrix();   // beta.*delta
        Eigen::MatrixXf beta2_gamma_i = beta2_i.array().cwiseProduct(gamma_i.array()).matrix(); // beta.^2.*delta
        Eigen::MatrixXf gamma2_i = gamma_i.array().pow(2.0).matrix(); // gamma.^2

        // coefficients
        B += (Ai * beta_i)(0);
        C += (Ai * (gamma_i + beta2_i/2.0))(0);
        D += (Ai * (delta_i + beta_gamma_i + beta3_i/6.0))(0);
        E += (Ai * (epsil_i + beta_delta_i +1/2.0*beta2_gamma_i + 1/2.0*gamma2_i + 1/24.0*beta4_i))(0);
    }

    // create polynomial coefficient vector
    Eigen::VectorXf p_coef(4);
    p_coef << 4.0*E,3.0*D,2.0*C,B;
    
    // solve polynomial roots
    Eigen::VectorXcf rc = poly_solver(p_coef);
    
    // find usable step size
    float temp_step = numeric_limits<float>::max();
    for(int i=0;i<rc.real().size();i++)
        if(rc(i,0).real()>0 && rc(i,0).real()<temp_step && rc(i,0).imag()==0)
            temp_step = rc(i,0).real();
    
    // if none of the roots are suitable, use min_step
    step = temp_step==numeric_limits<float>::max()? min_step:temp_step;

    // if step>0.8, just use 0.8 as step
    step = step>0.8 ? 0.8:step;
}

void rkhs_se3::set_pcd(string fixed_pth, string moving_pth){
    // load point clouds
    if(pcl::io::loadPCDFile<pcl::PointXYZRGB>(fixed_pth, fixed)==-1)
        PCL_ERROR("couldn't find ds fixed file");
    if(pcl::io::loadPCDFile<pcl::PointXYZRGB>(moving_pth, moving)==-1)
        PCL_ERROR("couldn't find ds moving file");

    // get total number of points
    num_fixed = fixed.width*fixed.height;
    num_moving = moving.width*moving.height;
    std::cout<<"num fixed: "<<num_fixed<<std::endl;
    std::cout<<"num moving: "<<num_moving<<std::endl;
    // convert pcl to eigen matrices
    cloud_x = fixed.getMatrixXfMap(3,8,0).transpose();
    cloud_y = moving.getMatrixXfMap(3,8,0).transpose();
    
    // initialization of parameters
    A.resize(num_fixed,num_moving); 
    R = Eigen::Matrix3f::Identity(3,3); // initialize rotation matrix to I
    T = Eigen::Vector3f::Zero();        // initialize translation matrix to zeros
    transform = Eigen::Affine3f::Identity();    // initialize transformation to I
}

void rkhs_se3::align(){
    // loop until MAX_ITER
    for(int k=0; k<MAX_ITER; k++){
        // std::cout<<"----------------------------"<<std::endl;
        // boost::timer::cpu_timer timer;
        // update transformation matrix
        update_tf();

        // apply transform to the point cloud
        pcl::transformPointCloud(moving, moved, transform);

        // convert transformed pcd into an eigen matrix
        cloud_y = moved.getMatrixXfMap(3,8,0).transpose();

        // compute omega and v
        compute_flow();

        // compute step size for integrating the flow
        compute_step_size();

        // stop if the step size is too small
        if(omega.norm()<eps && v.norm()<eps){
            iter = k;
            break;
        }

        // stacked omega and v for finding dtrans
        Eigen::VectorXf vec_joined(omega.size()+v.size());
        vec_joined << omega, v;

        // find the change of translation matrix dtrans
        Eigen::MatrixXf dtrans = Exp_SEK3(vec_joined, step);

        // extract dR and dT from dtrans
        Eigen::Matrix3f dR = dtrans.block<3,3>(0,0);
        Eigen::Vector3f dT = dtrans.block<3,1>(0,3);

        // calculate new R and T
        T = R * dT + T;
        R = R * dR;

        // if the se3 distance is smaller than eps2, break
        if(dist_se3(dR,dT)<eps_2){
            iter = k;
            break;
        }

        ell = (k>2)? 0.10:ell;
        ell = (k>9)? 0.06:ell;
        ell = (k>19)? 0.03:ell;

        // std::cout<<"one loop time: "<<timer.format()<<std::endl;
    }
    update_tf();
}