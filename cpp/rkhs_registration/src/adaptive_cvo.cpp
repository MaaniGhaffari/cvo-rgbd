/* ----------------------------------------------------------------------------
 * Copyright 2019, Tzu-yuan Lin <tzuyuan@umich.edu>, Maani Ghaffari <maanigj@umich.edu>
 * All Rights Reserved
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

/**
 *  @file   adaptive_cvo.cpp
 *  @author Tzu-yuan Lin, Maani Ghaffari 
 *  @brief  Source file for adapive contineuous visual odometry
 *  @date   September 20, 2019
 **/

#include "adaptive_cvo.hpp"

namespace acvo{

acvo::acvo():
    // initialize parameters
    init(false),           // initialization indicator
    ptr_fixed_fr(new frame),
    ptr_moving_fr(new frame),
    ptr_fixed_pcd(new point_cloud),
    ptr_moving_pcd(new point_cloud),
    ell_init(0.1),             // kernel characteristic length-scale
    ell(0.1),
    ell_min(0.0391),
    ell_max(0.15),
    dl(0),
    dl_step(0.3),
    min_dl_step(0.05),
    max_dl_step(1),
    sigma(0.1),             // kernel signal variance (set as std)      
    sp_thres(8.315e-3),     // kernel sparsification threshold       
    c(7.0),                 // so(3) inner product scale     
    d(7.0),                 // R^3 inner product scale
    c_ell(0.5),             // kernel characteristic length-scale for color kernel
    c_sigma(1),
    c_sp_thres(8.315e-3),
    MAX_ITER(2000),         // maximum number of iteration
    min_step(2*1.0e-1),     // minimum integration step
    eps(5*1.0e-5),          // threshold for stopping the function
    eps_2(1.0e-5),          // threshold for se3 distance
    R(Eigen::Matrix3f::Identity(3,3)), // initialize rotation matrix to I
    T(Eigen::Vector3f::Zero()),        // initialize translation matrix to zeros
    transform(Eigen::Affine3f::Identity()),    // initialize transformation to I
    prev_transform(Eigen::Affine3f::Identity()),
    accum_transform(Eigen::Affine3f::Identity())
{
}

acvo::~acvo(){
}

inline Eigen::VectorXcf acvo::poly_solver(const Eigen::VectorXf& coef){
    // extract order
    int order = coef.size()-1;
    Eigen::VectorXcf roots;
    
    // create M = diag(ones(n-1,1),-1)
    Eigen::MatrixXf M = Eigen::MatrixXf::Zero(order,order);
    M.bottomLeftCorner(order-1,order-1) = Eigen::MatrixXf::Identity(order-1,order-1);
    
    // M(1,:) = -p(2:n+1)./p(1)
    M.row(0) = -(coef/coef(0)).segment(1,order).transpose();

    // eigen(M) and get the answer
    roots = M.eigenvalues();

    return roots;
}

inline float acvo::dist_se3(const Eigen::Matrix3f& R, const Eigen::Vector3f& T){
    // create transformation matrix
    Eigen::Matrix4f temp_transform = Eigen::Matrix4f::Identity();
    temp_transform.block<3,3>(0,0)=R;
    temp_transform.block<3,1>(0,3)=T;
    
    // distance = frobenius_norm(logm(trans))
    float d = temp_transform.log().norm();
    
    return d;
}

inline void acvo::update_tf(){
    // transform = [R', -R'*T; 0,0,0,1]
    transform.matrix().block<3,3>(0,0) = R.transpose();
    transform.matrix().block<3,1>(0,3) = -R.transpose()*T;
}


void acvo::se_kernel(point_cloud* cloud_a, point_cloud* cloud_b, \
                         cloud_t* cloud_a_pos, cloud_t* cloud_b_pos,\
                         Eigen::SparseMatrix<float,Eigen::RowMajor>& A_temp){
    A_trip_concur.clear();
    const float s2= sigma*sigma;
    const float l = ell;

    // convert k threshold to d2 threshold (so that we only need to calculate k when needed)
    float d2_thres = -2.0*l*l*log(sp_thres/s2);
    float d2_c_thres = -2.0*c_ell*c_ell*log(c_sp_thres/c_sigma/c_sigma);

    /** 
     * kdtreeeeeeeeeeeeeeeeeeeee
     **/
    typedef KDTreeVectorOfVectorsAdaptor<cloud_t, float>  kd_tree_t;

	kd_tree_t mat_index(3 /*dim*/, (*cloud_b_pos), 10 /* max leaf */ );
	mat_index.index->buildIndex();

    // loop through points
    tbb::parallel_for(int(0),cloud_a->num_points,[&](int i){
    // for(int i=0; i<num_fixed; ++i){

        const float search_radius = d2_thres;
		std::vector<std::pair<size_t,float>>  ret_matches;

		nanoflann::SearchParams params;
		//params.sorted = false;

		const size_t nMatches = mat_index.index->radiusSearch(&(*cloud_a_pos)[i](0), search_radius, ret_matches, params);

        Eigen::Matrix<float,NUM_FEATURES,1> feature_a = cloud_a->features.row(i).transpose();

        // for(int j=0; j<num_moving; j++){
        for(size_t j=0; j<nMatches; ++j){
            int idx = ret_matches[j].first;
            float d2 = ret_matches[j].second;
            // d2 = (x-y)^2
            float k = 0;
            float ck = 0;
            float d2_color = 0;
            float a = 0;
            if(d2<d2_thres){
                Eigen::Matrix<float,NUM_FEATURES,1> feature_b = cloud_b->features.row(idx).transpose();
                d2_color = ((feature_a-feature_b).squaredNorm());

                if(d2_color<d2_c_thres){
                    k = s2*exp(-d2/(2.0*l*l));
                    ck = c_sigma*c_sigma*exp(-d2_color/(2.0*c_ell*c_ell));
                    a = ck*k;
                    if (a > sp_thres) A_trip_concur.push_back(Trip_t(i,idx,a));
                }
            }
        }
    });
    // }
    // form A
    A_temp.setFromTriplets(A_trip_concur.begin(), A_trip_concur.end());
    A_temp.makeCompressed();
}


void acvo::compute_flow(){
    // compute SE kernel for Axy
    se_kernel(ptr_fixed_pcd.get(),ptr_moving_pcd.get(),cloud_x,cloud_y,A);

    // compute SE kernel for Axx and Ayy
    se_kernel(ptr_fixed_pcd.get(),ptr_fixed_pcd.get(),cloud_x,cloud_x,Axx);
    se_kernel(ptr_moving_pcd.get(),ptr_moving_pcd.get(),cloud_y,cloud_y,Ayy);

    // some initialization of the variables
    omega = Eigen::Vector3f::Zero();
    v = Eigen::Vector3f::Zero();
    Eigen::Vector3d double_omega = Eigen::Vector3d::Zero(); // this is omega in double precision
    Eigen::Vector3d double_v = Eigen::Vector3d::Zero();
    dl = 0;
    tbb::spin_mutex omegav_lock;
    int sum = 0;

    float ell_3 = ell*ell*ell;

    // loop through points in cloud_x
    tbb::parallel_for(int(0),num_fixed,[&](int i){
    // for(int i=0; i<num_fixed; ++i){
        // initialize reused varaibles
        int num_non_zeros = A.innerVector(i).nonZeros();
        int num_non_zeros_xx = Axx.innerVector(i).nonZeros();
        int num_non_zeros_yy; 
        num_non_zeros_yy = (i<num_moving)?Ayy.innerVector(i).nonZeros():0;
        Eigen::MatrixXf Ai = Eigen::MatrixXf::Zero(1,num_non_zeros);
        Eigen::MatrixXf Axxi = Eigen::MatrixXf::Zero(1,num_non_zeros_xx);
        Eigen::MatrixXf Ayyi = Eigen::MatrixXf::Zero(1,num_non_zeros_yy);
        Eigen::MatrixXf cross_xy = Eigen::MatrixXf::Zero(num_non_zeros,3);
        Eigen::MatrixXf diff_yx = Eigen::MatrixXf::Zero(num_non_zeros,3);
        Eigen::MatrixXf diff_xx = Eigen::MatrixXf::Zero(num_non_zeros_xx,3);
        Eigen::MatrixXf diff_yy = Eigen::MatrixXf::Zero(num_non_zeros_yy,3);
        Eigen::MatrixXf sum_diff_yx_2 = Eigen::MatrixXf::Zero(num_non_zeros,1);
        Eigen::MatrixXf sum_diff_xx_2 = Eigen::MatrixXf::Zero(num_non_zeros_xx,1);
        Eigen::MatrixXf sum_diff_yy_2 = Eigen::MatrixXf::Zero(num_non_zeros_yy,1);
        Eigen::Matrix<double, 1, 3> partial_omega;
        Eigen::Matrix<double, 1, 3> partial_v;
        double partial_dl = 0;

        int j = 0;
        // loop through non-zero ids in ith row of A
        for(Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(A,i); it; ++it){
            int idx = it.col();
            Ai(0,j) = it.value();    // extract current value in A
            cross_xy.row(j) = ((*cloud_x)[i].transpose().cross((*cloud_y)[idx].transpose()));
            diff_yx.row(j) = ((*cloud_y)[idx]-(*cloud_x)[i]).transpose();
            sum_diff_yx_2(j,0) = diff_yx.row(j).squaredNorm();
            ++j;
        }
        j = 0; 
        for(Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(Axx,i); it; ++it){
            int idx = it.col();
            Axxi(0,j) = it.value();    // extract current value in A
            diff_xx.row(j) = ((*cloud_x)[idx]-(*cloud_x)[i]).transpose();
            sum_diff_xx_2(j,0) = diff_xx.row(j).squaredNorm();
            ++j;
        }
        if(i<num_moving){
            j = 0; 
            for(Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(Ayy,i); it; ++it){
                int idx = it.col();
                Ayyi(0,j) = it.value();    // extract current value in A
                diff_yy.row(j) = ((*cloud_y)[idx]-(*cloud_y)[i]).transpose();
                ++j;
            }
            // update dl from Ayy
            partial_dl += double((1/ell_3*Ayyi*sum_diff_yy_2)(0,0));
        }
        partial_omega = (1/c*Ai*cross_xy).cast<double>();
        partial_v = (1/d*Ai*diff_yx).cast<double>();

        // update dl from Axy
        partial_dl -= double(2*(1/ell_3*Ai*sum_diff_yx_2)(0,0));
        
        // update dl from Axx
        partial_dl += double((1/ell_3*Axxi*sum_diff_xx_2)(0,0));

        // sum them up
        omegav_lock.lock();
            double_omega += partial_omega.transpose();
            double_v += partial_v.transpose();
            dl += partial_dl;
        omegav_lock.unlock();
    });
    // }

    // if num_moving > num_fixed, update the rest of Ayy to dl 
    if(num_moving>num_fixed){
        tbb::parallel_for(int(num_fixed),num_moving,[&](int i){
            int num_non_zeros_yy = Ayy.innerVector(i).nonZeros();
            Eigen::MatrixXf Ayyi = Eigen::MatrixXf::Zero(1,num_non_zeros_yy);
            Eigen::MatrixXf diff_yy = Eigen::MatrixXf::Zero(num_non_zeros_yy,3);
            Eigen::MatrixXf sum_diff_yy_2 = Eigen::MatrixXf::Zero(num_non_zeros_yy,1);
            double partial_dl = 0;

            int j = 0; 
            for(Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(Ayy,i); it; ++it){
                int idx = it.col();
                Ayyi(0,j) = it.value();    // extract current value in A
                diff_yy.row(j) = ((*cloud_y)[idx]-(*cloud_y)[i]).transpose();
                sum_diff_yy_2(j,0) = diff_yy.row(j).squaredNorm();
                ++j;
            }
            partial_dl += double((1/ell_3*Ayyi*sum_diff_yy_2)(0,0));

            omegav_lock.lock();
                dl += partial_dl;
            omegav_lock.unlock();
        });
    }
    

    // update them to class-wide variables
    omega = double_omega.cast<float>();
    v = double_v.cast<float>();
    dl = dl/(Axx.nonZeros()+Ayy.nonZeros()-2*A.nonZeros());
}


void acvo::compute_step_size(){
    // compute skew matrix
    Eigen::Matrix3f omega_hat = skew(omega);
    
    // compute xi*z+v, xi^2*z+xi*v, xi^3*z+xi^2*v, xi^4*z+xi^3*v
    Eigen::MatrixXf xiz(num_moving,3);
    Eigen::MatrixXf xi2z(num_moving,3);
    Eigen::MatrixXf xi3z(num_moving,3);
    Eigen::MatrixXf xi4z(num_moving,3);
    Eigen::MatrixXf normxiz2(num_moving,1);
    Eigen::MatrixXf xiz_dot_xi2z(num_moving,1);
    Eigen::MatrixXf epsil_const(num_moving,1);

    tbb::parallel_for( int(0), num_moving, [&]( int j ){
        Eigen::Vector3f cloud_yi = (*cloud_y)[j];
        xiz.row(j) = omega.transpose().cross(cloud_yi.transpose())+v.transpose(); // (xi*z+v)
        xi2z.row(j) = (omega_hat*omega_hat*cloud_yi\
                            +(omega_hat*v)).transpose();    // (xi^2*z+xi*v)
        xi3z.row(j) = (omega_hat*omega_hat*omega_hat*cloud_yi\
                            +(omega_hat*omega_hat*v)).transpose();  // (xi^3*z+xi^2*v)
        xi4z.row(j) = (omega_hat*omega_hat*omega_hat*omega_hat*cloud_yi\
                            +(omega_hat*omega_hat*omega_hat*v)).transpose();    // (xi^4*z+xi^3*v)
        normxiz2(j,0) = xiz.row(j).squaredNorm();
        xiz_dot_xi2z(j,0) = (-xiz.row(j).dot(xi2z.row(j)));
        epsil_const(j,0) = xi2z.row(j).squaredNorm()+2*xiz.row(j).dot(xi3z.row(j));
    });

    // initialize coefficients
    float temp_coef = 1/(2.0*ell*ell);   // 1/(2*l^2)
    double B = 0;
    double C = 0;
    double D = 0;
    double E = 0;

    tbb::spin_mutex BCDE_lock;
    
    tbb::parallel_for(int(0),num_fixed,[&](int i){
        // loop through used index in ith row
        double Bi=0;
        double Ci=0;
        double Di=0;
        double Ei=0;

        for(Eigen::SparseMatrix<float,Eigen::RowMajor>::InnerIterator it(A,i); it; ++it){
            int idx = it.col();
            
            // diff_xy = x[i] - y[used_idx[j]]
            auto diff_xy = ((*cloud_x)[i] - (*cloud_y)[idx]);    
            // beta_i = -1/l^2 * dot(xiz,diff_xy)
            float beta_ij = (-2.0*temp_coef * xiz.row(idx)*diff_xy).value();
            // gamma_i = -1/(2*l^2) * (norm(xiz).^2 + 2*dot(xi2z,diff_xy))
            float gamma_ij = (-temp_coef * (normxiz2.row(idx)\
                            + 2.0*xi2z.row(idx)*diff_xy)).value();
            // delta_i = 1/l^2 * (dot(-xiz,xi2z) + dot(-xi3z,diff_xy))
            float delta_ij = (2.0*temp_coef * (xiz_dot_xi2z.row(idx)\
                            + (-xi3z.row(idx)*diff_xy))).value();
            // epsil_i = -1/(2*l^2) * (norm(xi2z).^2 + 2*dot(xiz,xi3z) + 2*dot(xi4z,diff_xy))
            float epsil_ij = (-temp_coef * (epsil_const.row(idx)\
                            + 2.0*xi4z.row(idx)*diff_xy)).value();

            float A_ij = it.value();

            Bi += double(A_ij * beta_ij);
            Ci += double(A_ij * (gamma_ij+beta_ij*beta_ij/2.0));
            Di += double(A_ij * (delta_ij+beta_ij*gamma_ij + beta_ij*beta_ij*beta_ij/6.0));
            Ei += double(A_ij * (epsil_ij+beta_ij*delta_ij+1/2.0*beta_ij*beta_ij*gamma_ij\
                        + 1/2.0*gamma_ij*gamma_ij + 1/24.0*beta_ij*beta_ij*beta_ij*beta_ij));
        }
        
        // sum them up
        BCDE_lock.lock();
            B+=Bi;
            C+=Ci;
            D+=Di;
            E+=Ei;
        BCDE_lock.unlock();
    });

    Eigen::VectorXf p_coef(4);
    p_coef << 4.0*float(E),3.0*float(D),2.0*float(C),float(B);
    
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

void acvo::transform_pcd(){
    tbb::parallel_for(int(0), num_moving, [&]( int j ){
        (*cloud_y)[j] = transform.linear()*ptr_moving_pcd->positions[j]+transform.translation();
    });
    
}

void acvo::write_pcl_point_cloud_to_disk(point_cloud* ptr_pcd, string pcd_pth, string pcd_dso_pth){
    pcl::io::savePCDFileASCII (pcd_pth, ptr_pcd->pcl_cloud);
    pcl::io::savePCDFileASCII (pcd_dso_pth, ptr_pcd->pcl_dso_cloud);

}

float acvo::function_inner_product(point_cloud* cloud_a, point_cloud* cloud_b){

    double sum_A=0;
    double sum = 0;

    // convert k threshold to d2 threshold (so that we only need to calculate k when needed)
    float d2_thres = -2.0*ell*ell*log(sp_thres/sigma/sigma);
    float d2_c_thres = -2.0*c_ell*c_ell*log(sp_thres/c_sigma/c_sigma);

    /** 
     * kdtreeeeeeeeeeeeeeeeeeeee
     **/
    typedef KDTreeVectorOfVectorsAdaptor<cloud_t, float>  kd_tree_t;
	kd_tree_t mat_index(3 /*dim*/, (cloud_b->positions), 10 /* max leaf */ );
	mat_index.index->buildIndex();

    tbb::spin_mutex sum_A_lock;
    // loop through points
    tbb::parallel_for(int(0),cloud_a->num_points,[&](int i){

        const float search_radius = d2_thres;
		std::vector<std::pair<size_t,float>>  ret_matches;
		nanoflann::SearchParams params;
		//params.sorted = false;
		const size_t nMatches = mat_index.index->radiusSearch(&(cloud_a->positions)[i](0), search_radius, ret_matches, params);
        Eigen::Matrix<float,NUM_FEATURES,1> feature_a = cloud_a->features.row(i).transpose();

        for(size_t j=0; j<nMatches; ++j){
            int idx = ret_matches[j].first;
            float d2 = ret_matches[j].second;
            float k = 0;
            float ck = 0;
            float d2_color = 0;
            float a = 0;
            if(d2<d2_thres){
                Eigen::Matrix<float,NUM_FEATURES,1> feature_b = cloud_b->features.row(idx).transpose();
                d2_color = ((feature_a-feature_b).squaredNorm());

                if(d2_color<d2_c_thres){
                    k = sigma*sigma*exp(-d2/(2.0*ell*ell));
                    ck = c_sigma*c_sigma*exp(-d2_color/(2.0*c_ell*c_ell));
                    a = ck*k;
                    if (a > sp_thres){
                        sum_A_lock.lock();
                            sum_A+=a;
                            sum += 1;
                        sum_A_lock.unlock();
                    }
                }
            }
        }
    });

    return float(sum_A/sum);
}

void acvo::set_pcd(const int dataset_seq ,const cv::Mat& RGB_img,const cv::Mat& dep_img, string pcd_pth, string pcd_dso_pth){

    // create pcd_generator class
    pcd_generator pcd_gen;
    pcd_gen.dataset_seq = dataset_seq;
    
    // if it's the first image
    if(init == false){
        std::cout<<"initializing acvo..."<<std::endl;
        pcd_gen.load_image(RGB_img,dep_img,ptr_fixed_fr.get());
        pcd_gen.create_pointcloud(0,ptr_fixed_fr.get(), ptr_fixed_pcd.get());
        std::cout<<"first pcd generated!"<<std::endl;
        init = true;
        // write_pcl_point_cloud_to_disk(ptr_fixed_pcd.get(),pcd_pth,pcd_dso_pth);
        return;
    }

    ptr_moving_fr.reset(new frame);
    ptr_moving_pcd.reset(new point_cloud);

    pcd_gen.load_image(RGB_img,dep_img,ptr_moving_fr.get());
    pcd_gen.create_pointcloud(0,ptr_moving_fr.get(), ptr_moving_pcd.get());
    // write_pcl_point_cloud_to_disk(ptr_moving_pcd.get(),pcd_pth,pcd_dso_pth);

    // get total number of points
    num_fixed = ptr_fixed_pcd->num_points;
    num_moving = ptr_moving_pcd->num_points;
    std::cout<<"num fixed: "<<num_fixed<<std::endl;
    std::cout<<"num moving: "<<num_moving<<std::endl;

    // extract cloud x and y
    cloud_x = &(ptr_fixed_pcd->positions);
    cloud_y = new std::vector<Eigen::Vector3f>(ptr_moving_pcd->positions);

    // initialization of parameters
    ell = ell_init;
    ell_max = 0.15;
    dl = 0;
    A_trip_concur.reserve(num_moving*20);
    A.resize(num_fixed,num_moving);
    Axx.resize(num_fixed,num_fixed);
    Ayy.resize(num_moving,num_moving);
    A.setZero();
    Axx.setZero();
    Ayy.setZero();
}



void acvo::align(){
    int n = tbb::task_scheduler_init::default_num_threads();
    // std::cout<<"num_thread: "<<n<<std::endl;

    // loop until MAX_ITER
    for(int k=0; k<MAX_ITER; k++){
        // update transformation matrix
        update_tf();

        // apply transform to the point cloud
        transform_pcd();
        
        // compute omega and v
        compute_flow();

        // compute step size for integrating the flow
        compute_step_size();

        // stop if the step size is too small
        if(omega.cast<double>().norm()<eps && v.cast<double>().norm()<eps){
            iter = k;
            std::cout<<"norm, omega: "<<omega.norm()<<", v: "<<v.norm()<<std::endl;
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
            std::cout<<"dist: "<<dist_se3(dR,dT)<<std::endl;
            break;
        }

        // update ell
        ell = ell + dl_step*dl;
        
        // set search range for ell
        if(ell>=ell_max){
            ell = ell_max*0.7;
            ell_max = ell_max*0.7;
        }
        ell = (ell<ell_min)? ell_min:ell;
        
    }
    prev_transform = transform.matrix();
    accum_transform = accum_transform * transform.matrix();
    update_tf();

    ptr_fixed_pcd = std::move(ptr_moving_pcd);

    delete cloud_y;
}

void acvo::run_cvo(const int dataset_seq,const cv::Mat& RGB_img,const cv::Mat& dep_img, string pcd_pth, string pcd_dso_pth){

    if(init == false){
        set_pcd(dataset_seq,RGB_img, dep_img, pcd_pth, pcd_dso_pth);
    }
    else{
        set_pcd(dataset_seq, RGB_img, dep_img, pcd_pth, pcd_dso_pth);
        align();
        
        std::cout<<"Total iterations: "<<iter<<std::endl;
        std::cout<<"RKHS-SE(3) Object Transformation Estimate: \n"<<transform.matrix()<<std::endl;
    }

}
}