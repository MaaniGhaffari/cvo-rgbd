/* ----------------------------------------------------------------------------
 * Copyright 2019, Tzu-yuan Lin <tzuyuan@umich.edu>, Maani Ghaffari <maanigj@umich.edu>
 * All Rights Reserved
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

/**
 *  @file   pcd_generator.cpp
 *  @author Tzu-yuan Lin, Maani Ghaffari 
 *  @brief  Source file for point cloud generator
 *  @date   September 20, 2019
 **/

#include "pcd_generator.hpp"


// using namespace cv;

namespace cvo{

    pcd_generator::pcd_generator():
    num_want(3000),
    dep_thres(20000),
    cam_info{1000,616.368,616.745,319.935,243.639}
    // ds_viewer(new pcl::visualization::PCLVisualizer ("CVO Pointcloud Visualization"))
    {
    }

    pcd_generator::~pcd_generator(){
        delete map;
    }

    void pcd_generator::make_pyramid(frame* ptr_fr){
        
        /** 
         * In this function we reference to Direct Sparse Odometry (DSO) by Engel et al.
         * https://github.com/JakobEngel/dso/blob/master/src/FullSystem/HessianBlocks.cpp#L128
         **/

        // initialize array for storing the pyramid
        int wl = ptr_fr->w;     // wl is the width of image at level l
        int hl = ptr_fr->h;     // hl is the height of image at level l
        for(int i=0; i<PYR_LEVELS; ++i){
            ptr_fr->dI_pyr[i] = new Eigen::Vector3f[wl*hl];
            ptr_fr->abs_squared_grad[i] = new float[wl*hl];
            wl /= 2;    // reduce the size of the image by 2 for each lower level
            hl /= 2;
        }

        ptr_fr->dI = ptr_fr->dI_pyr[0]; // dI = dI_pyr. Note: dI is a pointer so now dI and dI_ptr[0] point to the same place

        // extract intensity and flatten it into dI = dI_pyr[0]
        int h = ptr_fr->h;
        int w = ptr_fr->w;
        int _stride = ptr_fr->intensity.step;
        uint8_t *inten_val = ptr_fr->intensity.data;
        for(int i=0; i<h; ++i){
            for(int j=0; j<w; ++j){
                ptr_fr->dI[i*w+j][0] = inten_val[i*_stride+j];
            }
        }
        
        // initialize w and l to first level
        wl = w;
        hl = h;
        // start making pyramid, loop through different levels
        for(int lvl=0; lvl<PYR_LEVELS; ++lvl){
            
            // create a pointer point to dI at current level
            Eigen::Vector3f* dI_l = ptr_fr->dI_pyr[lvl];
        

            // create a pointer point to abs_squared_grad at current level
            float* abs_l = ptr_fr->abs_squared_grad[lvl];

            // if it's not the finest level, downsample 
            if(lvl>0){
                // create pointer to previous level
                int prev_lvl = lvl-1;
                int prev_wl = wl*2;
                Eigen::Vector3f* prev_dI_l = ptr_fr->dI_pyr[prev_lvl];

                // downsampling
                for(int y=0; y<hl; ++y){
                    for(int x=0; x<wl; ++x){
                        dI_l[x+y*wl][0] = 0.25f*(prev_dI_l[2*x   + 2*y*prev_wl][0] + \
												 prev_dI_l[2*x+1 + 2*y*prev_wl][0] + \
												 prev_dI_l[2*x   + 2*y*prev_wl+prev_wl][0] + \
												 prev_dI_l[2*x+1 + 2*y*prev_wl+prev_wl][0]);
                    }
                }
            }

            // calculate gradient
            // we skip the first row&col and the last row&col
            for(int idx=wl; idx<wl*(hl-1); ++idx){
                
                float dx = 0.5f*(dI_l[idx+1][0] - dI_l[idx-1][0]);
			    float dy = 0.5f*(dI_l[idx+wl][0] - dI_l[idx-wl][0]);

                // if it's not finite, set to 0
                if(!std::isfinite(dx)) dx=0;
			    if(!std::isfinite(dy)) dy=0;
                
                dI_l[idx][1] = dx;
			    dI_l[idx][2] = dy;

                // save current absolute gradient value (dx^2+dy^2) into ptr_fr->abs_squared_grad[lvl]
                abs_l[idx] = dx*dx+dy*dy;  
                ptr_fr->avg_abs_squared_grad += abs_l[idx];

            }

            ptr_fr->avg_abs_squared_grad/=(wl*(hl-1));

            // update level 
            wl/=2;
            hl/=2;
        }
    }

    void pcd_generator::select_point(frame* ptr_fr){
        
        int w = ptr_fr->w;
        int h = ptr_fr->h;
        
        make_pyramid(ptr_fr);   // create image pyramid

        map = new float[w*h];   // initialize the map for point selection

        dso::PixelSelector pixel_selector(w, h);    // create point selection class
        num_selected = pixel_selector.makeMaps(ptr_fr, map, num_want);
        
        // if the number of selected are less than expected, add results from Canny
        if(num_selected<num_want/3){
            int num_to_add = num_want - num_selected;
            int step = (h*w)/num_to_add;
            int block_size = 8;
            bool point_got = false;
            cv::Mat edge;
            cv::blur( ptr_fr->intensity, edge, cv::Size(3,3) );
            cv::Canny(edge, edge, 0, 25, 3);
            // cv::imshow("test",edge);
            for(int y=0; y<h; y+=block_size){
                for(int x=0; x<w; x+=block_size){
                    point_got = false;
                    for(int j=0; j<block_size; ++j){
                        for(int i=0; i<block_size; ++i){
                            if(edge.at<uchar>(cv::Point(x+i, y+j))!=0){
                                if(map[(y+j)*w+x+i]==0){
                                    map[(y+j)*w+x+i] = 1;
                                    point_got = true;
                                    break;
                                }
                            }
                        }
                        if(point_got==true){
                            break;
                        }
                    }
                }
            }
        }
    }

    void pcd_generator::visualize_selected_pixels(frame* ptr_fr){
        
        int h = ptr_fr->h;
        int w = ptr_fr->w;

        // visualize the selected pixels in image
        cv::Mat img_selected(480,640, CV_8UC1);
        cv::Mat img_map;
        cv::Mat img_selected_color;
        ptr_fr->image.copyTo(img_map);
        double dep_max;
        double dep_min;
        cv::minMaxLoc(ptr_fr->depth, &dep_min, &dep_max);
        std::cout<<dep_max<<std::endl;
        for(int y=0; y<h; ++y){
            for(int x=0; x<w; ++x){
                uint16_t dep = ptr_fr->depth.at<uint16_t>(cv::Point(x, y));
                if(map[y*w+x]==0 || dep==0 ||  isnan(dep)){
                    img_selected.at<uchar>(cv::Point(x,y)) = 0.0;
                }
                else{ //chose
                    float dep_temp = dep/dep_max*255;
                    img_selected.at<uchar>(cv::Point(x,y)) = int(dep_temp);
                }
            }
        }

        cv::applyColorMap(img_selected, img_selected_color, cv::COLORMAP_JET);

        for(int y=0; y<h; ++y){
            for(int x=0; x<w; ++x){
                if(img_selected.at<uchar>(cv::Point(x,y))==0){
                    float intensity_avg = (img_map.at<cv::Vec3b>(cv::Point(x, y)).val[0] + img_map.at<cv::Vec3b>(cv::Point(x, y)).val[1] + img_map.at<cv::Vec3b>(cv::Point(x, y)).val[2])/3.0;
                    img_map.at<cv::Vec3b>(cv::Point(x, y)).val[0] = intensity_avg;
                    img_map.at<cv::Vec3b>(cv::Point(x, y)).val[1] = intensity_avg;
                    img_map.at<cv::Vec3b>(cv::Point(x, y)).val[2] = intensity_avg;
                }
                else{ // selected point
                    img_map.at<cv::Vec3b>(cv::Point(x, y)).val[0] = img_selected_color.at<cv::Vec3b>(cv::Point(x, y)).val[0];
                    img_map.at<cv::Vec3b>(cv::Point(x, y)).val[1] = img_selected_color.at<cv::Vec3b>(cv::Point(x, y)).val[1];
                    img_map.at<cv::Vec3b>(cv::Point(x, y)).val[2] = img_selected_color.at<cv::Vec3b>(cv::Point(x, y)).val[2];
                }
            }
        }
        for(int y=0; y<h; ++y){
            for(int x=0; x<w; ++x){
                if(img_selected.at<uchar>(cv::Point(x,y))==0){
                }
                else{ // selected point
                    for(int c_y=y-1; c_y<y+2; c_y++){
                        for(int c_x=x-1; c_x<x+2; c_x++){
                            img_map.at<cv::Vec3b>(cv::Point(c_x, c_y)).val[0] = img_selected_color.at<cv::Vec3b>(cv::Point(x, y)).val[0];
                            img_map.at<cv::Vec3b>(cv::Point(c_x, c_y)).val[1] = img_selected_color.at<cv::Vec3b>(cv::Point(x, y)).val[1];
                            img_map.at<cv::Vec3b>(cv::Point(c_x, c_y)).val[2] = img_selected_color.at<cv::Vec3b>(cv::Point(x, y)).val[2];
                        }
                    }
                }
            }
        }

        cv::imshow("original image", ptr_fr->image);     // visualize original image
        cv::imshow("depth image", ptr_fr->depth);
        cv::imshow("selected image", img_map);      // visualize selected pixels
        cv::waitKey(0);

    }

    void pcd_generator::get_points_from_pixels(frame* ptr_fr, point_cloud* ptr_pcd){
        
        float scaling_factor = 5000;    // scaling factor for depth data
        float fx;  // focal length x
        float fy;  // focal length y
        float cx;  // optical center x
        float cy;  // optical center y

        // set camera parameters
        switch (dataset_seq){
        // 0: real sense camera
        case 0:
            cam_info.scaling_factor = 1000.0;
            cam_info.fx = 616.368;  // focal length x
            cam_info.fy = 616.745;  // focal length y
            cam_info.cx = 319.935;  // optical center x
            cam_info.cy = 243.639;  // optical center y
            break;
        case 1: // fr1
            cam_info.scaling_factor = 5000.0;    // scaling factor for depth data
            cam_info.fx = 517.3;  // focal length x
            cam_info.fy = 516.5;  // focal length y
            cam_info.cx = 318.6;  // optical center x
            cam_info.cy = 255.3;  // optical center y
            break;
        case 2: // fr2
            cam_info.scaling_factor = 5000.0;
            cam_info.fx = 520.9;  // focal length x
            cam_info.fy = 521.0;  // focal length y
            cam_info.cx = 325.1;  // optical center x
            cam_info.cy = 249.7;  // optical center y
            break;
        case 3: // fr3
            cam_info.scaling_factor = 5000.0;
            cam_info.fx = 535.4;  // focal length x
            cam_info.fy = 539.2;  // focal length y
            cam_info.cx = 320.1;  // optical center x
            cam_info.cy = 247.6;  // optical center y
            break;
        case 4: // kitti 15
            cam_info.scaling_factor = 2000.0;
            cam_info.fx = 718.856;
            cam_info.fy = 718.856;
            cam_info.cx = 607.1928;
            cam_info.cy = 185.2157;
            break;

        case 5: // kitti 05
            cam_info.scaling_factor = 2000.0;
            cam_info.fx = 707.0912;
            cam_info.fy = 707.0912;
            cam_info.cx = 601.8873;
            cam_info.cy = 183.1104;
            break;
        
        default:
            // default set to real sense
            cam_info.scaling_factor = 1000.0;
            cam_info.fx = 616.368;  // focal length x
            cam_info.fy = 616.745;  // focal length y
            cam_info.cx = 319.935;  // optical center x
            cam_info.cy = 243.639;  // optical center y
            break;
        }

        int h = ptr_fr->h;
        int w = ptr_fr->w;
        
        int idx = 0;
        Eigen::Vector3f temp_position;
        cv::Mat temp_cv_position;
        for(int y=0; y<h; ++y){
            for(int x=0; x<w; ++x){
                uint16_t dep = ptr_fr->depth.at<uint16_t>(cv::Point(x, y));
                // if the point is selected
                if(map[y*w+x]!=0 && dep!=0 && !isnan(dep)){
                    // construct depth
                    temp_position(2) = dep/cam_info.scaling_factor;
                    // construct x and y
                    temp_position(0) = (x-cam_info.cx) * temp_position(2) / cam_info.fx;
                    temp_position(1) = (y-cam_info.cy) * temp_position(2) / cam_info.fy;
                    
                    // add point to pcd
                    ptr_pcd->positions.push_back(temp_position);

                    ++idx;
                }
            }
        }
        num_selected = idx;

        
  
        
    }

    void pcd_generator::get_features(const int feature_type, frame* ptr_fr, point_cloud* ptr_pcd){

        int h = ptr_fr->h;
        int w = ptr_fr->w;
        ptr_pcd->features = Eigen::MatrixXf::Zero(num_selected,NUM_FEATURES);
        int idx = 0;

        if(feature_type == 0){      // HSV + gradient and normalized to 0~1
            for(int y=0; y<h; ++y){
                for(int x=0; x<w; ++x){
                    // if the point is selected
                    uint16_t dep = ptr_fr->depth.at<uint16_t>(cv::Point(x, y));
                    if(map[y*w+x]!=0 && dep!=0  && !isnan(dep)){
                        
                        // extract bgr value
                        ptr_pcd->features(idx,0) = ptr_fr->image_hsv.at<cv::Vec3b>(cv::Point(x, y)).val[0]/180.0; // b
                        ptr_pcd->features(idx,1) = ptr_fr->image_hsv.at<cv::Vec3b>(cv::Point(x, y)).val[1]/255.0; // g   
                        ptr_pcd->features(idx,2) = ptr_fr->image_hsv.at<cv::Vec3b>(cv::Point(x, y)).val[2]/255.0; // r  
                        
                        // ptr_pcd->features(idx,0) = ptr_fr->intensity.at<cv::Vec3b>(cv::Point(x, y)).val[0]/255.0;
                        
                        // extract gradient
                        ptr_pcd->features(idx,3) = ptr_fr->dI[y*w+x][1]/255.0*2;
                        ptr_pcd->features(idx,4) = ptr_fr->dI[y*w+x][2]/255.0*2;
                        
                        ++idx;
                    }
                }
            }
        }
        else if(feature_type == 1){     // RGB + gradient without normalize to 0~1       
            for(int y=0; y<h; ++y){
                for(int x=0; x<w; ++x){
                    // if the point is selected
                    uint16_t dep = ptr_fr->depth.at<uint16_t>(cv::Point(x, y));
                    if(map[y*w+x]!=0 && dep!=0  && !isnan(dep)){
                        
                        // extract bgr value
                        ptr_pcd->features(idx,0) = ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y)).val[0]; // b
                        ptr_pcd->features(idx,1) = ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y)).val[1]; // g   
                        ptr_pcd->features(idx,2) = ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y)).val[2]; // r  
                        
                        // ptr_pcd->features(idx,0) = ptr_fr->intensity.at<cv::Vec3b>(cv::Point(x, y)).val[0]/255.0;
                        
                        // extract gradient
                        ptr_pcd->features(idx,3) = ptr_fr->dI[y*w+x][1];
                        ptr_pcd->features(idx,4) = ptr_fr->dI[y*w+x][2];
                        
                        ++idx;
                    }
                }
            }
        }
    }

    void pcd_generator::load_image(const cv::Mat& RGB_img,const cv::Mat& dep_img,frame* ptr_fr){

        // load images                            
        ptr_fr->image = RGB_img;
        ptr_fr->depth = dep_img;

        cv::cvtColor(ptr_fr->image, ptr_fr->intensity, cv::COLOR_RGB2GRAY);
        cv::cvtColor(ptr_fr->image, ptr_fr->image_hsv, cv::COLOR_RGB2HSV);

        
        ptr_fr->h = ptr_fr->image.rows;
        ptr_fr->w = ptr_fr->image.cols;
    }

    void pcd_generator::create_pointcloud(const int feature_type, frame* ptr_fr, point_cloud* ptr_pcd){
        
        int h = ptr_fr->h;
        int w = ptr_fr->w;

        select_point(ptr_fr);

        // visualize_selected_pixels(ptr_fr);

        get_points_from_pixels(ptr_fr, ptr_pcd);

        get_features(feature_type, ptr_fr, ptr_pcd);

        // create_pcl_pointcloud(ptr_pcd);

        ptr_pcd->num_points = num_selected;

        // release memories
        for(int i=0;i<PYR_LEVELS;i++){
			delete[] ptr_fr->dI_pyr[i];
			delete[] ptr_fr->abs_squared_grad[i];
		}
    }

    void pcd_generator::create_pcl_pointcloud(frame* ptr_fr, point_cloud* ptr_pcd){

        /**
         *  for visualization
         **/

        int h = ptr_fr->h;
        int w = ptr_fr->w;
        // full point cloud
        ptr_pcd->pcl_cloud.width = h*w;
        ptr_pcd->pcl_cloud.height = 1;
        ptr_pcd->pcl_cloud.is_dense = false;
        ptr_pcd->pcl_cloud.points.resize (ptr_pcd->pcl_cloud.width * ptr_pcd->pcl_cloud.height);

        // for visualization
        ptr_pcd->pcl_dso_cloud.width = num_selected;
        ptr_pcd->pcl_dso_cloud.height = 1;
        ptr_pcd->pcl_dso_cloud.is_dense = false;
        ptr_pcd->pcl_dso_cloud.points.resize (ptr_pcd->pcl_dso_cloud.width * ptr_pcd->pcl_dso_cloud.height);
        int i = 0;
        int idx_vis = 0;
        for(int y=0; y<h; ++y){
            for(int x=0; x<w; ++x){
                uint16_t dep = ptr_fr->depth.at<uint16_t>(cv::Point(x, y));
                ptr_pcd->pcl_cloud.points[i].x = (x-cam_info.cx) * dep/cam_info.scaling_factor / cam_info.fx;
                ptr_pcd->pcl_cloud.points[i].y = (y-cam_info.cy) * dep/cam_info.scaling_factor / cam_info.fy;
                ptr_pcd->pcl_cloud.points[i].z = dep/cam_info.scaling_factor;

                ptr_pcd->pcl_cloud.points[i].r = ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y)).val[2]; // r
                ptr_pcd->pcl_cloud.points[i].g = ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y)).val[1]; // g  
                ptr_pcd->pcl_cloud.points[i].b = ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y)).val[0]; // b 

                if(map[y*w+x]!=0 && dep!=0 && !isnan(dep)){
                    // for visualization
                    ptr_pcd->pcl_dso_cloud.points[idx_vis].x = ptr_pcd->pcl_cloud.points[i].x;
                    ptr_pcd->pcl_dso_cloud.points[idx_vis].y = ptr_pcd->pcl_cloud.points[i].y;
                    ptr_pcd->pcl_dso_cloud.points[idx_vis].z = ptr_pcd->pcl_cloud.points[i].z;

                    ptr_pcd->pcl_dso_cloud.points[idx_vis].r = ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y)).val[2]; // r
                    ptr_pcd->pcl_dso_cloud.points[idx_vis].g = ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y)).val[1]; // g  
                    ptr_pcd->pcl_dso_cloud.points[idx_vis].b = ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y)).val[0]; // b 

                    ++idx_vis;
                }
                ++i;
            }
        }
    }
}