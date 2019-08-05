/* ----------------------------------------------------------------------------
 * Copyright 2019, Tzu-yuan Lin <tzuyuan@umich.edu>, Maani Ghaffari <maanigj@umich.edu>
 * All Rights Reserved
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

/**
 *  @file   pcd_generator.cpp
 *  @author Tzu-yuan Lin, Maani Ghaffari 
 *  @brief  Source file for point cloud generator
 *  @date   July 29, 2019
 **/

#include "pcd_generator.hpp"


// using namespace cv;

namespace cvo{

    pcd_generator::pcd_generator():
    num_want(3000)
    {

    }

    pcd_generator::~pcd_generator(){
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
                // abs_l[idx] = sqrt(dx*dx+dy*dy);    
            

            }

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
        
        int idx = 0;
        for(int y=0; y<h; ++y){
            for(int x=0; x<w; ++x){
                // if the point is selected
                if(map[y*w+x]!=0){
                    ++idx;
                }
            }
        }
    }

    void pcd_generator::visualize_selected_pixels(frame* ptr_fr){
        
        int h = ptr_fr->h;
        int w = ptr_fr->w;

        // visualize the selected pixels in image
        cv::Mat img_selected;
        ptr_fr->image.copyTo(img_selected);
        for(int y=0; y<h; ++y){
            for(int x=0; x<w; ++x){
                if(map[y*w+x]==0){
                    img_selected.at<cv::Vec3b>(cv::Point(x, y)).val[0] = 0;
                    img_selected.at<cv::Vec3b>(cv::Point(x, y)).val[1] = 0;
                    img_selected.at<cv::Vec3b>(cv::Point(x, y)).val[2] = 0;
                }
            }
        }
        // cv::imshow("original image", ptr_fr->image);     // visualize original image
        cv::imshow("selected image", img_selected);      // visualize selected pixels
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
        case 1:
            fx = 517.3;  // focal length x
            fy = 516.5;  // focal length y
            cx = 318.6;  // optical center x
            cy = 255.3;  // optical center y
            break;
        case 2:
            fx = 520.9;  // focal length x
            fy = 521.0;  // focal length y
            cx = 325.1;  // optical center x
            cy = 249.7;  // optical center y
            break;
        case 3:
            fx = 535.4;  // focal length x
            fy = 539.2;  // focal length y
            cx = 320.1;  // optical center x
            cy = 247.6;  // optical center y
            break;
        
        default:
            // default set to fr1
            fx = 517.3;  // focal length x
            fy = 516.5;  // focal length y
            cx = 318.6;  // optical center x
            cy = 255.3;  // optical center y
            break;
        }

        int h = ptr_fr->h;
        int w = ptr_fr->w;
        
        int idx = 0;
        Eigen::Vector3f temp_position;
        cv::Mat temp_cv_position;
        for(int y=0; y<h; ++y){
            for(int x=0; x<w; ++x){
                ushort dep = ptr_fr->depth.at<ushort>(cv::Point(x, y));
                // if the point is selected
                if(map[y*w+x]!=0 && dep!=0){
                    // construct depth
                    temp_position(2) = dep/scaling_factor;
                    // construct x and y
                    temp_position(0) = (x-cx) * temp_position(2) / fx;
                    temp_position(1) = (y-cy) * temp_position(2) / fy;
                    
                    // add point to pcd
                    ptr_pcd->positions.emplace_back(temp_position);
                    // eigen2cv(temp_position,temp_cv_position);
                    // temp_cv_position = temp_cv_position.t();
                    // ptr_pcd->cv_positions.push_back(temp_cv_position);
                    // ptr_pcd->cv_positions.row(j) = temp_cv_position.clone();
                    // calculate dot positions
                    // ptr_pcd->dot_positions.emplace_back(temp_position.norm());

                    ++idx;
                }
            }
        }
        

  
        // remove nan points
        // ptr_pcd->positions.conservativeResize(idx,3);
        num_selected = idx;
    }

    void pcd_generator::get_features(frame* ptr_fr, point_cloud* ptr_pcd){

        int h = ptr_fr->h;
        int w = ptr_fr->w;

        ptr_pcd->features = Eigen::MatrixXf::Zero(num_selected,5);
        int idx = 0;

        for(int y=0; y<h; ++y){
            for(int x=0; x<w; ++x){
                // if the point is selected
                if(map[y*w+x]!=0 && ptr_fr->depth.at<ushort>(cv::Point(x, y))!=0){
                    
                    // extract bgr value
                    ptr_pcd->features(idx,2) = ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y)).val[0]; // b 
                    ptr_pcd->features(idx,1) = ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y)).val[1]; // g   
                    ptr_pcd->features(idx,0) = ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y)).val[2]; // r

                    // extract gradient
                    ptr_pcd->features(idx,3) = ptr_fr->dI[y*w+x][1];
                    ptr_pcd->features(idx,4) = ptr_fr->dI[y*w+x][2];
                    
                    // ptr_pcd->features(idx,3) = 0.5f*(ptr_fr->image.at<cv::Vec3b>(cv::Point(x+1, y)).val[2]
                    //                                 -ptr_fr->image.at<cv::Vec3b>(cv::Point(x-1, y)).val[2]);
                    // ptr_pcd->features(idx,4) = 0.5f*(ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y+1)).val[2]
                    //                                 -ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y-1)).val[2]);
                    // ptr_pcd->features(idx,5) = 0.5f*(ptr_fr->image.at<cv::Vec3b>(cv::Point(x+1, y)).val[1]
                    //                                 -ptr_fr->image.at<cv::Vec3b>(cv::Point(x-1, y)).val[1]);
                    // ptr_pcd->features(idx,6) = 0.5f*(ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y+1)).val[1]
                    //                                 -ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y-1)).val[1]);
                    // ptr_pcd->features(idx,7) = 0.5f*(ptr_fr->image.at<cv::Vec3b>(cv::Point(x+1, y)).val[0]
                    //                                 -ptr_fr->image.at<cv::Vec3b>(cv::Point(x-1, y)).val[0]);
                    // ptr_pcd->features(idx,8) = 0.5f*(ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y+1)).val[0]
                    //                                 -ptr_fr->image.at<cv::Vec3b>(cv::Point(x, y-1)).val[0]);                                                                                                                        

                    ++idx;
                }
            }
        }
    }

    void pcd_generator::load_image(const string& RGB_pth,const string& dep_pth,frame* ptr_fr){

        // load images                            
        ptr_fr->image = cv::imread(RGB_pth);
        ptr_fr->depth = cv::imread(dep_pth,CV_LOAD_IMAGE_ANYDEPTH);

        cv::cvtColor(ptr_fr->image, ptr_fr->intensity, cv::COLOR_RGB2GRAY);
    
        ptr_fr->h = ptr_fr->image.rows;
        ptr_fr->w = ptr_fr->image.cols;
    }

    void pcd_generator::create_pointcloud(frame* ptr_fr, point_cloud* ptr_pcd){
        
        int h = ptr_fr->h;
        int w = ptr_fr->w;

        select_point(ptr_fr);

        // visualize_selected_pixels(ptr_fr);

        get_points_from_pixels(ptr_fr, ptr_pcd);

        get_features(ptr_fr, ptr_pcd);

        ptr_pcd->num_points = num_selected;
    }

    // void pcd_generator::write_pcd_to_disk(point_cloud* ptr_pcd, const string& folder){

    //     pcl::PointCloud<pcl::PointXYZRGB> cloud;
        
    //     cloud.width = num_selected;
    //     cloud.height = 1;
    //     cloud.is_dense = false;
    //     cloud.points.resize (cloud.width * cloud.height);


    //     for(int i=0; i<num_selected; ++i){

    //         cloud.points[i].x = ptr_pcd->positions(i,0);
    //         cloud.points[i].y = ptr_pcd->positions(i,1);
    //         cloud.points[i].z = ptr_pcd->positions(i,2);

    //         cloud.points[i].r = ptr_pcd->features(i,0);
    //         cloud.points[i].g = ptr_pcd->features(i,1);
    //         cloud.points[i].b = ptr_pcd->features(i,2);
    //     }

    //     pcl::io::savePCDFileASCII (folder, cloud);

    // }
}