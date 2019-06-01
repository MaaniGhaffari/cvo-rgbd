#include "rkhs_se3.hpp"
#include "pcl_visualizer.hpp"

void load_file_name(string assoc_pth, vector<string> &vstrRGBName);

int main(int argc, char** argv){
    
    // string dataset = "freiburg3_structure_notexture_far";
    string dataset = "freiburg1_desk";

    // downsampled pcd from tum rgbd dataset
    string ds_folder = "../../../../data/rgbd_dataset/" + dataset + "/pcd_ds/";
    string dense_folder = "../../../../data/rgbd_dataset/" + dataset + "/pcd_full/";
    string assoc_pth = "../../../../data/rgbd_dataset/" + dataset + "/assoc.txt";

    // create our registration class
    rkhs_se3 cvo;
    pcl_visualizer visualizer;

    // load associate file
    vector<string> vstrRGBName;     // vector for image names
    load_file_name(assoc_pth,vstrRGBName);
    int num_img = vstrRGBName.size();
    std::cout<<"num images: "<<num_img<<std::endl;

    // ofstream fPoseCsv;
    // fPoseCsv.open("cvo_poses.csv");
    // fPoseCsv << "frame1, frame2, tx, ty, tz, r11, r12, r13, r21, r22, r23, r31, r32, r33" << endl;

    // loop through all the images in associate files
    // change this for loop to desired number of images for alignment
    for(int i=0;i<1;i++){
        string fixed_pth = ds_folder + vstrRGBName[i] + ".pcd";
        string moving_pth = ds_folder + vstrRGBName[i+1] + ".pcd";
        string fixed_vis_pth = dense_folder + vstrRGBName[i] + ".pcd";
        string moving_vis_pth = dense_folder + vstrRGBName[i+1] + ".pcd";

        std::cout<<"----------------------"<<std::endl;
        std::cout<<"Processing frame "<<i<<std::endl;
        std::cout<<"Aligning " + vstrRGBName[i] + " and " + vstrRGBName[i+1] <<std::endl;

        boost::timer::cpu_timer timer;
        
        cvo.set_pcd(fixed_pth,moving_pth);
        cvo.align();
        
        std::cout<<"elapse time: "<<timer.format()<<std::endl;
        std::cout<<"Total iterations: "<<cvo.iter<<std::endl;
        std::cout<<"RKHS-SE(3) Object Transformation Estimate: \n"<<cvo.transform.matrix()<<std::endl;
        
        // log out transformation matrix for each frame
        // fPoseCsv << i << "," << i+1 << "," \
        //         << cvo.transform(0,3) << "," << cvo.transform(1,3) << "," << cvo.transform(2,3) << ","\
        //         << cvo.transform(0,0) << "," << cvo.transform(0,1) << "," << cvo.transform(0,2) << ","\
        //         << cvo.transform(1,0) << "," << cvo.transform(1,1) << "," << cvo.transform(1,2) << ","\
        //         << cvo.transform(2,0) << "," << cvo.transform(2,1) << "," << cvo.transform(2,2) << std::endl;

        visualizer.add_pcd_to_viewer(fixed_vis_pth,moving_vis_pth,cvo.transform);
    }

    // fPoseCsv.close();

    // DO NOT enable visualize() and visualize_unaligned_pcd() at the same time
    // visualizer.visualize();
    visualizer.visualize_unaligned_pcd();
}

void load_file_name(string assoc_pth, vector<string> &vstrRGBName){
    ifstream fAssociation;
    fAssociation.open(assoc_pth.c_str());
    while(!fAssociation.eof())
    {
        string s;
        getline(fAssociation,s);
        if(!s.empty())
        {
            stringstream ss;
            ss << s;
            string RGB;
            ss >> RGB;
            vstrRGBName.push_back(RGB);
        }
    }
    fAssociation.close();
}

