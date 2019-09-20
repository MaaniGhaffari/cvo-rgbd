#include "cvo.hpp"

void load_file_name(string assoc_pth, vector<string> &vstrRGBName, \
                    vector<string> &vstrRGBPth, vector<string> &vstrDepPth);

void load_img(cv::Mat& RGB_img, cv::Mat& dep_img, string& RGB_pth, string& dep_pth);

int main(int argc, char** argv){
    
    string folder = argv[1];
    int dataset_seq = std::stoi(argv[2]);

    // downsampled pcd from tum rgbd dataset
    string assoc_pth = folder + "assoc.txt";
    
    // create our registration class
    cvo::cvo cvo;

    // load associate file
    vector<string> vstrRGBName;     // vector for image names
    vector<string> vstrRGBPth;
    vector<string> vstrDepPth;
    load_file_name(assoc_pth,vstrRGBName,vstrRGBPth,vstrDepPth);
    int num_img = vstrRGBName.size();
    std::cout<<"num images: "<<num_img<<std::endl;
    // num_img = 100;

    // // export as quarternion
    ofstream fPoseQtTxt;
    fPoseQtTxt.open(folder+"cvo_poses_qt.txt");

    boost::timer::cpu_timer total_time;
    cv::Mat RGB_img;
    cv::Mat dep_img;
    // loop through all the images in associate files
    for(int i=0;i<num_img;i++){

        string RGB_pth = folder + vstrRGBPth[i];
        string dep_pth = folder + vstrDepPth[i];
        string pcd_save_pth = folder +"pcd/"+ vstrRGBName[i] +".pcd";
        string pcd_dso_save_pth = folder +"pcd_dso/"+ vstrRGBName[i] + ".pcd";
        
        if(cvo.init){
            std::cout<<"----------------------"<<std::endl;
            std::cout<<"Processing frame "<<i<<std::endl;
            std::cout<<"Aligning " + vstrRGBName[i-1] + " and " + vstrRGBName[i] <<std::endl;
        }

        boost::timer::cpu_timer timer;
        
        load_img(RGB_img, dep_img, RGB_pth, dep_pth);
        cvo.run_cvo(dataset_seq, RGB_img, dep_img, pcd_save_pth, pcd_dso_save_pth);
            
        std::cout<<"elapse time: "<<timer.format()<<std::endl;
        

        // log out quaternion
        if(cvo.init){

        // log out transformation matrix for each frame
            Eigen::Quaternionf q(cvo.accum_transform.matrix().block<3,3>(0,0));
            fPoseQtTxt<<vstrRGBName[i]<<" ";
            fPoseQtTxt<<cvo.accum_transform(0,3)<<" "<<cvo.accum_transform(1,3)<<" "<<cvo.accum_transform(2,3)<<" "; 
            fPoseQtTxt<<q.x()<<" "<<q.y()<<" "<<q.z()<<" "<<q.w()<<"\n";
        }
    }

    std::cout<<"======================="<<std::endl;
    std::cout<<"Total time for "<<num_img<<" frames is: "<<total_time.format()<<std::endl;
    std::cout<<"======================="<<std::endl;

    fPoseQtTxt.close();
}

void load_file_name(string assoc_pth, vector<string> &vstrRGBName, \
                    vector<string> &vstrRGBPth, vector<string> &vstrDepPth){
    std::ifstream fAssociation;
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
            string RGB_pth;
            ss >> RGB_pth;
            vstrRGBPth.push_back(RGB_pth);
            string dep;
            ss >> dep;
            string depPth;
            ss >> depPth;
            vstrDepPth.push_back(depPth);
        }
    }
    fAssociation.close();
}


void load_img(cv::Mat& RGB_img, cv::Mat& dep_img, string& RGB_pth, string& dep_pth){
    RGB_img = cv::imread(RGB_pth);
    dep_img = cv::imread(dep_pth,CV_LOAD_IMAGE_ANYDEPTH);
}