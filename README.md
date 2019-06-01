# cvo-rgbd
Continuous Direct Sparse Visual Odometry from RGB-D Images.

## MATLAB Examples
For a toy example of registration using Kinect data in MATLAB run `matlab/run_toy_example.m`

## Dependency
* ubuntu 16.04
* MATLAB
* PointCloudLibrary 1.4
* Eigen3
* Intel C++ compiler (For better speed performance)
* Intel MKL

## Environment Setting
* Add the following two source command to your ```.bashrc```, detailed instructions can be found [here](https://software.intel.com/en-us/articles/setting-up-the-build-environment-for-using-intel-c-or-fortran-compilers):

```
source opt/intel/mkl/bin/mklvars.sh intel64
source opt/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64
``` 

## Data
* We provide simple MATLAB code in ```util/generate_pointclouds.m``` to generate downsampled pointclouds from TUM RGBD dataset.

* Download TUM RGBD dataset from [here](https://vision.in.tum.de/data/datasets/rgbd-dataset/download) and unzip it into the corresponding folder in ```data/rgbd_dataset/```.
* Create ```pcd_full``` and ```pcd_ds``` folders for generating downsample pointclouds by MATLAB:
```
cd data/rgbd_dataset/<your_data>
mkdir pcd_full
mkdir pcd_ds
```
* Open ```util/generate_pointclouds.m``` with MATLAB.
* Modify ```line 87``` to your data path.
* Generate dense pcd files by indicating ```down_sample = false``` . (```line 76```)
* Generate downsampled pcd files with ```down_sample = true```.


## cpp Code
To compile the cpp code, type the command below:
``` 
cd dev/cpp/rkhs_se3_registration
mkdir build
cd build
```
If this is your first time compiling using intel compiler, set your cmake varaibles by the following command: ([learn more here](https://gitlab.kitware.com/cmake/community/wikis/FAQ#how-do-i-use-a-different-compiler))
```
cmake .. -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc
make
```
Next time, you only need to do the following:
```
cmake ..
make
```
Then an execute file ```testrkhs``` will be generated in build. To run it, simply type:
```
./testrkhs
```

## Citations
* Maani Ghaffari, William Clark, Anthony Bloch, Ryan M. Eustice, and Jessy W. Grizzle. "Continuous Direct Sparse Visual Odometry from RGB-D Images," in Proceedings of Robotics: Science and Systems, Freiburg, Germany, June 2019. https://arxiv.org/abs/1904.02266
```
@INPROCEEDINGS{MGhaffari-RSS-19, 
    AUTHOR    = {Maani Ghaffari AND William Clark AND Anthony Bloch AND Ryan M. Eustice AND Jessy W. Grizzle}, 
    TITLE     = {Continuous Direct Sparse Visual Odometry from RGB-D Images}, 
    BOOKTITLE = {Proceedings of Robotics: Science and Systems}, 
    YEAR      = {2019}, 
    ADDRESS   = {Freiburg, Germany}, 
    MONTH     = {June} 
} 
```
