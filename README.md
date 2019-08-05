# cvo-rgbd
Continuous Direct Sparse Visual Odometry from RGB-D Images.

## MATLAB Examples
For a toy example of registration using Kinect data in MATLAB run `matlab/run_toy_example.m`

## Dependency
* ubuntu 16.04
* C++ 11 or higher
* Eigen3
* OpenCV 3.0.0
* Intel C++ compiler (For better speed performance)
* Intel TBB
* OpenMP
* Boost (for timing only)

## Environment Setting
* Add the following two source command to your ```.bashrc```, detailed instructions can be found [here](https://software.intel.com/en-us/articles/setting-up-the-build-environment-for-using-intel-c-or-fortran-compilers):

```
source opt/intel/mkl/bin/mklvars.sh intel64
source opt/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64
``` 

## Data
* Download TUM RGBD dataset from [here](https://vision.in.tum.de/data/datasets/rgbd-dataset/download) and unzip it into the corresponding folder in ```data/rgbd_dataset/```.

## cpp Code
Modify ```line 11``` in ```main.cpp``` to your data folder.

To compile the cpp code, type the command below:
``` 
cd cpp/rkhs_se3_registration
mkdir build
cd build
```
If this is your first time compiling using intel compiler, type: 
```
source PATH_TO_YOUR_INTEL_COMPILER/linux/bin/compilervars.sh intel64
``` 

Set your cmake varaibles by the following command: ([learn more here](https://gitlab.kitware.com/cmake/community/wikis/FAQ#how-do-i-use-a-different-compiler))
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
