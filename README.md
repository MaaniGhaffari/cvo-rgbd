# CVO/Adaptive CVO
* <a href="https://arxiv.org/pdf/1904.02266.pdf" target="_blank">Continuous Direct Sparse Visual Odometry from RGB-D Images.</a>
* <a href="https://arxiv.org/pdf/1910.00713.pdf" target="_blank">Adaptive Continous Visual Odometry from RGB-D Images.</a>

## MATLAB Examples for CVO
For a toy example of registration using Kinect data in MATLAB run `matlab/run_toy_example.m`

## Dependency
* ubuntu 16.04
* C++ 11 or higher
* Eigen3
* OpenCV 3.0.0
* PCL 1.4 (For saving pcd files)
* Intel C++ compiler (For better speed performance)
* Intel TBB
* Boost (for timing only)

## Environment Setting
* Add the following two source command to your ```.bashrc```, detailed instructions can be found [here](https://software.intel.com/en-us/articles/setting-up-the-build-environment-for-using-intel-c-or-fortran-compilers):

```
source opt/intel/mkl/bin/mklvars.sh intel64
source opt/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64
``` 

## Data
* Download TUM RGBD dataset from [here](https://vision.in.tum.de/data/datasets/rgbd-dataset/download).
* Generate the association files using ```data/rgbd_dataset/rgbd_benchmark_tools/assoc.sh```.

## cpp Code
To compile the cpp code, type the command below:
``` 
cd cpp/rkhs_se3_registration
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
Then executable files ```cvo``` and ```adaptive_cvo``` will be generated in build.

To run cvo code: 
```
./cvo $path_to_data $tum_sequence_number(1 for fr1, 2 for fr2, 3 for fr3)
```
A txt file containing the trajectory, ```cvo_poses_qt.txt```, will be generated in your data folder.

To run adaptive cvo code:
```
./adaptive_cvo $path_to_data $tum_sequence_number(1 for fr1, 2 for fr2, 3 for fr3)
```
A txt file containing the trajectory, ```acvo_poses_qt.txt```, will be generated in your data folder.

An example of how to run the code is in the script folder.

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
* Tzu-Yuan Lin, William Clark, Ryan M. Eustice, Jessy W. Grizzle, Anthony Bloch, and Maani Ghaffari. "Adaptive Continuous Visual Odometry from RGB-D Images." arXiv preprint arXiv:1910.00713, 2019. https://arxiv.org/abs/1910.00713
```
@article{lin2019adaptive,
  title={Adaptive Continuous Visual Odometry from RGB-D Images},
  author={Lin, Tzu-Yuan and Clark, William and Eustice, Ryan M and Grizzle, Jessy W and Bloch, Anthony and Ghaffari, Maani},
  journal={arXiv preprint arXiv:1910.00713},
  year={2019}
}
```
