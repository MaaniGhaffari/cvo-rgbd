cmake_minimum_required(VERSION 3.5.1)
project(rkhs_se3)


# Check C++11 or C++0x support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX17)
   if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -march=native -mkl -vec -DMKL_ILP64 -I${MKLROOT}/include -g")
   else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -msse -msse2 -msse3 -msse4 -mavx2 -ftree-vectorize -g")
   endif()
   add_definitions(-DCOMPILEDWITHC17)
   message(STATUS "Using flag -std=c++17.")
elseif(COMPILER_SUPPORTS_CXX11)
   if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -mkl -vec -DMKL_ILP64 -I${MKLROOT}/include -g")
   else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -msse -msse2 -msse3 -msse4 -mavx2 -ftree-vectorize -g")
   endif()
   add_definitions(-DCOMPILEDWITHC11)
   message(STATUS "Using flag -std=c++11.")
elseif(COMPILER_SUPPORTS_CXX0X)
   if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -march=native -mkl -vec -DMKL_ILP64 -I${MKLROOT}/include -g")
   else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -msse -msse2 -msse3 -msse4 -mavx2 -ftree-vectorize -g")
   endif()
   add_definitions(-DCOMPILEDWITHC0X)
   message(STATUS "Using flag -std=c++0x.")
else()
   message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()


# find package
find_package(PCL 1.4 REQUIRED COMPONENTS common io filters visualization)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED filesystem program_options regex timer chrono)
#find_package(OpenMP)
#find_package(TBB REQUIRED)

#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

# include
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${PCL_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
)



link_directories(${PCL_LIBRARY_DIRS})

add_definitions(${PCL_DEFINITIONS})

# declaire a c++ library
add_library(rkhs_se3
  src/rkhs_se3.cpp
  src/LieGroup.cpp
  src/pcl_visualizer.cpp
)

# add executable
add_executable(testrkhs src/main.cpp)

# link libraries
target_link_libraries(testrkhs
  rkhs_se3
  ${PCL_LIBRARIES}
  ${EIGEN3_LIBS}
  ${Boost_GENERAL}
  ${Boost_LIBRARIES}
)


