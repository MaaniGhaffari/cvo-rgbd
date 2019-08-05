cmake_minimum_required(VERSION 3.5.1)
project(cvo)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_module")

# Check C++11 or C++0x support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX17)
      message(STATUS "Using C++ 17")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -march=native -qopenmp -O3 -vec")
   if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
      message(STATUS "Using intel compiler.")
   else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -msse -msse2 -msse3 -msse4 -mavx2 -ftree-vectorize")
      message(STATUS "Using other compiler.")
   endif()
   add_definitions(-DCOMPILEDWITHC17)
   message(STATUS "Using flag -std=c++17.")
elseif(COMPILER_SUPPORTS_CXX11)
   if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -qopenmp -O3 -vec")
   else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -msse -msse2 -msse3 -msse4 -mavx2 -ftree-vectorize")
   endif()
   add_definitions(-DCOMPILEDWITHC11)
   message(STATUS "Using flag -std=c++11.")
elseif(COMPILER_SUPPORTS_CXX0X)
   if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -march=native -qopenmp -O3 -vec")
   else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -msse -msse2 -msse3 -msse4 -mavx2 -ftree-vectorize")
   endif()
   add_definitions(-DCOMPILEDWITHC0X)
   message(STATUS "Using flag -std=c++0x.")
else()
   message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()


# find package
find_package(OpenCV 3.0.0 REQUIRED )
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED timer chrono)
find_package(OpenMP)
find_package(TBB REQUIRED)

# include
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty
  ${EIGEN3_INCLUDE_DIR}
)

# declaire a c++ library
add_library(cvo
  src/rkhs_se3.cpp
  src/LieGroup.cpp
  src/pcd_generator.cpp
  thirdparty/PixelSelector2.cpp
)

# add executable
add_executable(testrkhs src/main.cpp)

# link libraries
target_link_libraries(testrkhs
  cvo
  ${OpenCV_LIBS}
  ${EIGEN3_LIBS}
  ${Boost_GENERAL}
  ${Boost_LIBRARIES}
  tbb
)
