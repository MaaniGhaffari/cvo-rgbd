cmake_minimum_required(VERSION 3.5.1)
project(cvo_proj)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_module")

# Check C++11 or C++0x support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX17)
      message(STATUS "Using C++ 17")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -march=native -O3")
   if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp -vec")
      message(STATUS "Using intel compiler.")
   else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -msse -msse2 -msse3 -msse4 -mavx2 -ftree-vectorize")
      message(STATUS "Using other compiler.")
   endif()
   add_definitions(-DCOMPILEDWITHC17)
   message(STATUS "Using flag -std=c++17.")
elseif(COMPILER_SUPPORTS_CXX11)
   if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -mkl -vec -DMKL_ILP64 -I${MKLROOT}/include")
   else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -msse -msse2 -msse3 -msse4 -mavx2 -ftree-vectorize")
   endif()
   add_definitions(-DCOMPILEDWITHC11)
   message(STATUS "Using flag -std=c++11.")
elseif(COMPILER_SUPPORTS_CXX0X)
   if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -march=native -mkl -vec -DMKL_ILP64 -I${MKLROOT}/include")
   else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -msse -msse2 -msse3 -msse4 -mavx2 -ftree-vectorize")
   endif()
   add_definitions(-DCOMPILEDWITHC0X)
   message(STATUS "Using flag -std=c++0x.")
else()
   message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()


# find package
find_package(PCL 1.4 REQUIRED COMPONENTS common io filters visualization)
find_package(OpenCV 3.0.0 REQUIRED )
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED timer chrono)
find_package(OpenMP)
find_package(TBB REQUIRED)
# find_package(nanoflann)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

# include
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty
  ${EIGEN3_INCLUDE_DIR}
  ${PCL_INCLUDE_DIRS}
)

# declaire a c++ library
add_library(cvo_proj SHARED
  src/adaptive_cvo.cpp
  src/cvo.cpp
  src/LieGroup.cpp
  src/pcd_generator.cpp
  thirdparty/PixelSelector2.cpp
)
set_target_properties(cvo_proj PROPERTIES LIBRARY_OUTPUT_NAME cvo_proj)

# add executable
add_executable(adaptive_cvo src/adaptive_cvo_main.cpp)


# link libraries
target_link_libraries(adaptive_cvo
  cvo_proj
  ${PCL_LIBRARIES}
  ${OpenCV_LIBS}
  ${EIGEN3_LIBS}
  ${Boost_GENERAL}
  ${Boost_LIBRARIES}
  tbb
)

add_executable(cvo src/cvo_main.cpp)
target_link_libraries(cvo
  cvo_proj
  ${PCL_LIBRARIES}
  ${OpenCV_LIBS}
  ${EIGEN3_LIBS}
  ${Boost_GENERAL}
  ${Boost_LIBRARIES}
  tbb
)


export(PACKAGE cvo_proj)
set_target_properties(cvo_proj PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
get_property(cvo_include_dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
get_property(cvo_library_dirs TARGET cvo_proj PROPERTY LIBRARY_OUTPUT_DIRECTORY)
get_property(cvo_libraries TARGET cvo_proj PROPERTY LIBRARY_OUTPUT_NAME)

message("cvo_include_dirs: " ${cvo_include_dirs})
message("cvo_library_dirs: " ${cvo_library_dirs})
message("cvo_libraries: " ${cvo_libraries})

# Configure config file for local build tree
configure_file(cvoConfig.cmake.in
  "${PROJECT_BINARY_DIR}/cvoConfig.cmake" @ONLY)

message("PROJECT_BINARY_DIR: " ${PROJECT_BINARY_DIR})

