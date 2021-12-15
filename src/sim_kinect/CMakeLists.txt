set(
    sim_kinect_src
    SimKinect.cpp
)

find_package(OpenMP REQUIRED)
add_definitions(-fopenmp)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_OPENMP")
add_library(sim_kinect_lib ${sim_kinect_src})
target_link_libraries(sim_kinect_lib OpenMP::OpenMP_CXX)