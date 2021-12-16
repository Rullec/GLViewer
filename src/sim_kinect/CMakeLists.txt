set(
    sim_kinect_src
    SimKinect.cpp
)

find_package(OpenMP REQUIRED)
add_definitions(-fopenmp)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_OPENMP")
add_library(sim_kinect_lib ${sim_kinect_src})
target_link_libraries(sim_kinect_lib OpenMP::OpenMP_CXX)

pybind11_add_module(sim_kinect ./sim_kinect_noise_pybind_export.cpp)
target_link_libraries(sim_kinect PRIVATE ${OpenCV_LIBS} utils_lib sim_kinect_lib)

message("current source dir " ${CMAKE_SOURCE_DIR})

set_target_properties(sim_kinect PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR} LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})