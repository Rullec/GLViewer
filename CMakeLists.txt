cmake_minimum_required(VERSION 3.0.0)
project(KinectCap VERSION 0.1.0)

# custom cmake scripts
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(./cmake/configure_movecmd.cmake)
set (CMAKE_CXX_STANDARD 20)

# include path and lib path
include_directories(./includes)
include_directories(.)

# find OPENGL
find_package(OpenGL REQUIRED)
if(OPENGL_FOUND)
    include_directories(${OpenGL_INCLUDE_DIRS})
    link_directories(${OpenGL_LIBRARY_DIRS})
    add_definitions(${OpenGL_DEFINITIONS})
endif()

# find glfw3
find_package(GLFW3 REQUIRED)

# find static glew
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/includes) # for glew configuration
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib) # for glew configuration


set(GLEW_USE_STATIC_LIBS True)
# pkgconfig
# set(GLEW_VERBOSE True)
# find_package(GLEW REQUIRED)
# if(GLEW_FOUND)
#     link_directories(${GLEW_LIBRARY_DIRS})
#     include_directories(${GLEW_INCLUDE_DIRS})
# endif()
find_package(PkgConfig REQUIRED)
pkg_search_module(GLEW REQUIRED glew)
if(GLEW_FOUND)
	message("glEw3 lib: " "${GLEW_LIBRARIES}")
	message("glEw3 lib dir: " "${GLEW_LIBRARY_DIRS}")
	message("glEw3 include: " "${GLEW_INCLUDE_DIRS}")
	link_directories(${GLEW_LIBRARY_DIRS})
	include_directories(${GLEW_INCLUDE_DIRS})
endif()

# collect gl libs
set(GL_LIBS ${GLFW3_LIBRARY} ${GLEW_STATIC_LIBRARIES} ${OPENGL_LIBRARIES})

# add sources
add_subdirectory(./utils)
add_subdirectory(./cameras)

set(libs utils_lib)



add_executable(main ./main.cpp RenderCallback.cpp render.cpp)  
target_link_libraries(main ${GL_LIBS} cam_lib ${libs})