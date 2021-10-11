cmake_minimum_required(VERSION 3.0.0)
project(KinectCap VERSION 0.1.0)

# custom cmake scripts
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(./cmake/configure_movecmd.cmake)
set (CMAKE_CXX_STANDARD 20)

if(WIN32)
	include_directories(includes)
	link_directories(lib)
    add_definitions(-D _USE_MATH_DEFINES)
    add_definitions(-D _CRT_NO_VA_START_VALIDATION)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++17 /bigobj")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4061 /wd4244 /wd4246 /wd4305 /wd4267 /wd4711 /wd4710 /wd4514 /wd4477")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
endif()
# include path and lib path
include_directories(./includes)
include_directories(.)
include_directories(./src)

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
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})


if(UNIX)

find_package(PkgConfig REQUIRED)
pkg_search_module(GLEW REQUIRED glew)
elseif(WIN32)

find_package(GLEW REQUIRED glew)
endif()
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
add_subdirectory(./src/utils)
add_subdirectory(./src/cameras)
add_subdirectory(./src/geometries)
add_subdirectory(./src/restore)

set(libs utils_lib ${OpenCV_LIBS} geo_lib res_lib)


add_executable(main ./main.cpp src/RenderCallback.cpp src/render.cpp)
target_link_libraries(main ${GL_LIBS} cam_lib ${libs})


if(WIN32)
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        add_custom_command(
            TARGET main POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E rename
                ${CMAKE_CURRENT_BINARY_DIR}/Debug/main.exe
                ${CMAKE_SOURCE_DIR}/main.exe
        )
    else()
        add_custom_command(
            TARGET main POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E rename
                ${CMAKE_CURRENT_BINARY_DIR}/Release/main.exe
                ${CMAKE_SOURCE_DIR}/main.exe
        )
    endif()
else()
    add_custom_command(
        TARGET main POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E rename
            ${CMAKE_CURRENT_BINARY_DIR}/main        
            ${CMAKE_SOURCE_DIR}/main    
    )
endif()