set(
    render_src
    RenderCallback.cpp
    render.cpp
    RenderResource.cpp
    RenderImGui.cpp
)

add_library(render_lib ${render_src})
target_link_libraries(render_lib utils_lib sim_kinect_lib cam_lib res_lib imgui_lib ${GL_LIBS}) 