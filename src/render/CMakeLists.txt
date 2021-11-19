set(
    render_src
    RenderCallback.cpp
    Render.cpp
    RenderResource.cpp
    RenderImGui.cpp
)

add_library(render_lib ${render_src})