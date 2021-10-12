#include "RenderImGui.h"
#include "imgui.h"
#include "backends\imgui_impl_glfw.h"
#include "backends\imgui_impl_opengl3.h"
#include "GLFW/glfw3.h"
#include "GLFW/glfw3native.h"

void cRenderImGui::MouseMoveCallback(double xpos, double ypos)
{
    ImGuiIO &io = ImGui::GetIO();
    if (io.WantCaptureMouse)
        return;
    cRender::MouseMoveCallback(xpos, ypos);
}
void cRenderImGui::MouseButtonCallback(int but, int action, int mods)
{
    ImGuiIO &io = ImGui::GetIO();
    if (io.WantCaptureMouse)
        return;
    cRender::MouseButtonCallback(but, action, mods);
}
void cRenderImGui::KeyCallback(int key, int scancode, int action, int mods)
{
    ImGuiIO &io = ImGui::GetIO();
    if (io.WantCaptureKeyboard)
        return;
    cRender::KeyCallback(key, scancode, action, mods);
}
void cRenderImGui::ResizeCallback(int w, int h)
{
    cRender::ResizeCallback(w, h);
}
void cRenderImGui::ScrollCallback(double xoff, double yoff)
{
    ImGuiIO &io = ImGui::GetIO();
    if (io.WantCaptureMouse)
        return;
    cRender::ScrollCallback(xoff, yoff);
}

void cRenderImGui::InitGL()
{
    cRender::InitGL();
    // ------------------ add imgui code -------------------
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO &io = ImGui::GetIO();
    io.IniFilename = "";
    // (void)io;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsClassic();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(mWindow, true);
    const char *glsl_version = "#version 130";
    ImGui_ImplOpenGL3_Init(glsl_version);
}

void cRenderImGui::Update()
{
    cRender::Update();
    {

        // bool show_demo_window = true;
        // bool show_another_window = false;
        // ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
        {
            {
                ImVec2 init_window_size = ImVec2(300, 100);
                ImGui::SetNextWindowSize(init_window_size, ImGuiCond_Always);
                // if (mNeedToUpdateImGuiWindowPos == true)
                // {
                //     ImGui::SetNextWindowPos(ImVec2(float(gWindowWidth) - init_window_size.x, 0),
                //                             ImGuiCond_Always);
                //     mNeedToUpdateImGuiWindowPos = false;
                // }

                ImGuiWindowFlags window_flags = 0;
                // window_flags |= ImGuiWindowFlags_NoMove;
                // window_flags |= ImGuiWindowFlags_NoResize;
                bool open = false;
                bool *p_open = &open;

                ImGui::Begin("render setting", p_open, window_flags);
                ImGui::Text("hello imgui");
                ImGui::End();
            }
        }

        // Rendering
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    }
}
