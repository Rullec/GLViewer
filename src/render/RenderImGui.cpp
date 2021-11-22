#include "RenderImGui.h"
#include "imgui.h"
#include "backends\imgui_impl_glfw.h"
#include "backends\imgui_impl_opengl3.h"
#include "GLFW/glfw3.h"
#include "GLFW/glfw3native.h"
#include "cameras/CameraBase.h"

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

        // ImGui::ShowDemoWindow();
        // 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
        {
            {
                ImVec2 init_window_size = ImVec2(400, 200);
                ImGui::SetNextWindowSize(init_window_size, ImGuiCond_FirstUseEver);
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
                // ImGui::PushStyleColor(ImGuiCol_FrameBg, ImVec4(1, 0.7, 0.6, 1.00));
                // ImGui::Text("hello imgui");
                if (true == ImGui::Button("Reset camera"))
                {
                    mCam->ResetPos();
                }
                ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1, 0.7, 0.6, 1.00));

                // bool go = false, go1 = false;
                // ImGui::Checkbox("check1", &go);
                for (int i = 0; i < mRenderResources.size(); i++)
                {
                    auto res = mRenderResources[i];
                    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(
                                                             res->mColor[0],
                                                             res->mColor[1],
                                                             res->mColor[2],
                                                             1.0));
                    // bool & enable = ;
                    ImGui::PopStyleColor();
                    ImGui::Checkbox(res->mName.c_str(), &(mEnableRenderResource[i]));
                    if (mEnableRenderResource[i])
                    {
                        ImGui::SameLine();
                        std::string name = "adjust pos" + std::to_string(i);
                        ImGui::Checkbox(name.c_str(), &(mEnableTransformAdjust[i]));
                        // std::cout << mEnableTransformAdjust[i] << std::endl;
                        if (mEnableTransformAdjust[i] == true)
                        {
                            // ImGui::SameLine();
                            float tran_scale = 1e2;
                            float pos_val[3] = {
                                res->GetPos()[0] * tran_scale,
                                res->GetPos()[1] * tran_scale,
                                res->GetPos()[2] * tran_scale};
                            float y_rot = res->GetRotAxisAngle()[1];
                            std::string pos_name = "adjust Y [" + std::to_string(i) + "] [cm]";
                            std::string ang_name = "adjust Y [" + std::to_string(i) + "] [rad]";
                            // ImGui::DragFloat3(pos_name.c_str(), pos_val, 0.1, -20, 20);
                            ImGui::DragFloat(pos_name.c_str(), pos_val + 1, 0.1, -20, 20);
                            ImGui::DragFloat(ang_name.c_str(), &y_rot, 0.01, -3.14, 3.14);

                            res->SetPos(tVector3f(pos_val[0], pos_val[1], pos_val[2]) / tran_scale);
                            res->SetRotAxisAngle(tVector3f(0, y_rot, 0));
                        }
                    }
                    // ImGui::SameLine();
                }
                ImGui::PopStyleColor();
                ImGui::End();
            }
        }

        // Rendering
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    }
}
