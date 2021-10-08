#include "RenderCallback.h"
#include "render.h"
extern cRender * gRender;
void MouseMoveEventCallback(GLFWwindow *window, double xpos, double ypos)
{
    gRender->MouseMoveCallback(xpos, ypos);
}

void MouseButtonEventCallback(GLFWwindow *window, int button, int action, int mods)
{
    gRender->MouseButtonCallback(button, action, mods);
}

void ErrorCallback(int error, const char *description)
{
}

void KeyEventCallback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);

}

void ResizeEventCallback(GLFWwindow *window, int w, int h)
{
    // gScene->Update();
    glfwSwapBuffers(window);
}

void ScrollEventCallback(GLFWwindow *window, double xoff, double yoff)
{
    // std::cout << "scroll: x y = " << xoff << " " << yoff << std::endl;
}

// const char *vertexShaderSource = "#version 330 core\n"
//                                  "layout (location = 0) in vec3 aPos;\n"
//                                  "void main()\n"
//                                  "{\n"
//                                  "   gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
//                                  "}\0";
// const char *fragmentShaderSource = "#version 330 core\n"
//                                    "out vec4 FragColor;\n"
//                                    "void main()\n"
//                                    "{\n"
//                                    "   FragColor = vec4(1.0f, 0.5f, 0.2f, 1.0f);\n"
//                                    "}\n\0";
