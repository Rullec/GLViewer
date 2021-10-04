#include <iostream>
#include "utils/MathUtil.h"
#include "utils/GLUtil.h"
#include "render.h"

GLFWwindow *gWindow;
cRender * gRender;
int main()
{
    gRender = new cRender();
    gRender->Init();
    gWindow = gRender->GetWindow();
    // set up cameras
    while (!glfwWindowShouldClose(gWindow))
    {
        glfwPollEvents();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        gRender->Update();
        // glBindVertexArray(0); // no need to unbind it every time

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(gWindow);
        glfwPollEvents();
    }
    delete gRender;
    glfwTerminate();
}