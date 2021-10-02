#include <iostream>
#include "utils/MathUtil.h"
#include "utils/GLUtil.h"
#include "render.h"

GLFWwindow *gWindow;

int main()
{
    auto render = new cRender();
    render->Init();
    gWindow = render->GetWindow();
    // set up cameras
    while (!glfwWindowShouldClose(gWindow))
    {
        glfwPollEvents();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        render->Update();
        // glBindVertexArray(0); // no need to unbind it every time

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(gWindow);
        glfwPollEvents();
    }
    delete render;
    glfwTerminate();
}