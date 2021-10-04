#pragma once
#include "utils/GLUtil.h"
#include <iostream>
extern unsigned int gWindowWidth, gWindowHeight;
void MouseMoveEventCallback(GLFWwindow *window, double xpos, double ypos);

void MouseButtonEventCallback(GLFWwindow *window, int button, int action, int mods);
void ErrorCallback(int error, const char *description);
void KeyEventCallback(GLFWwindow *window, int key, int scancode, int action, int mods);
void ResizeEventCallback(GLFWwindow *window, int w, int h);
void ScrollEventCallback(GLFWwindow *window, double xoff, double yoff);
extern const char *vertexShaderSource;
extern const char *fragmentShaderSource;
