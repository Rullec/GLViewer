#pragma once
#include "render.h"

class cRenderImGui : public cRender
{
public:
    virtual void Update();
    virtual void MouseMoveCallback(double xpos, double ypos);
    virtual void MouseButtonCallback(int but, int action, int mods);
    virtual void KeyCallback(int key, int scancode, int action, int mods);
    virtual void ResizeCallback(int w, int h);
    virtual void ScrollCallback(double xoff, double yoff);

protected:
    virtual void InitGL();
};