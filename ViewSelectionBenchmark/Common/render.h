#pragma  once

class Render
{
public:
    Render(){}
    virtual ~Render(){}

    //called by Qt GLWidget
    virtual void initializeGL()=0;
    virtual void resizeGL(int width, int height)=0;
    virtual void paintGL()=0;
    virtual void add_trackball_quat(float x, float y, float xx, float yy)=0;
};
