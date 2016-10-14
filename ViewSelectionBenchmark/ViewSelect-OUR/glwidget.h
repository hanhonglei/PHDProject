#pragma once
#include "log.h"
#include <QGLWidget>
#include "render.h"

class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget(Render *r, QWidget *parent = 0);
    ~GLWidget();

    QSize sizeHint() const;
    bool selecting_parts_by_mouse;

    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);
    void mouseMoveEvent(QMouseEvent *event);
	void mouseMove();
    // for mouse drag event
    float startx;
    float starty;
    bool left_button_down;
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    Render *render;
};
