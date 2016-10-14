#pragma once

#include <QtGui>
//#include <QtOpenGL>
#include <QGLWidget>
#include "render.h"

class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget(Render *r, QWidget *parent = 0);
    ~GLWidget();

    QSize GLWidget::sizeHint() const;

protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);
    void mouseMoveEvent(QMouseEvent *event);

    Render *render;

public slots:
};
