#include "glwidget.h"

GLWidget::GLWidget(Render *r, QWidget *parent)
    : render(r),QGLWidget(parent)
{
    setMouseTracking(true);
}

GLWidget::~GLWidget()
{
}

QSize GLWidget::sizeHint() const
{
    return QSize(512,512);
    //return QSize(1024,1024);
}

int GLWidget::heightForWidth(int w)
{
    return w;
}

void GLWidget::initializeGL()
{
    render->initializeGL();
}

void GLWidget::paintGL()
{
    render->paintGL();
}

void GLWidget::resizeGL(int width, int height)
{
    render->resizeGL(width,height);
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    static float lastx=0;
    static float lasty=0;
    if (event->buttons() == Qt::LeftButton)
    {
        // rotate
        QSize sz = this->size();
        render->add_trackball_quat(
            (2.0*lastx - sz.width()) / sz.width(),
            (sz.height() - 2.0*lasty) / sz.height(),
            (2.0*event->x() - sz.width())    / sz.width(),
            (sz.height() - 2.0*event->y())    / sz.height());
        this->updateGL();
    }
    if (event->buttons() == Qt::RightButton)
    {
        // zoom
        QSize sz = this->size();
        render->zoom(
            (2.0*lastx - sz.width()) / sz.width(),
            (sz.height() - 2.0*lasty) / sz.height(),
            (2.0*event->x() - sz.width())    / sz.width(),
            (sz.height() - 2.0*event->y())    / sz.height());
        this->updateGL();
    }
    if (event->buttons() == Qt::MidButton)
    {
        QSize sz = this->size();
        render->pan(
            (2.0*lastx - sz.width()) / sz.width(),
            (sz.height() - 2.0*lasty) / sz.height(),
            (2.0*event->x() - sz.width())    / sz.width(),
            (sz.height() - 2.0*event->y())    / sz.height());
        this->updateGL();
    }
    lastx = event->x();
    lasty = event->y();
}