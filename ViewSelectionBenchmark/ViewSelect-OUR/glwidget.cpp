#include <QtGui>
#include <QtOpenGL>
#include "glwidget.h"

#include <iostream>

GLWidget::GLWidget(Render *r, QWidget *parent)
    : render(r),QGLWidget(parent)
{
    setMouseTracking(true);
    left_button_down = false;
    selecting_parts_by_mouse = false;
}

GLWidget::~GLWidget()
{
}

QSize GLWidget::sizeHint() const
{
    //return QSize(128,128);
    //return QSize(256,256);
    return QSize(512,512);
    //return QSize(1024,1024);
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
void GLWidget::mouseMove()
{
	QSize sz = this->size();
	render->add_trackball_quat(
		(2.0*256 - sz.width()) / sz.width(),
		(sz.height() - 2.0*256) / sz.height(),
		(2.0*257 - sz.width())    / sz.width(),
		(sz.height() - 2.0*256)    / sz.height());
	this->updateGL();
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    if (!selecting_parts_by_mouse)
    {
        // rotate the model
        static float lastx=0;
        static float lasty=0;
        if (event->buttons() == Qt::LeftButton)
        {
            QSize sz = this->size();
            render->add_trackball_quat(
                (2.0*lastx - sz.width()) / sz.width(),
                (sz.height() - 2.0*lasty) / sz.height(),
                (2.0*event->x() - sz.width())    / sz.width(),
                (sz.height() - 2.0*event->y())    / sz.height());
            this->updateGL();
			LOG(QString::number(sz.width()) +" " + QString::number(sz.height())+" " + QString::number(lastx) +" " + QString::number(lasty)+" " + QString::number(event->x())+" " + QString::number(event->y()));
			
        }
        lastx = event->x();
        lasty = event->y();
    }
    else
    {
        // set for mouse selection of different parts
        if (left_button_down && event->buttons() == Qt::LeftButton)
        {   // moving
            render->mouseDrag(startx,starty,event->x(),event->y()); // draw the drag rect
            this->updateGL();
        }
    }
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    if (selecting_parts_by_mouse && event->button() == Qt::LeftButton)
    {
        startx = event->x();
        starty = event->y();
        left_button_down = true;
    }
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    if (selecting_parts_by_mouse && event->button() == Qt::LeftButton)
    {
        left_button_down = false;
        render->mouseDrag(startx,starty,event->x(),event->y(),1); // end drag, start select
        this->updateGL();
		
    }
}