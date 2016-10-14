#pragma once

#include <QObject>
#include <QtOpenGL>

//#include "gl/glew.h"
//#pragma comment(lib,"glew32")

#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <ctime>

#include "render.h"
#include "mesh.h"
#include "trackball.h"
#include "glmatrix.h"

class MyRender : public QObject, public Render
{
    Q_OBJECT
public:
    static const int glcanvas_size = 512;

public:
    MyRender(); 
    ~MyRender();

    // virtual functions inherited from Render,
    // referred by GLWidget
    void initializeGL();
    void resizeGL(int width, int height);
    void paintGL();
    // for trackball manipulation
    void add_trackball_quat(float x, float y, float xx, float yy);
    void set_trackball_quat(float q[4]);
    void get_trackball_quat(float q[4]);
    void zoom(float x, float y, float xx, float yy);
    void pan(float x, float y, float xx, float yy);

public slots:
    // slots for manipulation
    bool load_model(const wstring filename);

private:
    // initiate and data processing functions
    bool init_mesh_vertex_array_pointer();
    void build_project_matrix();

private:    // data of MyRender
    Mesh *mesh;

    // for transform, rotation
    float trackball_quat[4];
    point eye_pos_to_bsphere_center;

    // controlling matrix
    float model_matrix[16];
    float view_matrix[16];
    float model_view_matrix[16];
    float project_matrix[16];
    float model_view_project_matrix[16];
};

