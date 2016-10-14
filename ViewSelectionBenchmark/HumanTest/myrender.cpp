#include "myrender.h"

MyRender::MyRender()
{
    // init mesh to Null
    mesh = 0;
}

MyRender::~MyRender()
{
    delete mesh;
}

void MyRender::initializeGL()
{
    //init GLEW lib
    //glewInit();

    static const GLfloat light0_color[4] = { 0.8f, 0.8f, 0.8f, 1.0f };
    static const GLfloat light1_color[4] = { 0.4f, 0.4f, 0.8f, 1.0f };
    static const GLfloat light0_pos[4]   = { 100.0f, 100.0f, 100.0f, 0.0f };
    static const GLfloat light1_pos[4]   = { -100.0f, 100.0f, 100.0f, 0.0f };

    glClearColor(1,1,1,1);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    // speedup
    glEnable(GL_DITHER);
    glShadeModel(GL_SMOOTH);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_FASTEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_FASTEST);

    // light
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  light0_color);
    glLightfv(GL_LIGHT1, GL_DIFFUSE,  light1_color);
    glLightfv(GL_LIGHT0, GL_POSITION, light0_pos);
    glLightfv(GL_LIGHT1, GL_POSITION, light1_pos);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);

    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    // init texture type
    //glEnable(GL_TEXTURE_RECTANGLE_ARB);
    glEnable(GL_TEXTURE_2D);
}

void MyRender::resizeGL(int width, int height)
{
    glViewport(0,0,width,height);
}

void MyRender::paintGL()
{
    glClearColor(1,1,1,1);
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    if (mesh == 0 || mesh->vertices.size() == 0)
    {
        return;
    }

    float trans1[16],rot1[16];  // translate and rotate matrix
    buildTranslateMatrix(
        -mesh->bsphere.center[0],
        -mesh->bsphere.center[1],
        -mesh->bsphere.center[2],
        trans1);
    sgi_trackball_space::build_row_major_rotmatrix(rot1,trackball_quat);
    //set viewpoint according to auto sampling parameter
    buildLookAtMatrix(
        // eye position
        eye_pos_to_bsphere_center[0],
        eye_pos_to_bsphere_center[1],
        eye_pos_to_bsphere_center[2],
        // look at center
        eye_pos_to_bsphere_center[0],
        eye_pos_to_bsphere_center[1],
        0,
        // upright direction
        0,1,0,
        view_matrix);
    //buildOrthoMatrix(
    //    -mesh->bsphere.r, mesh->bsphere.r,
    //    -mesh->bsphere.r, mesh->bsphere.r,
    //    2*mesh->bsphere.r, 4*mesh->bsphere.r,
    //    project_matrix);
    buildPerspectiveMatrix(
        45.0, 1.0,
        0.1*mesh->bsphere.r, 999*mesh->bsphere.r,
        project_matrix );
    multMatrix(model_matrix, rot1,trans1);
    multMatrix(model_view_matrix, view_matrix, model_matrix);
    multMatrix(model_view_project_matrix,project_matrix,model_view_matrix);

    // drawing the model and depth
    glMatrixMode(GL_PROJECTION);
    transposeMatrix(project_matrix);
    glLoadMatrixf(project_matrix);

    glMatrixMode(GL_MODELVIEW);
    transposeMatrix(model_view_matrix);
    glLoadMatrixf(model_view_matrix);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);

    //float mat_specular_o[] = { 0.6, 0.6, 0.6, 1.0 };
    //float low_shininess_o[] = { 2.0 };
    //glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular_o);
    //glMaterialfv(GL_FRONT, GL_SHININESS, low_shininess_o);

    static bool on_preload_first_model = true; // the first model presented, are cannot init the vertex array, since when init in load_model, the OpenGL contex has not been setup.
    if (on_preload_first_model)
    {
        init_mesh_vertex_array_pointer();
        on_preload_first_model = false;
    }
    glColor3f(0.3f,0.6f,0.15f);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glDrawElements(GL_TRIANGLES,(GLsizei)mesh->faces.size()*3,GL_UNSIGNED_INT,&mesh->faces[0]);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
}

bool MyRender::load_model(const wstring filename)
{
    Mesh *tm;    //temp mesh;
    tm = Mesh::read(filename);
    if(tm)
    {
        delete mesh;
        mesh = tm;

        mesh->compute_bsphere();
        mesh->compute_normals();
        init_mesh_vertex_array_pointer();

        //reset eye position to bounding sphere center, it's relative to (0,0,0) coords.
        eye_pos_to_bsphere_center = point(0,0,3*mesh->bsphere.r);
        sgi_trackball_space::trackball(trackball_quat,0,0,0,0); // and also trackball quat
        return true;
    }
    return false;
}


bool MyRender::init_mesh_vertex_array_pointer()
{
    if (!mesh) return false;

    if (mesh->vertices.size() != 0)
    {
        glVertexPointer(3,GL_FLOAT,0,&mesh->vertices[0]);

        //this used as texture coordinate to get principal curvature
        mesh->vertices_sequence.resize(mesh->vertices.size());
        for (vector<point>::size_type i=0; i<mesh->vertices.size(); i++)
        {
            mesh->vertices_sequence[i] = i;
        }
        glTexCoordPointer(1,GL_INT,0,&mesh->vertices_sequence[0]);
    }
    if (mesh->normals.size() != 0)
    {
        glNormalPointer(GL_FLOAT,0,&mesh->normals[0]);
    }
    return true;
}

void MyRender::add_trackball_quat(float x, float y, float xx, float yy)
{
    //add a new rotation operation to this trackball
    {
        float quat[4];
        sgi_trackball_space::trackball(quat,x,y,xx,yy);
        sgi_trackball_space::add_quats(quat, trackball_quat, trackball_quat);
    }
}

void MyRender::set_trackball_quat(float q[4])
{
    trackball_quat[0] = q[0];
    trackball_quat[1] = q[1];
    trackball_quat[2] = q[2];
    trackball_quat[3] = q[3];
}

void MyRender::get_trackball_quat(float q[4])
{
    q[0] = trackball_quat[0] ;
    q[1] = trackball_quat[1] ;
    q[2] = trackball_quat[2] ;
    q[3] = trackball_quat[3] ;
}

void MyRender::zoom(float x, float y, float xx, float yy)
{
    return; // we don't zoom in this app
    (void)x;
    (void)xx;
    float shift = 1+(yy-y)/2.0;
    eye_pos_to_bsphere_center[2] *= shift;
}

void MyRender::pan(float x, float y, float xx, float yy)
{
    return; // we don't pan in this app
    (void)x;
    (void)y;
    (void)xx;
    (void)yy;
    float shift_x = (xx-x)*mesh->bsphere.r;
    float shift_y = (yy-y)*mesh->bsphere.r;
    
    eye_pos_to_bsphere_center[0] -= shift_x;
    eye_pos_to_bsphere_center[1] -= shift_y;
}
