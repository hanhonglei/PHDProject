#include <gl/glew.h>
#include <gl/freeglut.h>

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include <boost/lexical_cast.hpp>
#include <boost/smart_ptr.hpp>

#include <vec.h>
#include "../common/com.h"
#include "../common/mesh.h"
#include "../Common/trackball.h"
#include "../Common/errorhandle.h"
#include "../Common/glmatrix.h"

//#define cimg_use_magick
//#define MAGICK_STATIC
//#include "CImgGL.h"

// verbose mode of mouse and keyboard events
#define MOUSE_KEYBOARD_VERBOSE 0

//  Initialization
void init();

//  Callback functions
void display(void);
void reshape(int w, int h);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void pmotion(int x, int y);
void keyboard(unsigned char key, int x, int y);
void special(int key, int x, int y);

//  Support Functions
void center_window();

//  define the window position on screen
int window_x = 0;
int window_y = 0;

//  variables representing the window size
int window_width = 512;
int window_height = 512;

//  variable representing the window title
char *window_title = "Human Data Statistic";

//  Tells whether to display the window full screen or not
int full_screen = 0;

//  Viewing parameters
Mesh *mesh = 0;

// for transform, rotation
float trackball_quat[4] = {0,0,0,0};
point eye_pos_to_bsphere_center;

// controlling matrix
float model_matrix[16];
float view_matrix[16];
float model_view_matrix[16];
float project_matrix[16];
float model_view_project_matrix[16];

// FBO common rendering
unsigned int fbo_common_displaying;
unsigned int texture_common_color_buffer;
unsigned int texture_common_depth_buffer;

// viewpoint candidates
const size_t VIEWPOINT_CANDIDATE_NUMBER = 200;
std::vector<point> viewpoint_candidates;

// global large allocation

// RGB -> derived from face's N.O.
__forceinline void get_color(int i, unsigned char v[3])
{
    if (i > 0x00FFFFFF)
    {
        // since RGB can only store 0x00FFFFFF colors
        std::cout<<"Color conversion out of range "<<i<<std::endl;
    }
    i += 1;
    v[0] = (i >> 16) & 0x000000FF;
    v[1] = (i >>  8) & 0x000000FF;
    v[2] = (i >>  0) & 0x000000FF;
}

struct file_size // functor for get file size
{
    unsigned int operator()(std::string fn)
    {
        FILE *pf = NULL;
        fopen_s(&pf, fn.c_str(), "rb");
        fseek(pf, 0, SEEK_END);
        unsigned int size = ftell(pf);
        fclose(pf);
        return size;
    }
};

template<class T>
void save_result(T idx)
{
	// Find first non-used filename
	FILE *f;
    std::string fn = std::string("buffer-") + boost::lexical_cast<string>(idx) + std::string(".ppm");
    fopen_s(&f, fn.c_str(), "wb");
	// Read pixels
	int width = window_width, height = window_height; 
	char *buf = new char[width*height*3];
    memset(buf,196,width*height*3);
	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, buf);

	// Flip top-to-bottom
	for(int i = 0; i < height/2; i++) {
		char *row1 = buf + 3 * width * i;
		char *row2 = buf + 3 * width *(height - 1 - i);
		for(int j = 0; j < 3 * width; j++)
			swap(row1[j], row2[j]);
	}

	// Write out file
	fprintf(f, "P6\n%d %d\n255\n", width, height);
	fwrite(buf, width*height*3, 1, f);
	fclose(f);
	delete [] buf;
}

void init_framebuffer_object()
{
    std::cout<<"Initiate FrameBuffer and Check Status...\n";
    /************************************************************************/
    /*       FBO AND BUFFER FOR COMMON RENDERING                            */
    /************************************************************************/

    glGenFramebuffersEXT(1, &fbo_common_displaying);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo_common_displaying);

    // attach color-buffer for displaying
    glGenTextures(1,&texture_common_color_buffer);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB,texture_common_color_buffer);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(
        GL_TEXTURE_RECTANGLE_ARB,
        0,
        GL_RGBA,
        window_width,
        window_height,
        0,
        GL_RGBA,
        GL_UNSIGNED_BYTE,
        NULL);
    glFramebufferTexture2DEXT(
        GL_FRAMEBUFFER_EXT,
        GL_COLOR_ATTACHMENT0_EXT,
        GL_TEXTURE_RECTANGLE_ARB,
        texture_common_color_buffer,
        0);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB,0);

    // attach depth buffer for displaying
    glGenTextures(1,&texture_common_depth_buffer);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB,texture_common_depth_buffer);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(
        GL_TEXTURE_RECTANGLE_ARB,
        0,
        GL_DEPTH_COMPONENT,
        window_width,
        window_height,
        0,
        GL_DEPTH_COMPONENT,
        GL_FLOAT,
        NULL);
    glFramebufferTexture2DEXT(
        GL_FRAMEBUFFER_EXT,
        GL_DEPTH_ATTACHMENT_EXT,
        GL_TEXTURE_RECTANGLE_ARB,
        texture_common_depth_buffer,
        0);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB,0);

    check_frame_buffer_object_status("Common rendering");
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
}

//-------------------------------------------------------------------------
//  Set OpenGL program initial state.
//-------------------------------------------------------------------------
void init()
{	
    glewInit();

	//  Set the frame buffer clear color to black. 
    glClearColor(1,1,1,1);
    glClearDepth(1.0);

    // depth and cull
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    // speedup
    glEnable(GL_DITHER);
    glShadeModel(GL_SMOOTH);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_FASTEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_FASTEST);

    // light
    static const GLfloat light0_color[4] = { 0.8f, 0.8f, 0.8f, 1.0f };
    static const GLfloat light1_color[4] = { 0.4f, 0.4f, 0.8f, 1.0f };
    static const GLfloat light0_pos[4]   = { 100.0f, 100.0f, 100.0f, 0.0f };
    static const GLfloat light1_pos[4]   = { -100.0f, 100.0f, 100.0f, 0.0f };
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  light0_color);
    glLightfv(GL_LIGHT1, GL_DIFFUSE,  light1_color);
    glLightfv(GL_LIGHT0, GL_POSITION, light0_pos);
    glLightfv(GL_LIGHT1, GL_POSITION, light1_pos);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);

    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    // init texture type
    glEnable(GL_TEXTURE_RECTANGLE_ARB);
    glEnable(GL_TEXTURE_2D);
}

//-------------------------------------------------------------------------
//  This function is passed to glutDisplayFunc in order to display 
//	OpenGL contents on the window.
//-------------------------------------------------------------------------
void display(void)
{
    // get result list
    const unsigned int TM = 45; //number of tested models
    std::vector<float> avg_viewpoint_diff;
    std::vector<float> avg_upright_diff;

    std::vector<float> stddev_viewpoint; // standard deviation
    std::vector<float> stddev_upright;

    for (unsigned int i=0; i<TM; i++)
    {
        ifstream rlist("rlist");
        // model i;
        // get trackball quats from testing results 
        std::vector<vec> viewpoints;
        std::vector<vec> uprights;
        while(!rlist.eof()) // result of each tester for model i;
        {
            std::string fn;
            std::getline(rlist,fn);
            if(fn == string("")) break;
            ifstream dat(fn.c_str());
            std::string dump;
            for (unsigned int j=0; j<i; j++)
            {
                // skip i lines
                std::getline(dat,dump);
            }
            // read quat data of i-th model
            float vquat[4];
            dat>>vquat[0]>>vquat[1]>>vquat[2]>>vquat[3];
            dat.close();

            // compute viewpoints and uprights direction
            float vrot[16];
            sgi_trackball_space::build_row_major_rotmatrix(vrot,vquat);
            float vrot_inverse[16];
            inverseMatrix_pm(vrot,vrot_inverse);

            float vlookat[3];
            float vlookat_orig[3] = {0,0,0};
            vlookat_orig[2] = 3*mesh->bsphere.r;
            multMatirxVecotr3(vlookat,vrot_inverse,vlookat_orig);
            float vupright[3];
            float vupright_orig[3] = {0,1,0};
            multMatirxVecotr3(vupright,vrot_inverse,vupright_orig);

            // normalize for later computing
            viewpoints.push_back(vec(vlookat[0],vlookat[1],vlookat[2]));
            uprights.push_back(vec(vupright[0],vupright[1],vupright[2]));

            normalize(viewpoints.back());
            normalize(uprights.back());
        }
        rlist.close();

        size_t ntester = viewpoints.size();
        std::vector<float> angles_viewpoint; 
        std::vector<float> angles_up;
        for (unsigned int j=0; j<ntester; j++)
        {
            for (unsigned int k=0; k<ntester; k++)
            {
                float vcos = viewpoints[j] DOT viewpoints[k];
                vcos = vcos>1?1:vcos; //clamp to [-1.0,1.0], maybe a little overflow to 1.0000001 or something.
                vcos = vcos<-1?-1:vcos;
                angles_viewpoint.push_back(acos(vcos));

                float ucos = uprights[j] DOT uprights[k];
                ucos = ucos>1?1:ucos; //clamp to [-1.0,1.0], maybe a little overflow to 1.0000001 or something.
                ucos = ucos<-1?-1:ucos;
                angles_up.push_back(acos(ucos));
            }
        }

        //get the most stable one;
        const float threshold_viewpoint = acos(-1.0f)/4.0f;
        const float threshold_up = acos(-1.0f)/4.0f;
        std::vector<int> similar_viewpoint;
        std::vector<int> similar_up;
        for (unsigned int j=0; j<ntester; j++)
        {
            int cnt_viewpoint =0;
            int cnt_up =0;
            for (unsigned int k=0; k<ntester; k++)
            {
                if (angles_viewpoint[j*ntester+k] < threshold_viewpoint)
                {
                    cnt_viewpoint++;
                }
                if (angles_up[j*ntester+k] < threshold_up)
                {
                    cnt_up++;
                }
            }
            similar_viewpoint.push_back(cnt_viewpoint);
            similar_up.push_back(cnt_up);
        }
        int ms_viewpoint = 0;
        int ms_no_viewpoint = 0;
        int ms_up = 0;
        int ms_no_up = 0;
        for (unsigned int j=0; j<ntester; j++)
        {
            if (similar_viewpoint[j]>ms_viewpoint)
            {
                ms_viewpoint = similar_viewpoint[j];
                ms_no_viewpoint = j;
            }
            if (similar_up[j]>ms_up)
            {
                ms_up = similar_up[j];
                ms_no_up = j;
            }
        }

        // compute the average difference to the most stable one
        float avg_viewpoint = 0;
        float avg_up = 0;
        for (unsigned int j=0; j<ntester; j++)
        {
            avg_viewpoint += angles_viewpoint[ms_no_viewpoint*ntester + j];
            avg_up += angles_up[ms_no_viewpoint*ntester + j];
        }
        avg_viewpoint /= ntester;
        avg_up /= ntester;
        avg_viewpoint_diff.push_back(avg_viewpoint);
        avg_upright_diff.push_back(avg_up);

        // compute standard deviation of difference
        float sd_viewpoint = 0;
        float sd_up = 0;
        for (unsigned int j=0; j<ntester; j++)
        {
            sd_viewpoint += sqr(angles_viewpoint[ms_no_viewpoint*ntester + j] - avg_viewpoint);
            sd_up += sqr(angles_up[ms_no_viewpoint*ntester + j] - avg_up);
        }
        sd_viewpoint /= ntester -1;
        sd_up /= ntester -1;
        sd_viewpoint = sqrt(sd_viewpoint);
        sd_up = sqrt(sd_up);
        stddev_viewpoint.push_back(sd_viewpoint);
        stddev_upright.push_back(sd_up);

        // save results to files
        // for model i;
        std::string sfn = std::string("rmodel-") + boost::lexical_cast<std::string>(i);
        std::ofstream osf(sfn.c_str());
        osf<<"For model:"<<i<<std::endl;
        osf<<"Viewpoint angle threshold:"<<std::endl<<threshold_viewpoint<<std::endl
            <<"Upright direction threshold:"<<std::endl<<threshold_up<<std::endl;
        osf<<"Viewpoint similar:"<<std::endl;
        for (unsigned int j=0; j<ntester; j++)
        {
            osf<<similar_viewpoint[j]<<' ';
        }
        osf<<std::endl;
        osf<<"The best one:"<<std::endl;
        osf<<ms_no_viewpoint<<std::endl;
        osf<<"Parameters:"<<std::endl;
        osf<<viewpoints[ms_no_viewpoint]<<std::endl;
        osf<<uprights[ms_no_viewpoint]<<std::endl;
        osf<<"Upright similar:"<<std::endl;
        for (unsigned int j=0; j<ntester; j++)
        {
            osf<<similar_up[j]<<' ';
        }
        osf<<std::endl;
        osf<<"The best one:"<<std::endl;
        osf<<ms_no_up<<std::endl;
        osf.close();

        // load model i;
        if(mesh!=0) delete mesh;
        std::wifstream mlist("mlist");
        std::wstring mname;
        for (unsigned int j=0; j<i+1; j++)
        {
            getline(mlist,mname);
        }
        mesh = Mesh::read(mname);
        mesh->compute_bsphere();
        mesh->compute_normals();
        // set eye position to bounding sphere center, it's relative to(0,0,0) coords.
        eye_pos_to_bsphere_center = point(0,0,3*mesh->bsphere.r);
        sgi_trackball_space::trackball(trackball_quat,0,0,0,0); // and also trackball quat

        // render and save the image
        // draw to frame buffer, and redraw to system framebuffer use texture;
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,fbo_common_displaying);
        glViewport(0,0,window_width,window_height);
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

        float trans1[16],rot1[16];  // translate and rotate matrix
        buildTranslateMatrix(
            -mesh->bsphere.center[0],
            -mesh->bsphere.center[1],
            -mesh->bsphere.center[2],
            trans1);
        sgi_trackball_space::build_row_major_rotmatrix(rot1,trackball_quat);
        //set viewpoint according to auto sampling parameter
        buildLookAtMatrix(
            //// eye position
            //eye_pos_to_bsphere_center[0],
            //eye_pos_to_bsphere_center[1],
            //eye_pos_to_bsphere_center[2],
            viewpoints[ms_no_viewpoint][0]*mesh->bsphere.r*3,
            viewpoints[ms_no_viewpoint][1]*mesh->bsphere.r*3,
            viewpoints[ms_no_viewpoint][2]*mesh->bsphere.r*3,
            // look at center
            0,
            0,
            0,
            // upright direction
            //0,1,0,
            uprights[ms_no_up][0],
            uprights[ms_no_up][1],
            uprights[ms_no_up][2],
            view_matrix);
        buildOrthoMatrix(
            -mesh->bsphere.r, mesh->bsphere.r,
            -mesh->bsphere.r, mesh->bsphere.r,
            2*mesh->bsphere.r, 4*mesh->bsphere.r,
            project_matrix);
        //buildPerspectiveMatrix(
        //    45.0, 1.0,
        //    0.1*mesh->bsphere.r, 999*mesh->bsphere.r,
        //    project_matrix );
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

        // draw according to color statistic
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        glEnable(GL_LIGHT1);
        glBegin(GL_TRIANGLES);
        glColor3f(0.3f,0.6f,0.15f);
        for (unsigned int j=0; j<mesh->faces.size(); j++)
        {
            glNormal3f(
                mesh->normals[mesh->faces[j][0]][0],
                mesh->normals[mesh->faces[j][0]][1],
                mesh->normals[mesh->faces[j][0]][2]);
            glVertex3f(
                mesh->vertices[mesh->faces[j][0]][0],
                mesh->vertices[mesh->faces[j][0]][1],
                mesh->vertices[mesh->faces[j][0]][2]);
            glNormal3f(
                mesh->normals[mesh->faces[j][1]][0],
                mesh->normals[mesh->faces[j][1]][1],
                mesh->normals[mesh->faces[j][1]][2]);
            glVertex3f(
                mesh->vertices[mesh->faces[j][1]][0],
                mesh->vertices[mesh->faces[j][1]][1],
                mesh->vertices[mesh->faces[j][1]][2]);
            glNormal3f(
                mesh->normals[mesh->faces[j][2]][0],
                mesh->normals[mesh->faces[j][2]][1],
                mesh->normals[mesh->faces[j][2]][2]);
            glVertex3f(
                mesh->vertices[mesh->faces[j][2]][0],
                mesh->vertices[mesh->faces[j][2]][1],
                mesh->vertices[mesh->faces[j][2]][2]);
        }
        glEnd();
        glFinish();

        // read back buffer
        boost::scoped_array<unsigned char> cb(new unsigned char[window_width*window_height*3]);
        glReadPixels(0, 0, window_width, window_height, GL_RGB, GL_UNSIGNED_BYTE, cb.get());
        //CImgGL<> res;
       // res.assign(window_width,window_height,cb.get());
        std::string rfn = std::string("rmodel-") + boost::lexical_cast<std::string>(i) + std::string(".ppm");
        //res.save(rfn.c_str());
    }

    std::ofstream avgs("avgs");
    for (unsigned int i=0; i<TM; i++)
    {
        avgs<<avg_viewpoint_diff[i]<<' ';
        avgs<<avg_upright_diff[i]<<' ';
        avgs<<stddev_viewpoint[i]<<' ';
        avgs<<stddev_upright[i]<<std::endl;
    }
    avgs.close();

    // system buffer
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
    glViewport(0,0,window_width,window_height);
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glutSwapBuffers();
}

//-------------------------------------------------------------------------
//  This function is passed to the glutReshapeFunc and is called 
//  whenever the window is resized.
//-------------------------------------------------------------------------
void reshape(int w, int h)
{
    //  Stay updated with the window width and height
    window_width = w;
    window_height = h;

    //  Reset viewport
    glViewport(0, 0, window_width, window_height);

    //  Print current width and height on the screen
    if(MOUSE_KEYBOARD_VERBOSE)
        std::cout
        <<"Window Width: "<<window_width<<std::endl
        <<"Window Height: "<<window_height<<std::endl;

    //glutHideWindow();
}

//-------------------------------------------------------------------------
//  This function is passed to the glutMouseFunc and is called 
//  whenever the mouse is clicked.
//-------------------------------------------------------------------------
void mouse(int button, int state, int x, int y)
{
    switch(button)
    {
        //  Left Button Clicked
    case GLUT_LEFT_BUTTON:

        switch(state)
        {
            //  Pressed 
        case GLUT_DOWN:
            if(MOUSE_KEYBOARD_VERBOSE)
                std::cout<<"Mouse Left Button Pressed(Down)..."<<std::endl;
            break;
            //  Released
        case GLUT_UP:
            if(MOUSE_KEYBOARD_VERBOSE)
                std::cout<<"Mouse Left Button Released(Up)..."<<std::endl;
            break;
        }

        break;

        //  Middle Button clicked
    case GLUT_MIDDLE_BUTTON:

        switch(state)
        {
            //  Pressed
        case GLUT_DOWN:
            if(MOUSE_KEYBOARD_VERBOSE)
                std::cout<<"Mouse Middle Button Pressed(Down)..."<<std::endl;
            break;
            //  Released
        case GLUT_UP:
            if(MOUSE_KEYBOARD_VERBOSE)
                std::cout<<"Mouse Middle Button Released(Up)..."<<std::endl;
            break;
        }

        break;

        //  Right Button Clicked
    case GLUT_RIGHT_BUTTON:

        switch(state)
        {
            //  Pressed
        case GLUT_DOWN:
            if(MOUSE_KEYBOARD_VERBOSE)
                std::cout<<"Mouse Right Button Pressed(Down)..."<<std::endl;
            break;
            //  Released
        case GLUT_UP:
            if(MOUSE_KEYBOARD_VERBOSE)
                std::cout<<"Mouse Right Button Released(Up)..."<<std::endl;
            break;
        }

        break;

    }
}

//-------------------------------------------------------------------------
//  This function is passed to the glutMotionFunc and is called 
//  whenever the mouse is dragged.
//-------------------------------------------------------------------------
void motion(int x, int y)
{
    //  Print the mouse drag position
    //if(MOUSE_KEYBOARD_VERBOSE)
    //    std::cout<<"Mouse Drag Position: "<<x<<' '<<y<<std::endl;
}

//-------------------------------------------------------------------------
//  This function is passed to the glutPassiveMotionFunc and is called 
//  whenever the mouse is moved.
//-------------------------------------------------------------------------
void pmotion(int x, int y)
{
    //  Print mouse move positopn
    //if(MOUSE_KEYBOARD_VERBOSE)
    //    std::cout<<"Mouse Move Position: "<<x<<' '<<y<<std::endl;
}

//-------------------------------------------------------------------------
//  This function is passed to the glutKeyboardFunc and is called 
//  whenever the user hits a key.
//-------------------------------------------------------------------------
void keyboard(unsigned char key, int x, int y)
{
    //  Print what key the user is hitting
    if(MOUSE_KEYBOARD_VERBOSE)
        std::cout<<"User is hitting the "<<key<<" key."<<std::endl;
    if(MOUSE_KEYBOARD_VERBOSE)
        std::cout<<"ASCII code is "<<unsigned int(key)<<std::endl;

    switch(key)
    {
        //  User hits A key
    case 'a':

        break;

        //  User hits Shift + A key
    case 'A':

        break;

        //  User hits Enter
    case '\r':
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"User is hitting the Return key."<<std::endl; 
        break;

        //  User hits Space
    case ' ':
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"User is hitting the Space key."<<std::endl; 
        break;

        //  User hits back space
    case 8:
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"User is hitting the Back Space key."<<std::endl; 
        break;

        //  User hits ESC key
    case 27:
        exit(1);
        break;
    }

    glutPostRedisplay();
}

//-------------------------------------------------------------------------
//  This function is passed to the glutSpecialFunc and is called 
//  whenever the user hits a special key.
//-------------------------------------------------------------------------
void special(int key, int x, int y)
{
    switch(key)
    {
    case GLUT_KEY_F1 :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"F1 function key."<<std::endl; 
        break;
    case GLUT_KEY_F2 :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"F2 function key."<<std::endl;  
        break;
    case GLUT_KEY_F3 :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"F3 function key."<<std::endl;  
        break;
    case GLUT_KEY_F4 :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"F4 function key."<<std::endl;  
        break;
    case GLUT_KEY_F5 :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"F5 function key."<<std::endl;  
        break;
    case GLUT_KEY_F6 :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"F6 function key."<<std::endl;  
        break;
    case GLUT_KEY_F7 :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"F7 function key."<<std::endl;  
        break;
    case GLUT_KEY_F8 :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"F8 function key."<<std::endl;  
        break;
    case GLUT_KEY_F9 :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"F9 function key."<<std::endl;  
        break;
    case GLUT_KEY_F10 :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"F10 function key."<<std::endl;  
        break;
    case GLUT_KEY_F11 :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"F11 function key."<<std::endl;  
        break;
    case GLUT_KEY_F12 :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"F12 function key."<<std::endl;  
        break;
    case GLUT_KEY_LEFT :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"Left directional key."<<std::endl;  
        break;
    case GLUT_KEY_UP :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"Up directional key."<<std::endl;  
        break;
    case GLUT_KEY_RIGHT :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"Right directional key."<<std::endl;  
        break;
    case GLUT_KEY_DOWN :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"Down directional key."<<std::endl;  
        break;
    case GLUT_KEY_PAGE_UP :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"Page up directional key."<<std::endl;  
        break;
    case GLUT_KEY_PAGE_DOWN :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"Page down directional key."<<std::endl;  
        break;
    case GLUT_KEY_HOME :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"Home directional key."<<std::endl;  
        break;
    case GLUT_KEY_END :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"End directional key."<<std::endl;  
        break;
    case GLUT_KEY_INSERT :
        if(MOUSE_KEYBOARD_VERBOSE)
            std::cout<<"Inset directional key."<<std::endl;  
        break;
    }

    glutPostRedisplay();
}

//-------------------------------------------------------------------------
//  This function sets the window x and y coordinates
//  such that the window becomes centered
//-------------------------------------------------------------------------
void center_window()
{
    window_x = (glutGet (GLUT_SCREEN_WIDTH) - window_width)/2;
    window_y = (glutGet (GLUT_SCREEN_HEIGHT) - window_height)/2;
}

//-------------------------------------------------------------------------
//  Program Main method.
//-------------------------------------------------------------------------
void main(int argc, char **argv)
{
    /************************************************************************/
    /*   Usage: ve.exe width height -cf vpc.dat -f target                   */
    /************************************************************************/
    //if(argc == 1)
    //{
    //    std::cout<<"Usage: ve.exe target_file -w width -h height -vc vcf\n"
    //        <<"target_file: model file\n"
    //        <<"-w width: screen resolutioin width\n"
    //        <<"-w height: screen resolutioin height\n"
    //        <<"-vc vcf: viewpoint candidate data file, when using static viewpoint candidates.\n"
    //        <<std::endl;
    //    return -1;
    //}

    // generate or load viewpoint candidates
    generate_viewpoint_candidates(VIEWPOINT_CANDIDATE_NUMBER, viewpoint_candidates);

    glutInit(&argc, argv);
    center_window();
    glutInitWindowSize(window_width, window_height);
    glutInitWindowPosition(window_x, window_y);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutCreateWindow(window_title);

    //  View in full screen if the full_screen flag is on
    if(full_screen) glutFullScreen ();

    //  Set OpenGL program initial state.
    init();

    // Set the callback functions
    glutDisplayFunc(display);
    glutReshapeFunc (reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutPassiveMotionFunc(pmotion);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special);

    // set run option after close the window
    //glutSetOption( GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION ) ;

    // init framebuffer
    init_framebuffer_object();

    // load model data and set selection conds.
    std::wstring fn;
    char *pfn = *(argv+1);
    while (*pfn != '\0') fn.push_back(*pfn++);
    mesh = Mesh::read(fn);
    mesh->compute_bsphere();
    mesh->compute_normals();
    // set eye position to bounding sphere center, it's relative to(0,0,0) coords.
    eye_pos_to_bsphere_center = point(0,0,3*mesh->bsphere.r);
    sgi_trackball_space::trackball(trackball_quat,0,0,0,0); // and also trackball quat

	//  Start GLUT event processing loop
	//glutMainLoop();
    glutMainLoopEvent();

    std::cout<<"I am back from glut main loop...\n";
}

