#include <gl/glew.h>
#include <gl/freeglut.h>

#include <iostream>
#include <vector>
#include <string>
#include <boost/smart_cast.hpp>
#include <boost/lexical_cast.hpp>

#include <vec.h>
#include "../common/com.h"
#include "../common/mesh.h"
#include "../Common/trackball.h"
#include "../Common/errorhandle.h"

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
char *window_title = "Viewpoint Entropy";

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
const size_t VIEPOINT_CANDIDATE_NUMBER = 3;
std::vector<point> viewpoint_candidates;

void save_result(int idx)
{
	// Find first non-used filename
	FILE *f;
    std::string fn = std::string("buffer") + boost::lexical_cast<string>(idx) + std::string(".ppm");
    fopen_s(&f, fn.c_str(), "wb");
	// Read pixels
	int width = window_width, height = window_height; 
	char *buf = new char[width*height*3];
//	glPixelStorei(GL_PACK_ALIGNMENT, 1);
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
	glClearColor(0.0, 0.0, 0.0, 0.0);
    glClearDepth(1.0);

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
    glEnable(GL_TEXTURE_RECTANGLE_ARB);
    glEnable(GL_TEXTURE_2D);
}

//-------------------------------------------------------------------------
//  This function is passed to glutDisplayFunc in order to display 
//	OpenGL contents on the window.
//-------------------------------------------------------------------------
void display(void)
{
    // draw to frame buffer, and redraw to system framebuffer use texture;
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,fbo_common_displaying);
    glViewport(0,0,window_width,window_height);

    for (unsigned int i=0; i<VIEPOINT_CANDIDATE_NUMBER; i++)
    {
        glClearColor(1,1,1,1);
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

        const float theta = 60.0;
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        gluLookAt(
            viewpoint_candidates[i][0], viewpoint_candidates[i][1], viewpoint_candidates[i][2],
            0,0,0,
            0,1,0);

        glColor3d(0.6,0.4,0.2);
        glEnable(GL_LIGHTING);
        glBegin(GL_TRIANGLES);
        glVertex3f(0,0,0);
        glVertex3f(1,0,0);
        glVertex3f(0,1,0);
        glEnd();
        glDisable(GL_LIGHTING);

        std::cout<<viewpoint_candidates[i][0]<<viewpoint_candidates[i][1]<<viewpoint_candidates[i][2];
        save_result(i);
    }

    // select the best one
    // render to system frame buffer
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
    glViewport(0,0,window_width,window_height);
    glClearColor(1,1,1,1);
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1,1,-1,1,-1,1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glBindTexture(GL_TEXTURE_RECTANGLE_ARB,texture_common_color_buffer);
    glColor3f(1,1,1);   // set color to original color
    glBegin(GL_QUADS);  // draw color buffer 0
    {
        glTexCoord2i(0, 0);                             glVertex3f(-1, -1, 0);
        glTexCoord2i(window_width, 0);                  glVertex3f(1, -1, 0);
        glTexCoord2i(window_width, window_height);      glVertex3f(1, 1, 0);
        glTexCoord2i(0, window_height);                 glVertex3f(-1, 1, 0);

        //glTexCoord2i(0, 0);                             glVertex3f(-1, -1, 0);
        //glTexCoord2i(window_width, 0);                  glVertex3f(0, -1, 0);
        //glTexCoord2i(window_width, window_height);      glVertex3f(0, 0, 0);
        //glTexCoord2i(0, window_height);                 glVertex3f(-1, 0, 0);
    }
    glEnd();
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB,0);

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
    generate_viewpoint_candidates(VIEPOINT_CANDIDATE_NUMBER, viewpoint_candidates);

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
    //mesh = Mesh::read(argv+1);
    //mesh->compute_bsphere();
    //mesh->compute_normals();

    //reset eye position to bounding sphere center, it's relative to(0,0,0) coords.
    //eye_pos_to_bsphere_center = point(0,0,3*mesh->bsphere.r);
    //sgi_trackball_space::trackball(trackball_quat,0,0,0,0); // and also trackball quat

	//  Start GLUT event processing loop
	//glutMainLoop();
    glutMainLoopEvent();

    std::cout<<"I am back from glut main loop...\n";
}

