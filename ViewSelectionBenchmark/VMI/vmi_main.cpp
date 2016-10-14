#include <gl/glew.h>
#include <gl/freeglut.h>

#include <iostream>
#include <vector>
#include <string>
#include <boost/serialization/smart_cast.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/smart_ptr.hpp>

#include <vec.h>
#include "../common/com.h"
#include "../common/mesh.h"
#include "../Common/trackball.h"
#include "../Common/errorhandle.h"
#include "../Common/glmatrix.h"

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
char *window_title = "Viewpoint Mutual Information";

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
char * model_no = 0;

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

    // init texture type
    glEnable(GL_TEXTURE_RECTANGLE_ARB);
    glEnable(GL_TEXTURE_2D);

    // disable anti-aliasing for accurate color projection
    glDisable(GL_LINE_SMOOTH);
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

    // buffer for save viewpoint entropy 
    float viewpoint_mutual_information[VIEWPOINT_CANDIDATE_NUMBER];
    memset(viewpoint_mutual_information,0,VIEWPOINT_CANDIDATE_NUMBER*sizeof(float));

    // saving p(o|v) buffer
    std::vector<std::vector<float> > povs;
    povs.resize(VIEWPOINT_CANDIDATE_NUMBER);
    for (unsigned int i=0; i<VIEWPOINT_CANDIDATE_NUMBER; i++)
    {
        povs[i].resize(mesh->faces.size());
    }

    clock_t start_t = clock();
    for (unsigned int iv=0; iv<VIEWPOINT_CANDIDATE_NUMBER; iv++)
    {
        // color statistic buffer, and color buffer read back.
        std::vector<unsigned int> csb;
        csb.resize(mesh->faces.size());

        // normalize the color statistic, for computing the probability
        std::vector<float> ncsb;
        ncsb.resize(mesh->faces.size());

        // color buffer
        boost::scoped_array<unsigned char> cb(new unsigned char[window_width*window_height*3]);

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
            // eye position
            //eye_pos_to_bsphere_center[0],
            //eye_pos_to_bsphere_center[1],
            //eye_pos_to_bsphere_center[2],
            viewpoint_candidates[iv][0]*3*mesh->bsphere.r,
            viewpoint_candidates[iv][1]*3*mesh->bsphere.r,
            viewpoint_candidates[iv][2]*3*mesh->bsphere.r,
            // look at center
            0,
            0,
            0,
            // upright direction
            0,1,0,
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
        unsigned char vc[3];
        glBegin(GL_TRIANGLES);
        for (unsigned int i=0; i<mesh->faces.size(); i++)
        {
            get_color(i,vc);
            glColor3ub(vc[0],vc[1],vc[2]);
            glVertex3f(
                mesh->vertices[mesh->faces[i][0]][0],
                mesh->vertices[mesh->faces[i][0]][1],
                mesh->vertices[mesh->faces[i][0]][2]);
            glVertex3f(
                mesh->vertices[mesh->faces[i][1]][0],
                mesh->vertices[mesh->faces[i][1]][1],
                mesh->vertices[mesh->faces[i][1]][2]);
            glVertex3f(
                mesh->vertices[mesh->faces[i][2]][0],
                mesh->vertices[mesh->faces[i][2]][1],
                mesh->vertices[mesh->faces[i][2]][2]);
        }
        glEnd();
        glFinish();

        //save_result(iv);

        // compute viewpoint entropy 
        // read back buffer
        glReadPixels(0, 0, window_width, window_height, GL_RGB, GL_UNSIGNED_BYTE, cb.get());
        for (int i=0; i<window_width*window_height; i++)
        {
            if (cb[i*3 + 0] == 255 &&
                cb[i*3 + 1] == 255 &&
                cb[i*3 + 2] == 255)
                continue;

            unsigned int idx =
                (cb[i*3 + 0] << 16) + 
                (cb[i*3 + 1] <<  8) + 
                (cb[i*3 + 2] <<  0) - 1;
            csb[idx] += 1;
        }

        // finally, compute the weight according to the view direction
        // and compute the probability
        for (unsigned int i=0; i<csb.size(); i++)
        {
            // skip backward faces
            if (csb[i] == 0) continue;

            // get i-th face center, relative to the transformed model center
            point  fc;
            fc[0] = (mesh->vertices[mesh->faces[i][0]][0] + 
                mesh->vertices[mesh->faces[i][1]][0] + 
                mesh->vertices[mesh->faces[i][2]][0])/3.0f - mesh->bsphere.center[0];
            fc[1] = (mesh->vertices[mesh->faces[i][0]][1] + 
                mesh->vertices[mesh->faces[i][1]][1] + 
                mesh->vertices[mesh->faces[i][2]][1])/3.0f - mesh->bsphere.center[1];
            fc[2] = (mesh->vertices[mesh->faces[i][0]][2] + 
                mesh->vertices[mesh->faces[i][1]][2] + 
                mesh->vertices[mesh->faces[i][2]][2])/3.0f - mesh->bsphere.center[2];
            vec view_direction = point(
                fc[0] - viewpoint_candidates[iv][0],
                fc[1] - viewpoint_candidates[iv][1],
                fc[2] - viewpoint_candidates[iv][2]);
            vec origin_direction = point(
                0.0f - viewpoint_candidates[iv][0],
                0.0f - viewpoint_candidates[iv][1],
                0.0f - viewpoint_candidates[iv][2]);

            // angle between the two directions
            normalize(view_direction); normalize(origin_direction);
            float cos_alpha = view_direction DOT origin_direction;

            // get the weight for transform area on viewing plane to viewing direction sphere
            // since each facet of the model is very little, we can use the approximation
            // it's to test cos_alpha^2 or cos_alpha^4 which is better.
            float w = cos_alpha * cos_alpha; // to direction perpendicular to view direction
            w *= cos_alpha * cos_alpha; // to view direction sphere

            // last, get the weighted area.
            ncsb[i] = w*csb[i];
        }

        // compute probability
        float sum_ncsb = 0;
        for (unsigned int i=0; i<ncsb.size(); i++)
        {
            sum_ncsb += ncsb[i];
        }

        // compute probability, and save to p(o|v) buffer
        for (unsigned int i=0; i<ncsb.size(); i++)
        {
            povs[iv][i] = ncsb[i]/sum_ncsb;
        }
    }

    // compute pos
    std::vector<float> pos;
    pos.resize(mesh->faces.size());
    for (unsigned int i=0; i<pos.size(); i++)
    {
        for (unsigned int j=0; j<VIEWPOINT_CANDIDATE_NUMBER; j++)
        {
            pos[i] += povs[j][i];
        }
        pos[i] /= VIEWPOINT_CANDIDATE_NUMBER;
    }

    // compute the viewpoint mutual information
    for (unsigned int i=0; i<VIEWPOINT_CANDIDATE_NUMBER; i++)
    {
        for (unsigned int j=0; j<pos.size(); j++)
        {
            if(povs[i][j] > 0 && pos[j] >0)
                viewpoint_mutual_information[i] += povs[i][j]*log(povs[i][j]/pos[j])/log(2.0f);
        }
        std::cout<<viewpoint_mutual_information[i]<<std::endl;
    }
    
    // select the best one
    float mve = viewpoint_mutual_information[0]; // minimum viewpoint entropy
    unsigned int bve_no = 0; //best viewpoint no.
    for (unsigned int i=0; i<VIEWPOINT_CANDIDATE_NUMBER; i++)
    {
        if (viewpoint_mutual_information[i] < mve)
        {
            mve = viewpoint_mutual_information[i];
            bve_no = i;
        }
    }
    std::cout<<"The best one of "<<model_no<<" is: "<<bve_no<<std::endl;

    clock_t end_t = clock();
    float tcost = (end_t - start_t)/(float)CLOCKS_PER_SEC;
    std::ofstream vmiresult("vmiresult",std::ios::out | std::ios::app);
    vmiresult<<tcost<<viewpoint_candidates[bve_no]<<std::endl;
    vmiresult.close();

    // render to system frame buffer
    //glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
    //glViewport(0,0,window_width,window_height);

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
        // eye position
        //eye_pos_to_bsphere_center[0],
        //eye_pos_to_bsphere_center[1],
        //eye_pos_to_bsphere_center[2],
        viewpoint_candidates[bve_no][0]*3*mesh->bsphere.r,
        viewpoint_candidates[bve_no][1]*3*mesh->bsphere.r,
        viewpoint_candidates[bve_no][2]*3*mesh->bsphere.r,
        // look at center
        0,
        0,
        0,
        // upright direction
        0,1,0,
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
    unsigned char vc[3];
    glBegin(GL_TRIANGLES);
    for (unsigned int i=0; i<mesh->faces.size(); i++)
    {
        get_color(i,vc);
        glColor3ub(vc[0],vc[1],vc[2]);
        glVertex3f(
            mesh->vertices[mesh->faces[i][0]][0],
            mesh->vertices[mesh->faces[i][0]][1],
            mesh->vertices[mesh->faces[i][0]][2]);
        glVertex3f(
            mesh->vertices[mesh->faces[i][1]][0],
            mesh->vertices[mesh->faces[i][1]][1],
            mesh->vertices[mesh->faces[i][1]][2]);
        glVertex3f(
            mesh->vertices[mesh->faces[i][2]][0],
            mesh->vertices[mesh->faces[i][2]][1],
            mesh->vertices[mesh->faces[i][2]][2]);
    }
    glEnd();
    glFinish();
    save_result(model_no);

    // system buffer
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
    glViewport(0,0,window_width,window_height);
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
    if (argc != 3)
    {
        std::cout<<"vmi.exe model-file model-no"<<std::endl;
        return;
    }
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
    model_no = *(argv+2);
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

