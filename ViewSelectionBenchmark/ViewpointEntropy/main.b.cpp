#include <gl/freeglut.h>
#include <cstdio>
#include <vector>
#include <string>
#include <boost/smart_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include "vec.h"

#include "../common/com.h"
#include "../common/mesh.h"
#include "../Common/trackball.h"

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
void centerOnScreen();
void drawObject();

//  define the window position on screen
int window_x;
int window_y;

//  variables representing the window size
int window_width = 512;
int window_height = 512;

//  variable representing the window title
char *window_title = "Viewpoint Entropy";

//  Tells whether to display the window full screen or not
//  Press Alt + Esc to exit a full screen.
int full_screen = 0;

//  Viewing parameters
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

// viewpoint candidates
const size_t VIEPOINT_CANDIDATE_NUMBER = 200;
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
//-------------------------------------------------------------------------
//  Set OpenGL program initial state.
//-------------------------------------------------------------------------
void init()
{	
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
    //glEnable(GL_TEXTURE_RECTANGLE_ARB);
    glEnable(GL_TEXTURE_2D);
}

//-------------------------------------------------------------------------
//  This function is passed to glutDisplayFunc in order to display 
//	OpenGL contents on the window.
//-------------------------------------------------------------------------
void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	//  Draw object
	drawObject();

	//  Swap contents of backward and forward frame buffers
	glutSwapBuffers();
}

//-------------------------------------------------------------------------
//  Draws our object.
//-------------------------------------------------------------------------
void drawObject()
{
    // draw...
    static int i = 0;
    const float theta = 60.0;
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glRotated(i*theta, 0,1,0);
    glColor3d(0.6,0.4,0.2);

    glEnable(GL_LIGHTING);
    glutSolidTeapot(0.5);
    save_result(i);
    printf("%d",i++);
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
	printf("Window Width: %d, Window Height: %d.\n", window_width, window_height);
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
					printf("Mouse Left Button Pressed(Down)...\n");
					break;
				//  Released
				case GLUT_UP:
					printf("Mouse Left Button Released(Up)...\n");
					break;
			}

			break;

		//  Middle Button clicked
		case GLUT_MIDDLE_BUTTON:
			
			switch(state)
			{
				//  Pressed
				case GLUT_DOWN:
					printf("Mouse Middle Button Pressed(Down)...\n");
					break;
				//  Released
				case GLUT_UP:
					printf("Mouse Middle Button Released(Up)...\n");
					break;
			}

			break;

		//  Right Button Clicked
		case GLUT_RIGHT_BUTTON:
			
			switch(state)
			{
				//  Pressed
				case GLUT_DOWN:
					printf("Mouse Right Button Pressed(Down)...\n");
					break;
				//  Released
				case GLUT_UP:
					printf("Mouse Right Button Released(Up)...\n");
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
	printf("Mouse Drag Position: %d, %d.\n", x, y);
}

//-------------------------------------------------------------------------
//  This function is passed to the glutPassiveMotionFunc and is called 
//  whenever the mouse is moved.
//-------------------------------------------------------------------------
void pmotion(int x, int y)
{
	//  Print mouse move positopn
	//printf("Mouse Move Position: %d, %d.\n", x, y);
}

//-------------------------------------------------------------------------
//  This function is passed to the glutKeyboardFunc and is called 
//  whenever the user hits a key.
//-------------------------------------------------------------------------
void keyboard(unsigned char key, int x, int y)
{
	//  Print what key the user is hitting
	printf("User is hitting the '%c' key.\n", key);
	printf("ASCII code is %d.\n", key);
	
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
			printf("User is hitting the Return key.\n"); 
			break;

		//  User hits Space
		case ' ':
			printf("User is hitting the Space key.\n"); 
			break;

		//  User hits back space
		case 8:
			printf("User is hitting the Back Space key.\n"); 
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
			printf("F1 function key.\n"); 
			break;
		case GLUT_KEY_F2 :
			printf("F2 function key. \n");  
			break;
		case GLUT_KEY_F3 :
			printf("F3 function key. \n");  
			break;
		case GLUT_KEY_F4 :
			printf("F4 function key. \n");  
			break;
		case GLUT_KEY_F5 :
			printf("F5 function key. \n");  
			break;
		case GLUT_KEY_F6 :
			printf("F6 function key. \n");  
			break;
		case GLUT_KEY_F7 :
			printf("F7 function key. \n");  
			break;
		case GLUT_KEY_F8 :
			printf("F8 function key. \n");  
			break;
		case GLUT_KEY_F9 :
			printf("F9 function key. \n");  
			break;
		case GLUT_KEY_F10 :
			printf("F10 function key. \n");  
			break;
		case GLUT_KEY_F11 :
			printf("F11 function key. \n");  
			break;
		case GLUT_KEY_F12 :
			printf("F12 function key. \n");  
			break;
		case GLUT_KEY_LEFT :
			printf("Left directional key. \n");  
			break;
		case GLUT_KEY_UP :
			printf("Up directional key. \n");  
			break;
		case GLUT_KEY_RIGHT :
			printf("Right directional key. \n");  
			break;
		case GLUT_KEY_DOWN :
			printf("Down directional key. \n");  
			break;
		case GLUT_KEY_PAGE_UP :
			printf("Page up directional key. \n");  
			break;
		case GLUT_KEY_PAGE_DOWN :
			printf("Page down directional key. \n");  
			break;
		case GLUT_KEY_HOME :
			printf("Home directional key. \n");  
			break;
		case GLUT_KEY_END :
			printf("End directional key. \n");  
			break;
		case GLUT_KEY_INSERT :
			printf("Inset directional key. \n");  
			break;
	}
	
	glutPostRedisplay();
}

//-------------------------------------------------------------------------
//  This function sets the window x and y coordinates
//  such that the window becomes centered
//-------------------------------------------------------------------------
void centerOnScreen()
{
	window_x = (glutGet(GLUT_SCREEN_WIDTH) - window_width)/2;
	window_y = (glutGet(GLUT_SCREEN_HEIGHT) - window_height)/2;
}

//-------------------------------------------------------------------------
//  Program Main method.
//-------------------------------------------------------------------------
int main(int argc, char **argv)
{
    /************************************************************************/
    /*   Usage: ve.exe width height -cf vpc.dat -f target                   */
    /************************************************************************/
    //if (argc == 1)
    //{
    //    std::cout<<"Usage: ve.exe target_file -w width -h height -vc vcf\n"
    //        <<"target_file: model file\n"
    //        <<"-w width: screen resolutioin width\n"
    //        <<"-w height: screen resolutioin height\n"
    //        <<"-vc vcf: viewpoint candidate data file, when using static viewpoint candidates.\n"
    //        <<std::endl;
    //    return -1;
    //}

    // initiate the contex
	glutInit(&argc, argv);

	//  Set the window x and y coordinates such that the 
	//  window becomes centered
	centerOnScreen();
	//  Connect to the windowing system + create a window
	//  with the specified dimensions and position
	//  + set the display mode + specify the window title.
	glutInitWindowSize(window_width, window_height);
	glutInitWindowPosition(window_x, window_y);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutCreateWindow(window_title);

    // set run option after close the window
    //glutSetOption( GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION ) ;

	//  View in full screen if the full_screen flag is on
	if(full_screen)
		glutFullScreen();

	//  Set OpenGL program initial state.
	init();

	// Set the callback functions
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutPassiveMotionFunc(pmotion);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(special);

    // load model data and set selection conds.
    generate_viewpoint_candidates(viewpoint_candidates);
    //mesh = Mesh::read(argv+1);
    //mesh->compute_bsphere();
    //mesh->compute_normals();

    //reset eye position to bounding sphere center, it's relative to(0,0,0) coords.
    //eye_pos_to_bsphere_center = point(0,0,3*mesh->bsphere.r);
    //sgi_trackball_space::trackball(trackball_quat,0,0,0,0); // and also trackball quat

	//  Start GLUT event processing loop
	glutMainLoop();
    //for(int i=0;i<6;i++)
    //{
    //    glutMainLoopEvent();
    //    glutPostRedisplay();
    //}

    std::cout<<"I am back from glut main loop...\n";
}

