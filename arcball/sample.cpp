// sample.cpp
//   sample program demonstrating use of the arcball.
//   see arcball.h for an explanation of what the
//   arcball_*** functions do.
//
//   must be compiled with glut
//
//   -- Bradley Smith, 5/05/2008

#include <GL/glut.h>
#include <cstdlib>
#include <cmath>
#include "arcball.h"

// =================
// Drawing the Scene
// =================

static float aspect_ratio = 1.0f;
static int width, height;

// scene parameters
const vec_arcball eye( 0.0f, 0.0f, -20.0f );
const vec_arcball centre( 0.0f, 0.0f, 0.0f );
const vec_arcball up( 0.0f, 1.0f, 0.0f );
const float SPHERE_RADIUS = 5.0f;
const int SPHERE_LAT_SLICES = 12;
const int SPHERE_LONG_SLICES = 24;
const int NUM_STARS = 256;
static vec_arcball star[NUM_STARS];

const float PI = 3.141592654f;

static void reset_view(int w, int h)
{
	width = w;
	height = h;
    aspect_ratio = (float) width / (float) height;

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective( 50.0f, aspect_ratio, 1.0f, 50.0f );
    gluLookAt(
        eye.x, eye.y, eye.z,
        centre.x, centre.y, centre.z,
        up.x, up.y, up.z );
	// set up the arcball using the current projection matrix
	arcball_setzoom( SPHERE_RADIUS, eye, up );

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity() ;
}

inline float randf()
{
	return ((1.0f / 127.f) * (((float)(rand() & 255)) - 128.0f)) ;
}

static void startup_scene()
{
	// generate random starfield
	for ( int i=0; i < NUM_STARS; ++i )
	{
		star[i] = vec_arcball( randf(), randf(), randf() ).unit() * 2.0f;
	}
}

static void shutdown_scene()
{
	// nothing to be done here
}

static void draw_stars()
{
	glBegin( GL_POINTS );
	glColor3f( 1.0f, 1.0f, 1.0f );
	for ( int i=0; i < NUM_STARS; ++i )
	{
		glVertex3f( star[i].x, star[i].y, star[i].z );
	}
	glEnd();
}

inline vec_arcball rotate_x( vec_arcball v, float sin_ang, float cos_ang )
{
	return vec_arcball(
	    v.x,
	    (v.y * cos_ang) + (v.z * sin_ang),
	    (v.z * cos_ang) - (v.y * sin_ang)
	    );
}

inline vec_arcball rotate_y( vec_arcball v, float sin_ang, float cos_ang )
{
	return vec_arcball(
	    (v.x * cos_ang) + (v.z * sin_ang),
	    v.y,
	    (v.z * cos_ang) - (v.x * sin_ang)
	    );
}

static void draw_sphere()
{
	const int NUM_FACE_COLOURS = 4;
	vec_arcball FACE_COLOURS[NUM_FACE_COLOURS] = {
		vec_arcball( 1.0f, 0.0f, 0.0f ),
		vec_arcball( 0.0f, 1.0f, 0.0f ),
		vec_arcball( 0.0f, 0.0f, 1.0f ),
		vec_arcball( 1.0f, 1.0f, 0.0f )
	};

	const float lat_angle = PI / (float)SPHERE_LAT_SLICES;
	const float long_angle = 2.0f * PI / (float)SPHERE_LONG_SLICES;

	const float sin_lat = sin( lat_angle );
	const float cos_lat = cos( lat_angle );
	const float sin_long = sin( long_angle );
	const float cos_long = cos( long_angle );
	
	glBegin( GL_QUADS );

	vec_arcball lat_0( 0.0f, SPHERE_RADIUS, 0.0 );
	
	for ( int y = 0; y < SPHERE_LAT_SLICES; ++y )
	{
		vec_arcball lat_1 = rotate_x( lat_0, sin_lat, cos_lat );

		vec_arcball long_0_0 = lat_0;
		vec_arcball long_1_0 = lat_1;
		
		for ( int x = 0; x < SPHERE_LONG_SLICES; ++x )
		{
			vec_arcball long_0_1 = rotate_y( long_0_0, sin_long, cos_long );
			vec_arcball long_1_1 = rotate_y( long_1_0, sin_long, cos_long );
			
			vec_arcball colour = FACE_COLOURS[(x+y)%NUM_FACE_COLOURS];
			glColor4f( colour.x, colour.y, colour.z, 1.0f );
				
			glVertex3f( long_0_0.x, long_0_0.y, long_0_0.z );
			glVertex3f( long_1_0.x, long_1_0.y, long_1_0.z );
			glVertex3f( long_1_1.x, long_1_1.y, long_1_1.z );
			glVertex3f( long_0_1.x, long_0_1.y, long_0_1.z );

			long_0_0 = long_0_1;
			long_1_0 = long_1_1;
		}
		
		lat_0 = lat_1;
	}
	
	glEnd();

	// Something simpler:
    //glColor3f(1.0f,0.0f,0.0f);
	//glutWireSphere(SPHERE_RADIUS,SPHERE_LONG_SLICES,SPHERE_LAT_SLICES);
}

static void click_scene(int x, int y)
{
	int invert_y = (height - y) - 1; // OpenGL viewport coordinates are Cartesian
	arcball_start(x,invert_y);
}

static void drag_scene(int x, int y)
{
	int invert_y = (height - y) - 1;
	arcball_move(x,invert_y);
}

static void draw_scene()
{
	// stars: to simulate infinite distance,
	//        translate the sphere of stars to the eye
	//        and perform the arcball rotation around the eye
	//        (also disable depth test so they become the background)
	glPushMatrix();
	glDisable( GL_DEPTH_TEST );
	glTranslatef( eye.x, eye.y, eye.z );
    arcball_rotate();
	draw_stars();
	glEnable( GL_DEPTH_TEST );
	glPopMatrix();

	// now render the regular scene under the arcball rotation about 0,0,0
	// (generally you would want to render everything here)
    arcball_rotate();
	draw_sphere();
}

// ==============
// GLUT Callbacks
// ==============

static void resize(int w, int h)
{
	reset_view(w,h);
}

static void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    draw_scene();

    glutSwapBuffers();
}

static void key(unsigned char key, int x, int y)
{
    switch (key) 
    {
        case 27 : 
        case 'q':
			shutdown_scene();
            exit(0);
            break;
		default :
			break;
    }
    glutPostRedisplay();
}

static void idle()
{
    glutPostRedisplay();
}

static void mouse_button(int button, int state, int x, int y)
{
	if ( state == GLUT_DOWN ) click_scene(x,y);
}

static void mouse_motion(int x, int y)
{
	// glutMotionFunc only called when a mouse button is held down
	drag_scene(x,y);
}

// ====
// Main
// ====

int main(int argc, char ** argv)
{
	arcball_reset();
	startup_scene();

    glutInit(&argc, argv);
    glutInitWindowSize(400,400);
    glutInitWindowPosition(10,10);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("Arcball Sample");

    glutReshapeFunc(resize);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutIdleFunc(idle);
    glutMouseFunc(mouse_button);
    glutMotionFunc(mouse_motion);

    glClearColor(0.0f,0.0f,0.0f,1.0f);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    glutMainLoop();

    return 0; // never reached
}

// end of file
