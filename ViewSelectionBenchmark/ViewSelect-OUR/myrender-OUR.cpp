#include <GL/glew.h>
#pragma comment(lib,"glew32.lib")

#include <GL/glut.h>
#pragma comment(lib, "glut32.lib")

#include <ctime>

#include "myrender-OUR.h"
#include "trackball.h"
#include "glmatrix.h"
#include <algorithm>
#include <boost/random.hpp>
#include <atlstr.h> 
#include "gfx/gfx.h"
#include <math.h>

#pragma warning(disable:4244) //damn stupid warning

const float PI = 3.141592654f;

struct StdDeviations
{
	float standard_deviations;
	int idx;
};

bool stdDeviationsLess(const StdDeviations &a,const StdDeviations &b)
{
	return a.standard_deviations < b.standard_deviations;
}

StdDeviations stdDeviations[VIEWPOINT_SAMPLES_NUM];


bool operator==(Pos2D const &_a,Pos2D const &_b)
{
	if (_b.x==_a.x&&_b.y==_a.y) return true;
	else return false;
}

bool operator==(Pos3D const &_a,Pos3D const &_b)
{
	if (_a.x==_b.x&&_a.y==_b.y&&_a.z==_b.z) return true;
	else return false;
}

// colors for showing segmentation
unsigned char seg_colors[][3] = {    // 18 kinds of color from mspaint stand colors
    {255, 224,   8},
    {255, 192,  16},
    {255, 128,  32},
    {255, 96 ,  64},
    {255, 64 ,  96},
    {255, 32 , 128},
    {255, 16 , 192},
    {255, 8  , 224},

    {128, 224,   8},
    {128, 192,  16},
    {128, 128,  32},
    {128, 96 ,  64},
    {128, 64 ,  96},
    {128, 32 , 128},
    {128, 16 , 192},
    {128, 8  , 224},

    { 64, 224,   8},
    { 64, 192,  16},
    { 64, 128,  32},
    { 64, 96 ,  64},
    { 64, 64 ,  96},
    { 64, 32 , 128},
    { 64, 16 , 192},
    { 64, 8  , 224}
};

MyRender_OUR::MyRender_OUR()
{
    sgi_trackball_space::trackball(trackball_quat, 0,0,0,0);
    model_loaded = false;
    showing_type = SHOWING_MODEL;

    showing_primary_skeleton_number_2d = 0;
    showing_primary_skeleton_number_3d = 0;
    showing_final_skeleton_number_2d = 0;
    showing_final_skeleton_number_3d = 0;

    WOS_2D = 16;
    WOS_3D = 16;

    mesh_viewplane = 0;
    mesh_viewing_shpere = 0;
    current_standard_deviation = 0;
    adaptive_box_size = 1;
    number_of_histogram_intervals = 512;

    selecting_parts_by_mouse = false;

	// 颜色信息
	int i=0;
	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.515600;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.531300;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.546900;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.562500;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.578100;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.593800;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.609400;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.625000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.640600;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.656300;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.671900;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.687500;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.703100;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.718800;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.734400;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.750000;
	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.765600;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.781300;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.796900;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.812500;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.828100;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.843800;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.859400;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.875000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.890600;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.906300;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.921900;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.937500;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.953100;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.968800;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.984400;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	1.000000;
	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.015600;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.031300;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.046900;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.062500;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.078100;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.093800;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.109400;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.125000;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.140600;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.156300;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.171900;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.187500;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.203100;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.218800;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.234400;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.250000;	LUT_CourbureClust[i++]=	1.000000;
	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.265600;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.281300;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.296900;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.312500;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.328100;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.343800;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.359400;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.375000;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.390600;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.406300;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.421900;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.437500;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.453100;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.468800;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.484400;	LUT_CourbureClust[i++]=	1.000000;
	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.500000;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.515600;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.531300;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.546900;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.562500;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.578100;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.593800;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.609400;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.625000;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.640600;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.656300;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.671900;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.687500;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.703100;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.718800;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.734400;	LUT_CourbureClust[i++]=	1.000000;
	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.750000;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.765600;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.781300;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.796900;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.812500;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.828100;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.843800;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.859400;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.875000;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.890600;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.906300;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.921900;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.937500;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.953100;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.968800;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.984400;	LUT_CourbureClust[i++]=	1.000000;
	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.015600;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	1.000000;		LUT_CourbureClust[i++]=	0.031300;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.984400;		LUT_CourbureClust[i++]=	0.046900;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.968800;		LUT_CourbureClust[i++]=	0.062500;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.953100;		LUT_CourbureClust[i++]=	0.078100;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.937500;		LUT_CourbureClust[i++]=	0.093800;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.921900;		LUT_CourbureClust[i++]=	0.109400;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.906300;		LUT_CourbureClust[i++]=	0.125000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.890600;		LUT_CourbureClust[i++]=	0.140600;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.875000;		LUT_CourbureClust[i++]=	0.156300;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.859400;		LUT_CourbureClust[i++]=	0.171900;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.843800;		LUT_CourbureClust[i++]=	0.187500;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.828100;		LUT_CourbureClust[i++]=	0.203100;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.812500;		LUT_CourbureClust[i++]=	0.218800;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.796900;		LUT_CourbureClust[i++]=	0.234400;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.781300;
	LUT_CourbureClust[i++]=	0.250000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.765600;		LUT_CourbureClust[i++]=	0.265600;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.750000;		LUT_CourbureClust[i++]=	0.281300;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.734400;		LUT_CourbureClust[i++]=	0.296900;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.718800;		LUT_CourbureClust[i++]=	0.312500;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.703100;		LUT_CourbureClust[i++]=	0.328100;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.687500;		LUT_CourbureClust[i++]=	0.343800;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.671900;		LUT_CourbureClust[i++]=	0.359400;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.656300;		LUT_CourbureClust[i++]=	0.375000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.640600;		LUT_CourbureClust[i++]=	0.390600;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.625000;		LUT_CourbureClust[i++]=	0.406300;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.609400;		LUT_CourbureClust[i++]=	0.421900;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.593800;		LUT_CourbureClust[i++]=	0.437500;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.578100;		LUT_CourbureClust[i++]=	0.453100;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.562500;		LUT_CourbureClust[i++]=	0.468800;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.546900;		LUT_CourbureClust[i++]=	0.484400;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.531300;
	LUT_CourbureClust[i++]=	0.500000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.515600;		LUT_CourbureClust[i++]=	0.515600;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.500000;		LUT_CourbureClust[i++]=	0.531300;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.484400;		LUT_CourbureClust[i++]=	0.546900;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.468800;		LUT_CourbureClust[i++]=	0.562500;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.453100;		LUT_CourbureClust[i++]=	0.578100;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.437500;		LUT_CourbureClust[i++]=	0.593800;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.421900;		LUT_CourbureClust[i++]=	0.609400;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.406300;		LUT_CourbureClust[i++]=	0.625000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.390600;		LUT_CourbureClust[i++]=	0.640600;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.375000;		LUT_CourbureClust[i++]=	0.656300;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.359400;		LUT_CourbureClust[i++]=	0.671900;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.343800;		LUT_CourbureClust[i++]=	0.687500;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.328100;		LUT_CourbureClust[i++]=	0.703100;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.312500;		LUT_CourbureClust[i++]=	0.718800;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.296900;		LUT_CourbureClust[i++]=	0.734400;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.281300;
	LUT_CourbureClust[i++]=	0.750000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.265600;		LUT_CourbureClust[i++]=	0.765600;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.250000;		LUT_CourbureClust[i++]=	0.781300;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.234400;		LUT_CourbureClust[i++]=	0.796900;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.218800;		LUT_CourbureClust[i++]=	0.812500;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.203100;		LUT_CourbureClust[i++]=	0.828100;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.187500;		LUT_CourbureClust[i++]=	0.843800;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.171900;		LUT_CourbureClust[i++]=	0.859400;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.156300;		LUT_CourbureClust[i++]=	0.875000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.140600;		LUT_CourbureClust[i++]=	0.890600;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.125000;		LUT_CourbureClust[i++]=	0.906300;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.109400;		LUT_CourbureClust[i++]=	0.921900;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.093800;		LUT_CourbureClust[i++]=	0.937500;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.078100;		LUT_CourbureClust[i++]=	0.953100;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.062500;		LUT_CourbureClust[i++]=	0.968800;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.046900;		LUT_CourbureClust[i++]=	0.984400;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.031300;
	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.015600;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.984400;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.968800;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.953100;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.937500;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.921900;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.906300;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.890600;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.875000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.859400;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.843800;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.828100;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.812500;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.796900;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.781300;	LUT_CourbureClust[i++]=	0.000000;
	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.765600;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.750000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.734400;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.718800;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.703100;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.687500;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.671900;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.656300;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.640600;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.625000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.609400;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.593800;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.578100;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.562500;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.546900;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.531300;	LUT_CourbureClust[i++]=	0.000000;
	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.515600;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.500000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.484400;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.468800;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.453100;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.437500;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.421900;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.406300;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.390600;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.375000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.359400;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.343800;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.328100;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.312500;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.296900;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.281300;	LUT_CourbureClust[i++]=	0.000000;
	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.265600;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.250000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.234400;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.218800;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.203100;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.187500;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.171900;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.156300;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.140600;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.125000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.109400;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.093800;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.078100;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.062500;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.046900;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.031300;	LUT_CourbureClust[i++]=	0.000000;
	LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.015600;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	1.000000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.984400;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.968800;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.953100;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.937500;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.921900;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.906300;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.890600;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.875000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.859400;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.843800;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.828100;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.812500;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.796900;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.781300;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;
	LUT_CourbureClust[i++]=	0.765600;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.750000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.734400;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.718800;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.703100;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.687500;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.671900;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.656300;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.640600;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.625000;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.609400;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.593800;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.578100;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.562500;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.546900;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;		LUT_CourbureClust[i++]=	0.531300;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;
	LUT_CourbureClust[i++]=	0.515600;	LUT_CourbureClust[i++]=	0.000000;	LUT_CourbureClust[i++]=	0.000000;


}

void MyRender_OUR::generate_viewpoint_candidates()
{                                            
    //initiate viewpoint candidate
    meta_viewpoint_candidates.clear();

    time_t tm;
    srand((unsigned)time(&tm));  //set seed;


	/*float r = mesh3d->bsphere.r;   //radius of bounding sphere
	r=1.0;
	// sampling a round 
	for (int i=-90; i<360-90; i++)
	{
		float alpha = i/360.0f*2*PI;        //angle with positive-x-axis
		point p;
		p[0] = r*cos(alpha);
		p[2] = 0;
		p[1] = r*sin(alpha);
		meta_viewpoint_candidates.push_back(p);
	}

	return ;*/

#if 0
    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        float alpha = (rand()%360)/360.0f*2*PI;        //angle with positive-x-axis
        float beta = (rand()%360)/360.0f*2*PI - PI;  //angle with xy-plane
        point p;
        p[0] = cos(beta)*sin(alpha);
        p[1] = cos(beta)*cos(alpha);
        p[2] = sin(beta);
        meta_viewpoint_candidates.push_back(p);
    }
    //for (unsigned int i=0; i<30; i++)
    //{
    //    for (unsigned int j=0; j<60; j++)
    //    {
    //        float alpha = i*6/180.0*PI;
    //        float theta = j*6/180.0*PI;
    //        point p;
    //        p[0] = cos(alpha)*sin(theta);
    //        p[1] = cos(alpha)*cos(theta);
    //        p[2] = sin(alpha);
    //        meta_viewpoint_candidates.push_back(p);
    //    }
    //}
#else
    // generation of two series of random number, which are independent
    // when sampling ,use it to multiple bounding sphere radius

    boost::mt19937 rng_a;                 // produces randomness out of thin air
    boost::uniform_01<> dist_a;
    boost::variate_generator<boost::mt19937&, boost::uniform_01<> >
        rand_a(rng_a, dist_a);             // glues randomness with mapping
    rng_a.seed(rand());
    //rng_a.seed(1);
    std::vector<float> ra;
    for (unsigned int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        ra.push_back(rand_a());
    }

    boost::mt11213b rng_b;                 // produces randomness out of thin air
    boost::uniform_01<> dist_b;
    boost::variate_generator<boost::mt11213b&, boost::uniform_01<> >
        rand_b(rng_b, dist_b);             // glues randomness with mapping
    rng_b.seed(rand());
    //rng_b.seed(8);
    std::vector<float> rb;
    for (unsigned int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        rb.push_back(rand_b());
    }
    float theta = 0;
    float phi = 0;
    point p;
    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        theta = 2*acos(sqrt(1-ra[i]));
        phi = 2*PI*rb[i];

        p[0] = sin(theta)*cos(phi);
        p[1] = sin(theta)*sin(phi);
        p[2] = cos(theta);
        meta_viewpoint_candidates.push_back(p);
    }
#endif
}

MyRender_OUR::~MyRender_OUR()
{
	if (mesh_viewing_shpere!=NULL)
		delete mesh_viewing_shpere;
	if (mesh_viewplane!=NULL)
		delete mesh_viewplane;
	vector<ViewPoint*>::iterator iter;
	for (iter = viewpoint_candidates.begin(); iter != viewpoint_candidates.end(); ++iter)
		delete *iter;
	viewpoint_candidates.clear();

	//if (f != NULL)
	//	free(f);

	// common displaying
	glDeleteTextures(1,&texture_common_depth_buffer);
	glDeleteTextures(1,&texture_common_color_buffer);
	glDeleteFramebuffersEXT(1,&fbo_common_displaying);

	// model space curvature
	glDeleteTextures(1, &texture_model_space_curvature);
	glDeleteRenderbuffersEXT(1, &rbo_model_space_curvature_depth_buffer);
	glDeleteFramebuffersEXT(1, &fbo_model_space_curvature);

	// radial curvature
	glDeleteTextures(1,&texture_radial_curvature);
	glDeleteTextures(1,&texture_primary_curvature_1_direction_and_value);
	glDeleteTextures(1,&texture_primary_curvature_2_direction_and_value);
	glDeleteRenderbuffersEXT(1,&rbo_radial_curvature_depth_buffer);
	glDeleteFramebuffersEXT(1,&fbo_radial_curvature);

	// view dependent curvature
	glDeleteTextures(1,&texture_view_dependent_curvature);
	glDeleteRenderbuffersEXT(1,&rbo_view_dependent_curvature_depth_buffer);
	glDeleteFramebuffersEXT(1,&fbo_view_dependent_curvature); 

	// programs
	cgDestroyProgram(cg_vprogram_radial_curvature);
	cgDestroyProgram(cg_fprogram_radial_curvature);
	cgDestroyProgram(cg_fprogram_view_dependent_curvature);
	cgDestroyProgram(cg_vprogram_view_dependent_curvature);
	cgDestroyContext(cg_context);
}

void MyRender_OUR::initializeGL()
{
    // init GLEW
    if(glewInit() != GLEW_OK)
        int test = 1;// LOG("GLEW initialize error...\n");
    else
        int test = 1;// LOG("GLEW initialized...\n");

    // initialize framebuffer object
    // initialize Cg context and Cg programs
    init_framebuffer_object();
    init_cg_context_and_programs();

    // example lighting
    static const GLfloat light0_color[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
    static const GLfloat light1_color[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
    static const GLfloat light2_color[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
    static const GLfloat light3_color[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
    static const GLfloat light4_color[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
    static const GLfloat light0_pos[4]   = { 999.0f, 999.0f, 999.0f, 0.0f };
    static const GLfloat light1_pos[4]   = { -999.0f, 999.0f, 999.0f, 0.0f };
    static const GLfloat light2_pos[4]   = { 9.0f, 9.0f, 9.0f, 1.0f };
    static const GLfloat light3_pos[4]   = { -999.0f, -999.0f, 999.0f, 0.0f };
    static const GLfloat light4_pos[4]   = { .0f, .0f, -999.0f, 0.0f };

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
    glLightfv(GL_LIGHT0, GL_POSITION, light0_pos);
    glEnable(GL_LIGHT0);

    glLightfv(GL_LIGHT1, GL_DIFFUSE,  light1_color);
    glLightfv(GL_LIGHT1, GL_POSITION, light1_pos);
    glEnable(GL_LIGHT1);

    //glLightfv(GL_LIGHT2, GL_SPECULAR,  light2_color);
    //glLightfv(GL_LIGHT2, GL_POSITION, light2_pos);
    //glEnable(GL_LIGHT2);
    //
    //glLightfv(GL_LIGHT3, GL_DIFFUSE,  light3_color);
    //glLightfv(GL_LIGHT3, GL_POSITION, light3_pos);
    //glEnable(GL_LIGHT3);

    //glLightfv(GL_LIGHT4, GL_DIFFUSE,  light4_color);
    //glLightfv(GL_LIGHT4, GL_POSITION, light4_pos);
    //glEnable(GL_LIGHT4);

    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    // init texture type
    glEnable(GL_TEXTURE_RECTANGLE_ARB);
    glEnable(GL_TEXTURE_2D);
}

void MyRender_OUR::resizeGL(int width, int height)
{
    glViewport(0,0,width,height);
}

void MyRender_OUR::Draw2D()
{

	// 2d condition

		memset(drawing2d.get(),0,3*sizeof(float)*model2d->width*model2d->height);
		// normalize the potential field for rendering
		float maxv=0;
		for (unsigned int i=0; i<model2d->width*model2d->height; i++)
		{
			if (potential2d[i] > maxv)
			{
				maxv = potential2d[i];
			}
		}
		if (maxv != 0)
		{
			for (unsigned int i=0; i<model2d->width*model2d->height; i++)
			{
				drawing2d[i*3] = potential2d[i]/maxv;
			}
		}

		switch(showing_type)
		{
		case SHOWING_MODEL:
			for (unsigned int i=0; i<model2d->width*model2d->height; i++)
			{
				if (model2d->data[i*3 + 0] != 255 ||
					model2d->data[i*3 + 1] != 255 ||
					model2d->data[i*3 + 2] != 255)
					drawing2d[3*i + 0] = 1.0f;
			}
			break;
		case SHOWING_WHOLE_SKELETON:
			if(final_skls2d.size() == 0) break;
			// mix with final skeleton
			for (unsigned int i=0; i<final_skls2d.size(); i++)
			{
				for (unsigned int j=0; j<final_skls2d[i].size(); j++)
				{
					unsigned int x = final_skls2d[i][j][0];
					unsigned int y = final_skls2d[i][j][1];
					drawing2d[3*(x+y*model2d->width)+1] = 1.0f;
				}
			}
			// mix with corresponding model segmentations
			for (unsigned int i=0; i<segmentation2d.size(); i++)
			{
				for (unsigned int j=0; j<segmentation2d[i].size(); j++)
				{
					unsigned int x = segmentation2d[i][j][0];
					unsigned int y = segmentation2d[i][j][1];
					drawing2d[3*(x+y*model2d->width)+2] = 1.0f;
				}
			}
			break;
		case SHOWING_PRIMARY_SKELETON:
			if(primary_skls2d.size() == 0) break;
			for (unsigned int i=1; i<primary_skls2d[showing_primary_skeleton_number_2d].size(); i++)
			{
				unsigned int x = primary_skls2d[showing_primary_skeleton_number_2d][i][0];
				unsigned int y = primary_skls2d[showing_primary_skeleton_number_2d][i][1];
				drawing2d[3*(x+y*model2d->width)+1] = 1.0f;
			}
			break;
		case SHOWING_FINAL_SKELETON:
			if(final_skls2d.size() == 0) break;
			//int test = 1;// LOG(showing_final_skeleton_number_2d<<": \n";
			for (unsigned int i=1; i<final_skls2d[showing_final_skeleton_number_2d].size(); i++)
			{
				unsigned int x = final_skls2d[showing_final_skeleton_number_2d][i][0];
				unsigned int y = final_skls2d[showing_final_skeleton_number_2d][i][1];
				drawing2d[3*(x+y*model2d->width)+0] = 1.0f;
				drawing2d[3*(x+y*model2d->width)+1] = 1.0f;
				drawing2d[3*(x+y*model2d->width)+2] = 0.0f;
				//int test = 1;// LOG(final_skls2d[showing_final_skeleton_number_2d][i];
			}
			for (unsigned int i=1; i<segmentation2d[showing_final_skeleton_number_2d].size(); i++)
			{
				unsigned int x = segmentation2d[showing_final_skeleton_number_2d][i][0];
				unsigned int y = segmentation2d[showing_final_skeleton_number_2d][i][1];
				drawing2d[3*(x+y*model2d->width)+2] = 1.0f;
			}
			//int test = 1;// LOG("\n";
			break;
		}
		// drawing the result
		glDrawPixels(model2d->width, model2d->height, GL_RGB, GL_FLOAT, drawing2d.get());

}


void MyRender_OUR::DrawModel()
{
	if(mesh3d->faces.size() == 0) return;

	float mat_specular_o[] = { 0.0, 0.0, 0.0, 1.0 };
	float low_shininess_o[] = { 0.0 };

	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular_o);
	glMaterialfv(GL_FRONT, GL_SHININESS, low_shininess_o);

	//glColor3f(1.0f, 0.0f, 0.0f);
	//glPointSize(20);
	//glBegin(GL_POINTS);
	//glVertex3f(0, 0, 0);
	//glEnd();


	switch(showing_type)
	{
	case /*VISIBLE_TEST*/SHOWING_MODEL:
	case SHOWING_VIEW_SPHERE_MAP:
		glColor3f(0.5f, 0.6f, 0.7f);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular_o);
		glMaterialfv(GL_FRONT, GL_SHININESS, low_shininess_o);
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);

		glDrawElements(GL_TRIANGLES,(GLsizei)mesh3d->faces.size()*3,GL_UNSIGNED_INT,&mesh3d->faces[0]);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		break;
	case SHOWING_MESH_SALIENCY:
		if (mesh3d->saliency.size() == 0) break;
		glColorPointer(3, GL_FLOAT, 0, &mesh3d->mapped_color[0]);
		glEnableClientState(GL_COLOR_ARRAY);
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);

		glDrawElements(GL_TRIANGLES,(GLsizei)mesh3d->faces.size()*3,GL_UNSIGNED_INT,&mesh3d->faces[0]);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);
		//glBegin(GL_TRIANGLES);
		//glColor3f(1.0,0,0);
		//for (unsigned int i=0; i<mesh3d->faces.size(); i++)
		//{
		//	for (int j=0; j<3; j++)
		//	{
		//		std::vector<unsigned int>::size_type v = mesh3d->faces[i][j];
		//		glNormal3f(
		//			mesh3d->normals[v][0],
		//			mesh3d->normals[v][1],
		//			mesh3d->normals[v][2]);
		//		glColor3f(
		//			mesh3d->mapped_color[3*v + 0],
		//			mesh3d->mapped_color[3*v + 1],
		//			mesh3d->mapped_color[3*v + 2]);
		//		glVertex3f(
		//			mesh3d->vertices[v][0],
		//			mesh3d->vertices[v][1],
		//			mesh3d->vertices[v][2]);
		//	}
		//}
		//glEnd();

		break;
	case SHOWING_SEGMENTATION:
		if (mesh3d->segmentation.size() != mesh3d->faces.size())
			return;
		if(selecting_parts_by_mouse)
		{
			glDrawPixels(fbo_size,fbo_size,
				GL_RGB,GL_UNSIGNED_BYTE,parts_color_for_rendering);
			break;
		}
		glBegin(GL_TRIANGLES);
		for (unsigned int i=0; i<mesh3d->faces.size(); i++)
		{
			int segn = mesh3d->segmentation[i];
			glColor3ubv(seg_colors[segn%(sizeof(seg_colors)/sizeof(seg_colors[0]))]);
			glNormal3f(
				mesh3d->normals[mesh3d->faces[i][0]][0],
				mesh3d->normals[mesh3d->faces[i][0]][1],
				mesh3d->normals[mesh3d->faces[i][0]][2]);
			glVertex3f(
				mesh3d->vertices[mesh3d->faces[i][0]][0],
				mesh3d->vertices[mesh3d->faces[i][0]][1],
				mesh3d->vertices[mesh3d->faces[i][0]][2]);
			glNormal3f(
				mesh3d->normals[mesh3d->faces[i][1]][0],
				mesh3d->normals[mesh3d->faces[i][1]][1],
				mesh3d->normals[mesh3d->faces[i][1]][2]);
			glVertex3f(
				mesh3d->vertices[mesh3d->faces[i][1]][0],
				mesh3d->vertices[mesh3d->faces[i][1]][1],
				mesh3d->vertices[mesh3d->faces[i][1]][2]);
			glNormal3f(
				mesh3d->normals[mesh3d->faces[i][2]][0],
				mesh3d->normals[mesh3d->faces[i][2]][1],
				mesh3d->normals[mesh3d->faces[i][2]][2]);
			glVertex3f(
				mesh3d->vertices[mesh3d->faces[i][2]][0],
				mesh3d->vertices[mesh3d->faces[i][2]][1],
				mesh3d->vertices[mesh3d->faces[i][2]][2]);
		}
		glEnd();
		break;

	case /*SHOWING_MODEL*/VISIBLE_TEST:
		glPushAttrib(GL_ALL_ATTRIB_BITS );
		glDisable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		glDisable(GL_BLEND);

		// 三个颜色通道都使用，每种颜色均分，得到对应的三角形编号
		// 颜色1不使用
		int num_t = mesh3d->faces.size()/3 + 1;

		glBegin(GL_TRIANGLES);
		GLfloat t_color[3];
		for (unsigned int i=0; i<mesh3d->faces.size(); i++)
		{
			//memset(t_color,1.0f,  3*sizeof(GLfloat));
			for(int n = 0; n < 3; n++)
				t_color[n] = 1.0f;
			int s = i/num_t;

			t_color[s] = GLfloat(i-s*num_t)/(num_t);

			glColor3fv(t_color);			
			
			//glPassThrough(i);
			//glLoadName(i);										// Assign Object A Name (ID)
			glVertex3f(
				mesh3d->vertices[mesh3d->faces[i][0]][0],
				mesh3d->vertices[mesh3d->faces[i][0]][1],
				mesh3d->vertices[mesh3d->faces[i][0]][2]);
			glVertex3f(
				mesh3d->vertices[mesh3d->faces[i][1]][0],
				mesh3d->vertices[mesh3d->faces[i][1]][1],
				mesh3d->vertices[mesh3d->faces[i][1]][2]);
			glVertex3f(
				mesh3d->vertices[mesh3d->faces[i][2]][0],
				mesh3d->vertices[mesh3d->faces[i][2]][1],
				mesh3d->vertices[mesh3d->faces[i][2]][2]);
		}
		glEnd();
		glPopAttrib();
		break;
	}
}

void MyRender_OUR::DrawVDCurvature()
{

		if(mesh_viewplane != 0) delete mesh_viewplane;
		mesh_viewplane = GetViewplane();
		mesh_viewplane->compute_curvatures();
		//reset the canvas background
		memset(RGB_mapped_from_curvature,96,
			3*fbo_size*fbo_size*sizeof(unsigned char));
		for (int i=0; i<fbo_size*fbo_size; i++)
		{   // map the curvatures to color space for displaying
			// just for view, so we only handle curvature < 255
			int idx = mesh_viewplane_idx[i];
			if (idx == -1) continue; // not covered

			float mc = (mesh_viewplane->curv1[idx] + mesh_viewplane->curv2[idx])/2.0;
			mc*=1000;

			if (mc > 255)
			{
				RGB_mapped_from_curvature[i*3+0] = 255;
				RGB_mapped_from_curvature[i*3+1] = 0;
				RGB_mapped_from_curvature[i*3+2] = 0;
			}
			else if (mc > 0 && mc <= 255)
			{
				RGB_mapped_from_curvature[i*3+0] = 255;
				RGB_mapped_from_curvature[i*3+1] = 255 - int(mc);
				RGB_mapped_from_curvature[i*3+2] = 255 - int(mc);
			}
			else if (mc < 0 && mc >= -255)
			{
				RGB_mapped_from_curvature[i*3+0] = 255 + int(mc);
				RGB_mapped_from_curvature[i*3+1] = 255;
				RGB_mapped_from_curvature[i*3+2] = 255 + int(mc);
			}
			else if (mc < -255)
			{
				RGB_mapped_from_curvature[i*3+0] = 0;
				RGB_mapped_from_curvature[i*3+1] = 255;
				RGB_mapped_from_curvature[i*3+2] = 0;
			}
		}
		// render to system framebuffer
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
		//glViewport(0,0,fbo_size,fbo_size);
		glClearColor(0,0,0,0);
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		// set to system framebuffer viewport
		glDrawPixels(fbo_size,fbo_size,GL_RGB,GL_UNSIGNED_BYTE,RGB_mapped_from_curvature);

}

void MyRender_OUR::DrawViewSphere()
{
	//////////////////////////////////////////////////////////////////////////
	glEnable(GL_CULL_FACE); // backface culling
	glPointSize(8);
	glDisable(GL_LIGHTING);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	//glColor4f(0.2f, 0.6, 0.1, 0.5);

	//glEnable(GL_POLYGON_OFFSET_FILL);
	//glPolygonOffset(1.0, 1.0);

	//glBegin(GL_TRIANGLES);
	//for(int i = 0; i < vf.size(); i++)
	//{
	//	//glColor4f(viewpoint_candidates.at(vf[i].p1)->color[0],viewpoint_candidates.at(vf[i].p1)->color[1],viewpoint_candidates.at(vf[i].p1)->color[2],0.5f);
	//	glVertex3fv(viewpoint_candidates.at(vf[i].p1)->pos);
	//	//glColor4f(viewpoint_candidates.at(vf[i].p2)->color[0],viewpoint_candidates.at(vf[i].p2)->color[1],viewpoint_candidates.at(vf[i].p2)->color[2],0.5f);
	//	glVertex3fv(viewpoint_candidates.at(vf[i].p2)->pos);
	//	//glColor4f(viewpoint_candidates.at(vf[i].p3)->color[0],viewpoint_candidates.at(vf[i].p3)->color[1],viewpoint_candidates.at(vf[i].p3)->color[2],0.5f);
	//	glVertex3fv(viewpoint_candidates.at(vf[i].p3)->pos);
	//}
	//glEnd();

	//glDisable(GL_CULL_FACE); // backface culling

	glDisable(GL_POLYGON_OFFSET_FILL);
	glDisable(GL_BLEND);

	glColor4f(0.0f, 0.0, 0.0, 0.5);
	glPolygonMode( GL_FRONT_AND_BACK , GL_LINE);
	glBegin(GL_TRIANGLES);
	for(int i = 0; i < vf.size(); i++)
	{
		//glColor3fv(viewpoint_candidates.at(vf[i].p1)->color);
		glVertex3fv(viewpoint_candidates.at(vf[i].p1)->pos);
		//glColor3fv(viewpoint_candidates.at(vf[i].p2)->color);
		glVertex3fv(viewpoint_candidates.at(vf[i].p2)->pos);
		//glColor3fv(viewpoint_candidates.at(vf[i].p3)->color);
		glVertex3fv(viewpoint_candidates.at(vf[i].p3)->pos);
	}
	glEnd();

	glPolygonMode( GL_FRONT , GL_POINT);
	glColor4f(1.0f, 0.0, 0.0, 0.5);
	glBegin(GL_TRIANGLES);
	for(int i = 0; i < vf.size(); i++)
	{
		//glColor4f(viewpoint_candidates.at(vf[i].p1)->color[0],viewpoint_candidates.at(vf[i].p1)->color[1],viewpoint_candidates.at(vf[i].p1)->color[2],0.5f);
		glVertex3fv(viewpoint_candidates.at(vf[i].p1)->pos);
		//glColor4f(viewpoint_candidates.at(vf[i].p2)->color[0],viewpoint_candidates.at(vf[i].p2)->color[1],viewpoint_candidates.at(vf[i].p2)->color[2],0.5f);
		glVertex3fv(viewpoint_candidates.at(vf[i].p2)->pos);
		//glColor4f(viewpoint_candidates.at(vf[i].p3)->color[0],viewpoint_candidates.at(vf[i].p3)->color[1],viewpoint_candidates.at(vf[i].p3)->color[2],0.5f);
		glVertex3fv(viewpoint_candidates.at(vf[i].p3)->pos);
	}
	glEnd();
}

void MyRender_OUR::DrawVoxel()
{

		float trans1[16],rot1[16];  // translate and rotate matrix
		buildTranslateMatrix(
			-model3d->sx/2.0f,
			-model3d->sy/2.0f,
			-model3d->sz/2.0f,
			trans1);
		sgi_trackball_space::build_row_major_rotmatrix(rot1,trackball_quat);
		//set viewpoint according to auto sampling parameter
		buildLookAtMatrix(
			0.0f,0.0f,3.0f*model3d->sz,
			0.0f,0.0f,0.0f,
			0.0f,1.0f,0.0f,
			view_matrix);
		buildOrthoMatrix(
			-model3d->sx/2.0, model3d->sx/2.0,
			-model3d->sy/2.0, model3d->sy/2.0,
			model3d->sz*2.0, model3d->sz*4.0,
			project_matrix);
		multMatrix(model_matrix, rot1,trans1);
		multMatrix(model_view_matrix, view_matrix, model_matrix);
		multMatrix(model_view_project_matrix,project_matrix,model_view_matrix);

		glClearColor(1,1,1,1);
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		glMatrixMode(GL_PROJECTION);
		transposeMatrix(project_matrix);
		glLoadMatrixf(project_matrix);

		glMatrixMode(GL_MODELVIEW);
		transposeMatrix(model_view_matrix);
		glLoadMatrixf(model_view_matrix);

		switch(showing_type)
		{
		case SHOWING_MODEL_VOXELS:
			// drawing the model and depth
			glBegin(GL_POINTS);
			glColor3f(0,0.4f,0.8f);
			for (int i=0; i<model3d->sx*model3d->sy*model3d->sz; i++)
			{
				if (model3d->data[i] != 0)
				{
					glVertex3f(
						(i%(model3d->sx*model3d->sy))%model3d->sx,
						(i%(model3d->sx*model3d->sy))/model3d->sx,
						i/(model3d->sx*model3d->sy));
				}
			}
			glEnd();
			break;
		case SHOWING_WHOLE_SKELETON:
			if(final_skls3d.size() == 0) break;
			glBegin(GL_LINE_STRIP);
			for (unsigned int i=0; i<final_skls3d.size(); i++)
			{
				glColor3f(rand()%100/100.0,rand()%100/100.0,rand()%100/100.0);
				for (unsigned int j=0; j<final_skls3d[i].size(); j++)
				{
					glVertex3f(
						final_skls3d[i][j][0],
						final_skls3d[i][j][1],
						final_skls3d[i][j][2]);
				}
			}
			glEnd();
			break;
		case SHOWING_PRIMARY_SKELETON:
			if(primary_skls3d.size() == 0 ) break;
			glBegin(GL_LINE_STRIP);
			glColor3f(0,1,0.5);
			for (unsigned int i=0; i<primary_skls3d[showing_primary_skeleton_number_3d].size(); i++)
			{
				glVertex3f(
					primary_skls3d[showing_primary_skeleton_number_3d][i][0],
					primary_skls3d[showing_primary_skeleton_number_3d][i][1],
					primary_skls3d[showing_primary_skeleton_number_3d][i][2]);
			}
			glEnd();
			break;
		case SHOWING_FINAL_SKELETON:
			if(mesh_seg3d.size() == 0) break;
			glBegin(GL_TRIANGLES);
			glColor3f(1,0.5,0.5);
			for (unsigned int i=0; i<mesh_seg3d[showing_final_skeleton_number_3d]->faces.size(); i++)
			{
				unsigned int v0 = mesh_seg3d[showing_final_skeleton_number_3d]->faces[i][0];
				unsigned int v1 = mesh_seg3d[showing_final_skeleton_number_3d]->faces[i][1];
				unsigned int v2 = mesh_seg3d[showing_final_skeleton_number_3d]->faces[i][2];

				glVertex3f(
					mesh_seg3d[showing_final_skeleton_number_3d]->vertices[v0][0],
					mesh_seg3d[showing_final_skeleton_number_3d]->vertices[v0][1],
					mesh_seg3d[showing_final_skeleton_number_3d]->vertices[v0][2]);
				glVertex3f(
					mesh_seg3d[showing_final_skeleton_number_3d]->vertices[v1][0],
					mesh_seg3d[showing_final_skeleton_number_3d]->vertices[v1][1],
					mesh_seg3d[showing_final_skeleton_number_3d]->vertices[v1][2]);
				glVertex3f(
					mesh_seg3d[showing_final_skeleton_number_3d]->vertices[v2][0],
					mesh_seg3d[showing_final_skeleton_number_3d]->vertices[v2][1],
					mesh_seg3d[showing_final_skeleton_number_3d]->vertices[v2][2]);
			}
			glEnd();
			break;
		}
}
void MyRender_OUR::paintGL()                                                                                   
{
	glClearColor(1,1,1,1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

    if (!model_loaded) return;

	if (dimension == D2)
		Draw2D();
	else
	{
		switch (showing_type)
		{
    /************************************************************************/
    /*    Mesh Representation                                               */
    /************************************************************************/
		case SHOWING_VIEW_SPHERE_MAP:

		case SHOWING_MODEL:
		case SHOWING_SEGMENTATION:
		case SHOWING_MESH_SALIENCY:
		case VISIBLE_TEST:

			glLoadIdentity();
			//glTranslatef(-mesh3d->bsphere.center[0],
			//	-mesh3d->bsphere.center[1],
			//	-mesh3d->bsphere.center[2]);
			gluLookAt( eye_pos_to_bsphere_center[0],
				eye_pos_to_bsphere_center[1],
				eye_pos_to_bsphere_center[2],
				0,0,0,
				0,1,0);
			float rot1[16]; 
			sgi_trackball_space::build_row_major_rotmatrix(rot1,trackball_quat);
			glMultMatrixf(rot1);

			DrawModel();
			if (showing_type == SHOWING_VIEW_SPHERE_MAP)
			{
				/************************************************************************/
				/*  Viewing sphere map to color                                         */
				/************************************************************************/
				DrawViewSphere();
			}

			break;
    /************************************************************************/
    /*   View Selection Condition                                           */
    /************************************************************************/
		case SHOWING_VIEW_DEPENDENT_CURVATURE:
			DrawVDCurvature();
			break;
    /************************************************************************/
    /*  Voxel Representation                                                */
    /************************************************************************/
		case SHOWING_MODEL_VOXELS:
		case SHOWING_WHOLE_SKELETON:
		case SHOWING_PRIMARY_SKELETON: 
		case SHOWING_FINAL_SKELETON:
			DrawVoxel();
			break;
		}
	}
}

bool MyRender_OUR::load_model_2d(std::wstring fn)
{
    // load fn
    model2d = boost::shared_ptr<PPM>(new PPM());
    model2d->load(fn);
    dimension = D2;
    model_loaded = true;

    // allocate potential and drawing memory
    potential2d = boost::shared_array<float>(new float[model2d->width * model2d->height]);
    memset(potential2d.get(),0,sizeof(float)*model2d->width*model2d->height);
    drawing2d = boost::shared_array<float>(new float[model2d->width*model2d->height*3]);
    memset(drawing2d.get(),0,sizeof(float)*model2d->width*model2d->height);
    primary_skls2d.clear();
    final_skls2d.clear();
    segmentation2d.clear();

    return true;
}

bool MyRender_OUR::compute_potential_2d()
{
    memset(potential2d.get(),0,sizeof(float)*model2d->width*model2d->height);

    // statistic points on the model
    unsigned int N=0;
    for (unsigned int i=0; i<model2d->width*model2d->height; i++)
    {
        if (model2d->data[3*i] != 255 ||
            model2d->data[3*i+1] != 255 ||
            model2d->data[3*i+2] != 255)
        {
            N++;
        }
    }
    // data_s stores the models' point positions
    boost::scoped_array<unsigned int> data_s(new unsigned int[2*N]);
    unsigned int s=0;
    for (unsigned int i=0; i<model2d->width*model2d->height; i++)
    {
        if (model2d->data[3*i] != 255 ||
            model2d->data[3*i+1] != 255 ||
            model2d->data[3*i+2] != 255)
        {
            data_s[2*s] = i%model2d->width;
            data_s[2*s+1] = i/model2d->width;
            s++;
        }
    }

    // compute the potential
    // TODO: this part can be computed by CUDA
    for (unsigned int i=0; i<model2d->width*model2d->height; i++)
    {
        unsigned int x = i%model2d->width;
        unsigned int y = i/model2d->width;

        double p = 0.0; // potential field accumulation

        for(unsigned int j = 0; j<N; j++)
        {
            unsigned int xx = data_s[2*j];
            unsigned int yy = data_s[2*j+1];

            int dx = x - xx;
            int dy = y - yy;

            double dist = sqrt(double(dx*dx + dy*dy));
            if (dist != 0)
            {
                p += pow(1.0/dist, 4.0);
            }
        }
        potential2d[i] = float(p);
    }
    return true;
}


// skeleton extraction procedure
// if current point is the smallest, on the model, on the edge, return false, vise versa true;
__forceinline bool potential_field_trace_2d(
                              PPM *model,
                              float *field,
                              Pos2D &pos,
                              Pos2D &pos_next)
                              
{
    unsigned int width = model->width;
    unsigned int height = model->height;

    // on the edge
    if (pos[0] == 0 ||
        pos[0] == width -1 ||
        pos[1] == 0 ||
        pos[1] == height -1)
    {
        pos_next = pos;
        return false;
    }

    // on the model
    if (model->data[(width*pos[1] + pos[0])*3] != 255 ||
        model->data[(width*pos[1] + pos[0])*3 + 1] != 255 ||
        model->data[(width*pos[1] + pos[0])*3 + 2] != 255)
    {
        pos_next = pos;
        return false;
    }

#define IDX2D(idxx, idxy) (idxx + (idxy)*model->width)

    // find the smallest one as the next one
    // start point
    float v = field[IDX2D(pos[0],pos[1])];
    pos_next = Pos2D(pos[0], pos[1]);

    float delta_v = -999999.0;  // the difference operator
    for (int i=-1; i<2; i++)
    {
        for (int j=-1; j<2; j++)
        {
            float delta_v_temp = v - field[IDX2D(pos[0]+i,pos[1]+j)];
            // compute direction difference equation, the find the gradient
            if (i*i + j*j == 2)
            {
                delta_v_temp *= 1.0f/1.414f;
            }

            if (delta_v_temp > delta_v)
            {
                delta_v = delta_v_temp;
                pos_next = Pos2D(pos[0]+i, pos[1]+j);
            }
        }
    }

    // current position has the smallest value
    if (pos == pos_next)
    {
        return false;
    }
    return true;
#undef IDX2D
}

bool MyRender_OUR::extract_primary_skeleton_2d()
{
    primary_skls2d.clear();

    // statistic seeds, i.e. voxels near the model 
    std::vector<Pos2D> seeds;
    for (unsigned int i=0; i<model2d->width; i++)
    {
        for (unsigned int j=0; j<model2d->height; j++)
        {
            // the voxel on the edge is not considered as a seed, so just ignore it
            // or on the model, deal with in the same way.
            if (i==0 || i==model2d->width-1 || j==0 || j==model2d->height-1 ||
                model2d->data[(i+j*model2d->width)*3 + 0] != 255 ||
                model2d->data[(i+j*model2d->width)*3 + 1] != 255 ||
                model2d->data[(i+j*model2d->width)*3 + 2] != 255)
            {
                // do nothing
            }
            else
            {
                // check all the near voxels
                for (int ii=-1; ii<2; ii++)
                {
                    for (int jj=-1; jj<2; jj++)
                    {
                        unsigned int idx = (i+ii) + (j+jj)*model2d->width;
                        if (model2d->data[idx*3 + 0] != 255 ||
                            model2d->data[idx*3 + 1] != 255 ||
                            model2d->data[idx*3 + 2] != 255)
                        {
                            // if (i,j) has a neighbor on the model, then, add to the seeds set 
                            seeds.push_back(Pos2D(i,j));
                            goto CurrentPosAdded;   // no need to check other neighbors
                        }
                    }
                }// check the nearby voxels
            }
            CurrentPosAdded:;
        }
    }
	// LOG("Seeds number: " + QString::number(seeds.size() ));

    // search start from all the seeds
    // primary skeleton branches not ended on the edge
    for (unsigned int i=0; i<seeds.size(); i++)
    {
        Pos2D seed = seeds[i];
        Pos2D next;
        std::vector<Pos2D> tskl;    // temporal primary skeleton branch

        tskl.push_back(seed);
        while(potential_field_trace_2d(model2d.get(), potential2d.get(), seed, next))
        {
            tskl.push_back(next);
            seed = next;
        }

        // the last skeleton point, to check if it's on the model border
        unsigned int endx = tskl[tskl.size() - 1][0];
        unsigned int endy = tskl[tskl.size() - 1][1];

        if (tskl.size() != 1 &&
            endx != 0 && endx != model2d->width - 1 &&
            endy != 0 && endy != model2d->height -1)
        {
            primary_skls2d.push_back(tskl);
        }
    }

    // reveres all the primary skeletons
    for (unsigned int i=0; i<primary_skls2d.size(); i++)
    {
        std::reverse(primary_skls2d[i].begin(), primary_skls2d[i].end());
    }

	// LOG("Useful primary skeletons number: " + QString::number(primary_skls2d.size()));
    return true;
}

bool MyRender_OUR::build_skeleton_tree_2d()
{
    // clear the root
    for (std::map<Pos2D, Skl2DNode*, Pos2DLess>::iterator i = skls_tree2d.children_pos.begin();
        i != skls_tree2d.children_pos.end(); i++)
    {
        delete  i->second;
    }
    skls_tree2d.children.clear();
    skls_tree2d.children_pos.clear();

    // build the tree
    Skl2DNode *cur_p;
    for (unsigned int i=0; i<primary_skls2d.size(); i++)
    {
        cur_p = &skls_tree2d;
        std::vector<Pos2D>::iterator ins_p = primary_skls2d[i].begin();
        while (ins_p != primary_skls2d[i].end())
        {
            if (cur_p->children.find(*ins_p) == cur_p->children.end())
            {
                Skl2DNode *tsn = new Skl2DNode();
                tsn->pos = *ins_p;

                cur_p->children[*ins_p] = 1;
                cur_p->children_pos[*ins_p] = tsn;
            }
            else
            {
                cur_p->children[*ins_p] += 1;
            }
            cur_p = cur_p->children_pos[*ins_p];
            ins_p++;
        }
    }

    // trim the tree according to a threshold
    // broad first travel the tree and trim the nodes whose weight < WOS
    std::stack<Skl2DNode*> trim_stack;
    trim_stack.push(&skls_tree2d);
    while (!trim_stack.empty())
    {
        Skl2DNode *t = trim_stack.top();
        trim_stack.pop();
        std::vector<Pos2D> to_be_trimed;

        for (std::map<Pos2D, unsigned int, Pos2DLess>::iterator i = t->children.begin();
            i != t->children.end(); i++)
        {
            Pos2D p = i->first;
            unsigned int w = i->second;

            if ( w > WOS_2D)
            {
                trim_stack.push(t->children_pos[p]);
            }
            else
            {
                to_be_trimed.push_back(p);
            }
        }

        for (unsigned int j=0; j<to_be_trimed.size(); j++)
        {
            t->children.erase(to_be_trimed[j]);
            t->children_pos.erase(to_be_trimed[j]);
        }
    }
    return true;
}

bool MyRender_OUR::extract_final_skeleton_2d()
{
    final_skls2d.clear();

    //find the skeletons path
    std::vector<Skl2DNode*> travel_stack; //used as a stack;
    std::vector<Pos2D> temp_skl;
    travel_stack.push_back(&skls_tree2d);

    // find all the path from root to a leaf as a skeleton
    while (!travel_stack.empty())
    {
        Skl2DNode *t = *travel_stack.rbegin(); // get the last one

        if (t->children_pos.size() == 0)
        {
            // it's a leaf, save the skeleton
            temp_skl.clear();
            temp_skl.push_back((*travel_stack.rbegin())->pos);
            travel_stack.pop_back();

            // save the none-branching path on the leaf node
            while(!travel_stack.empty() && (*travel_stack.rbegin())->children_pos.size() == 1)
            {
                temp_skl.push_back((*travel_stack.rbegin())->pos);
                travel_stack.pop_back();
            }
            if (!travel_stack.empty())
            {
                temp_skl.push_back((*travel_stack.rbegin())->pos);
            }

            // delete the root node (0,0,0)
            if (!temp_skl.empty() && (*temp_skl.rbegin()) == Pos2D(0,0))
            {
                temp_skl.pop_back();
            }

            // not empty then save
            // since the 
            if (!temp_skl.empty())
            {
                final_skls2d.push_back(temp_skl);
            }
        }
        else
        {
            bool all_visited = true; // is there any unvisited child of t?
            for (std::map<Pos2D, Skl2DNode*, Pos2DLess>::iterator i = t->children_pos.begin();
                i != t->children_pos.end(); i++)
            {
                Pos2D pos = i->first;
                Skl2DNode *child = i->second;

                if (child->flag == false)
                {
                    // not visited
                    all_visited = false;

                    child->flag = true;
                    travel_stack.push_back(child);

                    break; // visit the children first 
                }
            }
            // all the child has been visited
            if (all_visited)
            {
                temp_skl.clear();
                // the tail node has been saved before
                travel_stack.pop_back();

                // save the none-branching path on the leaf node
                while(!travel_stack.empty() && (*travel_stack.rbegin())->children_pos.size() == 1)
                {
                    temp_skl.push_back((*travel_stack.rbegin())->pos);
                    travel_stack.pop_back();
                }
                if (!travel_stack.empty())
                {
                    temp_skl.push_back((*travel_stack.rbegin())->pos);
                }

                // delete the root node
                if (!temp_skl.empty() && (*temp_skl.rbegin()) == Pos2D(0,0))
                {
                    temp_skl.pop_back();
                }
                // not empty then save
                if (!temp_skl.empty())
                {
                    final_skls2d.push_back(temp_skl);
                }
            }
        }
    }

    skeleton_computed_2d = true;
    int test = 1;// LOG("We totally get " + QString::number(final_skls2d.size()) + " branches!");

#if 0
    for (unsigned int i=0; i<final_skls2d.size(); i++)
    {
        int test = 1;// LOG("Final skeleton "<<i<<" : ");
        for (unsigned int j=0; j<final_skls2d[i].size(); j++)
        {
            int test = 1;// LOG(final_skls2d[i][j];
        }
        int test = 1;// LOG(std::endl;
    }
#endif

    return true;
}

/************************************************************************/
/*     pskl : primary_skl, fskl : final_skl                             */
/*     return the first collision pos of the two skl                    */
/*     pos ref to the start pos of pskl.begin()                         */
/************************************************************************/
__forceinline int skeleton_collision_pos_2d(std::vector<Pos2D> pskl, std::vector<Pos2D> fskl)
{
    std::vector<Pos2D>::iterator fit = fskl.begin();

    while(
        fit != fskl.end() &&
        std::find(pskl.begin(), pskl.end(), *fit) == pskl.end())
    {
        ++fit;
    }

    if (fit == fskl.end())
    {
        return INT_MAX;
    }
    else
    {
        return int(std::find(pskl.begin(), pskl.end(), *fit) - pskl.begin());
    }
}

bool MyRender_OUR::segmentation_by_skeleton_2d()
{
    segmentation2d.clear();

    // the segment the model
    // init
    for (unsigned int i=0; i<final_skls2d.size(); i++)
    {
        segmentation2d.push_back(std::vector<Pos2D>());
    }
    // let primary skeleton starts from seed points, like final_skls...
    for (unsigned int i=0; i<primary_skls2d.size(); i++)
    {
        std::reverse(primary_skls2d[i].begin(), primary_skls2d[i].end());
    }
    // segment
    boost::scoped_array<int> fp(new int[final_skls2d.size()]);
    for (unsigned int i=0; i<primary_skls2d.size(); i++)
    {
        for (unsigned int j=0; j<final_skls2d.size(); j++)
        {
            fp[j] = skeleton_collision_pos_2d(primary_skls2d[i], final_skls2d[j]); 
        }

        // find max fp
        int min_pos = INT_MAX;
        unsigned int skl_no = 0;
        for (unsigned int j=0; j<final_skls2d.size(); j++)
        {
            if (fp[j] < min_pos)
            {
                min_pos = fp[j];
                skl_no = j;
            }
        }
        segmentation2d[skl_no].push_back(*primary_skls2d[i].begin());
    }
    return true;
}

void MyRender_OUR::ResetEyePos()
{
	// 设定一个恰好能观察到模型的位置 [5/1/2012 Han]
    eye_pos_to_bsphere_center = point(0,0, 3*mesh3d->bsphere.r/sin(3.14159*(45/2.0)/180)/*3*mesh3d->bsphere.r*/);
}

bool MyRender_OUR::load_model_3d(std::wstring fn)
{
    model_filename_3d = fn;
    mesh3d = boost::shared_ptr<Mesh>(Mesh::read(fn));
    whole_mesh_buff = mesh3d;
    dimension = D3;
    model_loaded = true;
	
    mesh3d->compute_bsphere();

	mesh3d->align_center();	// 将所有顶点进行偏移，以便使模型中心和原点对齐
	
    //int test = 1;// LOG(mesh3d->bsphere.center[0]<<' '
    //    <<mesh3d->bsphere.center[1]<<' '
    //    <<mesh3d->bsphere.center[2]<<'\n';

    mesh3d->compute_curvatures();
    
    // init for computing view dependent curvature
    init_mesh_vertex_array_pointer();
	
	ResetEyePos();

    model3d = boost::shared_ptr<VolumeData>(new VolumeData());
    model3d->voxelize(mesh3d.get());

    potential3d = boost::shared_array<float>(new float[model3d->sx*model3d->sy*model3d->sz]);
    memset(potential3d.get(), 0, sizeof(float)*model3d->sx*model3d->sy*model3d->sz);

    primary_skls3d.clear();
    final_skls3d.clear();
    segmentation3d.clear();
    reset_trackball();


    //generating viewpoint candidates
    //viewpoint_candidates.clear();
    //viewpoint_candidates.resize(VIEWPOINT_SAMPLES_NUM);
    //generate_viewpoint_candidates();

    //float r = mesh3d->bsphere.r;   //radius of bounding sphere
    //for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    //{
    //    viewpoint_candidates[i][0] = 3*r*meta_viewpoint_candidates[i][0];
    //    viewpoint_candidates[i][1] = 3*r*meta_viewpoint_candidates[i][1];
    //    viewpoint_candidates[i][2] = 3*r*meta_viewpoint_candidates[i][2];
    //}
	// 使用均匀划分的候选视点
	EqualDistributeViews(10, 10, 3*mesh3d->bsphere.r);

    return true;
}

bool MyRender_OUR::load_segmentation(std::wstring fn)
{
    mesh3d->read_segmentation(fn);

    // number of segmentations
    unsigned int nos=0;
    for (unsigned int i=0; i<mesh3d->segmentation.size(); i++)
    {
        if (mesh3d->segmentation[i] > nos)
        {
            nos = mesh3d->segmentation[i];
        }
    }
    return true;
}

bool MyRender_OUR::compute_potential_3d()
{
    memset(potential3d.get(),0,sizeof(float)*model3d->sx*model3d->sy*model3d->sz);

    // statistic points on the model
    unsigned int N=0;
    for (int i=0; i<model3d->sx*model3d->sy*model3d->sz; i++)
    {
        if (model3d->data[i] != 0)
        {
            N++;
        }
    }
    // data_s stores the models' point positions
    boost::scoped_array<unsigned int> data_s(new unsigned int[3*N]);
    unsigned int s=0;
    for (int i=0; i<model3d->sx*model3d->sy*model3d->sz; i++)
    {
        if (model3d->data[i] != 0)
        {
            data_s[3*s] = (i%(model3d->sx*model3d->sy))%model3d->sx;
            data_s[3*s+1] = (i%(model3d->sx*model3d->sy))/model3d->sx;
            data_s[3*s+2] = i/(model3d->sx*model3d->sy);
            s++;
        }
    }

    // compute potential
    for (int i=0; i<model3d->sx*model3d->sy*model3d->sz; i++)
    {
        unsigned int x = (i%(model3d->sx*model3d->sy))%model3d->sx;
        unsigned int y = (i%(model3d->sx*model3d->sy))/model3d->sx;
        unsigned int z = i/(model3d->sx*model3d->sy);

        double p = 0.0; // potential field accumulation

        for(unsigned int j = 0; j<N; j++)
        {
            unsigned int xx = data_s[3*j];
            unsigned int yy = data_s[3*j+1];
            unsigned int zz = data_s[3*j+2];

            int dx = x - xx;
            int dy = y - yy;
            int dz = z - zz;

            double dist = sqrt((double)(dx*dx + dy*dy + dz*dz));
            if (dist != 0)
            {
                p += pow(1.0/dist, 4.0);
            }
        }
        potential3d[i] = float(p);
    }

    return true;
}

__forceinline bool potential_field_trace_3d(
                              VolumeData *model,
                              float *field,
                              Pos3D &pos,
                              Pos3D &pos_next)
                              
{
    unsigned int sx= model->sx;
    unsigned int sy= model->sy;
    unsigned int sz= model->sz;

    // on the edge
    if (pos[0] == 0 ||
        pos[0] == sx-1 ||
        pos[1] == 0 ||
        pos[1] == sy-1 ||
        pos[2] == 0 ||
        pos[2] == sz-1)
    {
        pos_next = pos;
        return false;
    }

    // on the model
    if (model->data[pos[0] + pos[1]*sx + pos[2]*sx*sy] != 0)
    {
        pos_next = pos;
        return false;
    }

#define IDX3D(idxx, idxy, idxz) (idxx + (idxy)*model->sx + (idxz)*model->sx*model->sy)
    // find the smallest one as the next one
    // start point
    float v = field[IDX3D(pos[0],pos[1],pos[2])];
    pos_next = Pos3D(pos[0], pos[1], pos[2]);

    float delta_v = -999999.0;  // the difference operator
    for (int i=-1; i<2; i++)
    {
        for (int j=-1; j<2; j++)
        {
            for (int k=-1; k<2; k++)
            {
                float delta_v_temp = v - field[IDX3D(pos[0]+i,pos[1]+j,pos[2]+k)];
                //int test = 1;// LOG("("<<i<<','<<j<<','<<k<<")"<<delta_v_temp<<std::endl;
                if (i*i + j*j + k*k == 2)
                {
                    delta_v_temp *= 1.0f/1.414f;
                }
                if (i*i + j*j + k*k == 3)
                {
                    delta_v_temp *= 1.0f/1.732f;
                }
                if (delta_v_temp > delta_v)
                {
                    delta_v = delta_v_temp;
                    pos_next = Pos3D(pos[0]+i, pos[1]+j, pos[2]+k);
                }
            }
        }
    }

    // current position has the smallest value
    if (pos == pos_next)
    {
        return false;
    }
    return true;
#undef IDX3D
}

bool MyRender_OUR::extract_primary_skeleton_3d()
{
    primary_skls3d.clear();

    // statistic seeds, i.e. voxels near the model 
    std::vector<Pos3D> seeds;
    for (int i=0; i<model3d->sx; i++)
    {
        for (int j=0; j<model3d->sy; j++)
        {
            for (int k=0; k<model3d->sz; k++)
            {
                // the voxel on the edge is not considered as a seed, so just ignore it
                // or on the model, deal with in the same way.
                if (i==0 || i==model3d->sx-1 ||
                    j==0 || j==model3d->sy-1 ||
                    k==0 || k==model3d->sz-1 ||
                    model3d->data[i + j*model3d->sx + k*model3d->sx*model3d->sy] != 0)
                {
                    // do nothing
                }
                else
                {
                    // check all the near voxels
                    for (int ii=-1; ii<2; ii++)
                    {
                        for (int jj=-1; jj<2; jj++)
                        {
                            for (int kk=-1; kk<2; kk++)
                            {
                                unsigned int idx = (i+ii) + (j+jj)*model3d->sx + (k+kk)*model3d->sx*model3d->sy;
                                if (model3d->data[idx] != 0)
                                {
                                    // if (i,j) has a neighbor on the model, then, add to the seeds set 
                                    seeds.push_back(Pos3D(i,j,k));
                                    goto CurrentPosAdded;   // no need to check other neighbors
                                }
                            }
                        }
                    }// check the nearby voxels
                }
                CurrentPosAdded:;
            }
        }
    }
    // LOG("Seeds number: "+QString::number(seeds.size()));

    // search start from all the seeds
    // primary skeleton branches not ended on the edge
    for (unsigned int i=0; i<seeds.size(); i++)
    //for (unsigned int i=593; i<594; i++)
    {
        Pos3D seed = seeds[i];
        Pos3D next;
        std::vector<Pos3D> tskl;    // temporal primary skeleton branch

        tskl.push_back(seed);
        while(potential_field_trace_3d(model3d.get(), potential3d.get(), seed, next))
        {
            tskl.push_back(next);
            seed = next;
        }


        // the last skeleton point, to check if it's on the model border
        unsigned int endx = tskl[tskl.size() - 1][0];
        unsigned int endy = tskl[tskl.size() - 1][1];
        unsigned int endz = tskl[tskl.size() - 1][2];

        if (tskl.size() != 1 &&
            endx != 0 && endx != model3d->sx - 1 &&
            endy != 0 && endy != model3d->sy - 1 &&
            endz != 0 && endz != model3d->sz - 1)
        {
            primary_skls3d.push_back(tskl);
        }
    }

    // reveres all the primary skeletons
    for (unsigned int i=0; i<primary_skls3d.size(); i++)
    {
        std::reverse(primary_skls3d[i].begin(), primary_skls3d[i].end());
    }

    // LOG("Useful primary skeletons number: "+QString::number(primary_skls3d.size()));
    return true;
}

bool MyRender_OUR::build_skeleton_tree_3d()
{
    // clear the root
    for (std::map<Pos3D, Skl3DNode*, Pos3DLess>::iterator i = skls_tree3d.children_pos.begin();
        i != skls_tree3d.children_pos.end(); i++)
    {
        delete  i->second;
    }
    skls_tree3d.children.clear();
    skls_tree3d.children_pos.clear();

    // build the tree
    Skl3DNode *cur_p;
    for (unsigned int i=0; i<primary_skls3d.size(); i++)
    {
        cur_p = &skls_tree3d;
        std::vector<Pos3D>::iterator ins_p = primary_skls3d[i].begin();
        while (ins_p != primary_skls3d[i].end())
        {
            if (cur_p->children.find(*ins_p) == cur_p->children.end())
            {
                Skl3DNode *tsn = new Skl3DNode();
                tsn->pos = *ins_p;

                cur_p->children[*ins_p] = 1;
                cur_p->children_pos[*ins_p] = tsn;
            }
            else
            {
                cur_p->children[*ins_p] += 1;
            }
            cur_p = cur_p->children_pos[*ins_p];
            ins_p++;
        }
    }

    // child of rooter
    //int test = 1;// LOG("Root:";
    //for (std::map<Pos3D, unsigned int, Pos3DLess>::iterator i = skls_tree_root.children.begin();
    //    i != skls_tree_root.children.end(); i++)
    //{
    //    int test = 1;// LOG(i->first<<std::endl;
    //}

    // trim the tree according to a threshold
    // broad first travel the tree and trim the nodes whose weight < WOS
    std::stack<Skl3DNode*> trim_stack;
    trim_stack.push(&skls_tree3d);
    while (!trim_stack.empty())
    {
        Skl3DNode *t = trim_stack.top();
        trim_stack.pop();
        std::vector<Pos3D> to_be_trimed;

        for (std::map<Pos3D, unsigned int, Pos3DLess>::iterator i = t->children.begin();
            i != t->children.end(); i++)
        {
            Pos3D p = i->first;
            unsigned int w = i->second;

            if ( w > WOS_3D)
            {
                trim_stack.push(t->children_pos[p]);
            }
            else
            {
                to_be_trimed.push_back(p);
            }
        }

        for (unsigned int j=0; j<to_be_trimed.size(); j++)
        {
            t->children.erase(to_be_trimed[j]);
            t->children_pos.erase(to_be_trimed[j]);
        }
    }
    return true;
}

bool MyRender_OUR::extract_final_skeleton_3d()
{
    final_skls3d.clear();

    //find the skeletons path
    std::vector<Skl3DNode*> travel_stack; //used as a stack;
    std::vector<Pos3D> temp_skl;
    travel_stack.push_back(&skls_tree3d);

    // find all the path from root to a leaf as a skeleton
    while (!travel_stack.empty())
    {
        Skl3DNode *t = *travel_stack.rbegin(); // get the last one

        if (t->children_pos.size() == 0)
        {
            // it's a leaf, save the skeleton
            temp_skl.clear();
            temp_skl.push_back((*travel_stack.rbegin())->pos);
            travel_stack.pop_back();

            // save the none-branching path on the leaf node
            while(!travel_stack.empty() && (*travel_stack.rbegin())->children_pos.size() == 1)
            {
                temp_skl.push_back((*travel_stack.rbegin())->pos);
                travel_stack.pop_back();
            }
            if (!travel_stack.empty())
            {
                temp_skl.push_back((*travel_stack.rbegin())->pos);
            }

            // delete the root node (0,0,0)
            if (!temp_skl.empty() && (*temp_skl.rbegin()) == Pos3D(0,0,0))
            {
                temp_skl.pop_back();
            }
            // not empty then save
            if (!temp_skl.empty())
            {
                final_skls3d.push_back(temp_skl);
            }
        }
        else
        {
            bool all_visited = true; // is there any unvisited child of t?
            for (std::map<Pos3D, Skl3DNode*, Pos3DLess>::iterator i = t->children_pos.begin();
                i != t->children_pos.end(); i++)
            {
                Pos3D pos = i->first;
                Skl3DNode *child = i->second;

                if (child->flag == false)
                {
                    // not visited
                    all_visited = false;

                    child->flag = true;
                    travel_stack.push_back(child);

                    break; // visit the children first 
                }
            }
            // all the child has been visited
            if (all_visited)
            {
                temp_skl.clear();
                // the tail node has been saved before
                travel_stack.pop_back();

                // save the none-branching path on the leaf node
                while(!travel_stack.empty() && (*travel_stack.rbegin())->children_pos.size() == 1)
                {
                    temp_skl.push_back((*travel_stack.rbegin())->pos);
                    travel_stack.pop_back();
                }
                if (!travel_stack.empty())
                {
                    temp_skl.push_back((*travel_stack.rbegin())->pos);
                }

                // delete the root node
                if (!temp_skl.empty() && (*temp_skl.rbegin()) == Pos3D(0,0,0))
                {
                    temp_skl.pop_back();
                }
                // not empty then save
                if (!temp_skl.empty())
                {
                    final_skls3d.push_back(temp_skl);
                }
            }
        }
    }

    int test = 1;// LOG("We totally get "+QString::number(final_skls3d.size())+" branches!");
    return true;
}

/************************************************************************/
/*     pskl : primary_skl, fskl : final_skl                             */
/*     return the first collision pos of the two skl                    */
/*     pos ref to the start pos of pskl.begin()                         */
/************************************************************************/
__forceinline int skeleton_collision_pos_3d(std::vector<Pos3D> pskl, std::vector<Pos3D> fskl)
{
    std::vector<Pos3D>::iterator fit = fskl.begin();

    while(
        fit != fskl.end() &&
        std::find(pskl.begin(), pskl.end(), *fit) == pskl.end())
    {
        ++fit;
    }

    if (fit == fskl.end())
    {
        return INT_MAX;
    }
    else
    {
        return int(std::find(pskl.begin(), pskl.end(), *fit) - pskl.begin());
    }
}

bool MyRender_OUR::segmentation_by_skeleton_3d()
{
    segmentation3d.clear();

    // the segment the model
    for (unsigned int i=0; i<final_skls3d.size(); i++)
    {
        segmentation3d.push_back(std::vector<Pos3D>());
    }
    // let primary skeleton starts from seed points, like final_skls...
    for (unsigned int i=0; i<primary_skls3d.size(); i++)
    {
        std::reverse(primary_skls3d[i].begin(), primary_skls3d[i].end());
    }
    // segment
    boost::scoped_array<int> fp(new int[final_skls3d.size()]);
    for (unsigned int i=0; i<primary_skls3d.size(); i++)
    {
        for (unsigned int j=0; j<final_skls3d.size(); j++)
        {
            fp[j] = skeleton_collision_pos_3d(primary_skls3d[i], final_skls3d[j]); 
        }

        // find max fp
        int min_pos = INT_MAX;
        unsigned int skl_no = 0;
        for (unsigned int j=0; j<final_skls3d.size(); j++)
        {
            if (fp[j] < min_pos)
            {
                min_pos = fp[j];
                skl_no = j;
            }
        }
        if (min_pos != INT_MAX) //  a single branch from the root
        {
            segmentation3d[skl_no].push_back(*primary_skls3d[i].begin());
        }
    }

    // map to model voxels
    model_seg3d.resize(segmentation3d.size());
    for (unsigned int i=0; i<segmentation3d.size(); i++)
    {
        for (unsigned int j=0; j<segmentation3d[i].size(); j++)
        {
            // check is there any voxel adjacent segmentation3d[i][j] is on the model and are not visited
            // if so, add to the corresponding model segmentation part
            // and set to 0, mark as visited
            for (int x=-1; x<2; x++)
            {
                for (int y=-1; y<2; y++)
                {
                    for (int z=-1; z<2; z++)
                    {
                        int idx =
                            segmentation3d[i][j][0] + x +
                            segmentation3d[i][j][1]*model3d->sx + y +
                            segmentation3d[i][j][2]*model3d->sx*model3d->sy + z;
                        if (model3d->data[idx] > 0)
                        {
                            model_seg3d[i].push_back(Pos3D(
                                segmentation3d[i][j][0] + x,
                                segmentation3d[i][j][1] + y,
                                segmentation3d[i][j][2] + z));
                            model3d->data[idx] = 0;
                        }
                    }
                }
            }// int x, y, z -1,0,1
        }
    }

    int c=0;
    for (int i=0; i<model3d->sx*model3d->sy*model3d->sz; i++)
    {
        if (model3d->data[i] > 0)
        {
            c++;
        }
    }
    int test = 1;// LOG("\nUnhandled voxels : "+QString::number(c));
    return true;
}

bool MyRender_OUR::map_to_mesh_segmentation()
{
    if (mesh_seg3d.size() > 0)
    {
        for (unsigned int i=0; i<mesh_seg3d.size(); i++)
        {
            delete mesh_seg3d[i];
        }
    }
    mesh_seg3d.resize(model_seg3d.size());
    for (unsigned int i=0; i<model_seg3d.size(); i++)
    {
        mesh_seg3d[i] = new Mesh();
        mesh_seg3d[i]->vertices = mesh3d->vertices;
    }

    std::vector<char> faces;
    faces.resize(mesh3d->faces.size(),0);
    float lx = mesh3d->bsphere.center[0] - mesh3d->bsphere.r; // low bounder of x
    float ux = mesh3d->bsphere.center[0] + mesh3d->bsphere.r; // up bounder of x
    float dx = ux - lx; // length of the bounding box along x axis
    float ly = mesh3d->bsphere.center[1] - mesh3d->bsphere.r;
    float uy = mesh3d->bsphere.center[1] + mesh3d->bsphere.r;
    float dy = uy - ly;
    float lz = mesh3d->bsphere.center[2] - mesh3d->bsphere.r;
    float uz = mesh3d->bsphere.center[2] + mesh3d->bsphere.r;
    float dz = uz - lz;
    float vs = (dx+dy+dz)/VolumeData::GRANULARITY;  // voxel size

    for (unsigned int i=0; i<mesh3d->faces.size(); i++)
    {
        // group each face to certain mesh_seg3d[i];
        // we use a simple criterion, mesh3d->faces[i][0] \belong voxel \belong mesh_seg3d[x]
        // then we set mesh3d->faces[i] \belong mesh_seg3d[x]
        unsigned int idx = mesh3d->faces[i][0];
        unsigned int ix = unsigned int((mesh3d->vertices[idx][0] - lx) / vs);
        unsigned int iy = unsigned int((mesh3d->vertices[idx][1] - ly) / vs);
        unsigned int iz = unsigned int((mesh3d->vertices[idx][2] - lz) / vs);           
        Pos3D v(ix,iy,iz);
        // find v in which model_seg3d[x]
        for (unsigned int j=0; j<model_seg3d.size(); j++)
        {
            std::vector<Pos3D>::iterator pos = find(model_seg3d[j].begin(),model_seg3d[j].end(),v);
            if (pos != model_seg3d[j].end())
            {
                mesh_seg3d[j]->faces.push_back(mesh3d->faces[i]);
                break;
            }
        }
    }
    return true;
}

void MyRender_OUR::add_trackball_quat(float x, float y, float xx, float yy)
{
    if (!selecting_parts_by_mouse)
    {
        //add a new rotation operation to this trackball
        float quat[4];
        sgi_trackball_space::trackball(quat,x,y,xx,yy);
        sgi_trackball_space::add_quats(quat, trackball_quat, trackball_quat);
    }
}

void MyRender_OUR::reset_trackball()
{
    sgi_trackball_space::trackball(trackball_quat, 0, 0, 0, 0);
}


bool MyRender_OUR::save_potential_3d(std::wstring fn)
{
    std::ofstream f(fn.c_str(), std::ios::binary|std::ios::out);

    // length of model filename, and save it,
    int s = model_filename_3d.size()*sizeof(wchar_t);
    f.write((char*)&s, sizeof(int));
    // save the corresponding file name
    f.write((char*)model_filename_3d.c_str(), s);
    // save the potential
    f.write((char*)potential3d.get(), sizeof(float)*model3d->sx*model3d->sy*model3d->sz);

    f.close();
    return true;
}

bool MyRender_OUR::load_potential_3d(std::wstring fn)
{
    std::fstream f(fn.c_str(), std::ios::binary|std::ios::in);

    // read length of model filename
    int s=0;
    f.read((char*)&s, sizeof(int));
    // read model filename
    boost::scoped_array<char> mn(new char[s]);
    f.read(mn.get(), s);
    std::wstring wmn = std::wstring((wchar_t*)mn.get(), s/sizeof(wchar_t));
    // read model file
    load_model_3d(wmn);
    // read potential data
    potential3d = boost::shared_array<float>(new float[model3d->sx*model3d->sy*model3d->sz]);
    f.read((char*)potential3d.get(),sizeof(float)*model3d->sx*model3d->sy*model3d->sz);
    f.close();
    return true;
}


bool MyRender_OUR::compute_mesh_saliency_3d()
{
    mesh3d->compute_mesh_saliency(true);


	// test:save some some rendered saliency mesh [3/19/2012 Han]
	// 保存得到视点对应的渲染图像 [2/24/2012 Han]
	//time_t tm;
	//srand((unsigned)time(&tm));  //set seed;
	//Image_Type ba = showing_type;
	//showing_type = SHOWING_MESH_SALIENCY;
	//wchar_t fn[128];
	//for ( int i = 0; i < 4; i++)
	//{
	//	swprintf( fn,   L"%d_%d.ppm ", mesh3d->faces.size(),rand());
	//	eye_pos_to_bsphere_center = viewpoint_candidates[i*viewpoint_candidates.size()/5];
	//	PPM *ppmWriter = new PPM;
	//	ppmWriter->width = fbo_size;
	//	ppmWriter->height = fbo_size;
	//	ppmWriter->version = "P6";
	//	ppmWriter->data = new unsigned char [fbo_size*fbo_size*3];
	//	paintGL();
	//	//	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	//	glReadPixels(0, 0, fbo_size, fbo_size, GL_RGB, GL_UNSIGNED_BYTE, ppmWriter->data);
	//	ppmWriter->save(fn);

	//	delete ppmWriter;
	//}
	//showing_type = ba;

    return true;
}

bool MyRender_OUR::init_mesh_vertex_array_pointer()
{
    // no model has been loaded
    if (!model_loaded) return false;

    if (mesh3d->vertices.size() != 0)
    {
        glVertexPointer(3,GL_FLOAT,0,&mesh3d->vertices[0]);

        //this used as texture coordinate to get principal curvature
        mesh3d->vertices_sequence.resize(mesh3d->vertices.size());
        for (vector<point>::size_type i=0; i<mesh3d->vertices.size(); i++)
        {
            mesh3d->vertices_sequence[i] = i;
        }
        glTexCoordPointer(1,GL_INT,0,&mesh3d->vertices_sequence[0]);

		
    }
    if (mesh3d->normals.size() != 0)
    {
        glNormalPointer(GL_FLOAT,0,&mesh3d->normals[0]);
    }
    return true;
}

bool MyRender_OUR::init_mesh_curvature_texture()
{
    if (!model_loaded)
    {
        int test = 1;// LOG("mesh is empty, no curvature texture initialized.\n");
        return false;
    }

    mesh3d->compute_curvatures();
    int width = (int)floor(sqrt(float(mesh3d->vertices.size())))+1;
    width_of_primary_curvature_texture = width;
    float *pd = new float[4*width*width];
    if (!pd)
    {
        int test = 1;// LOG("Not enough memory, make more money, and upgrade your machine!\n");
        return false;
    }

    // use GL_ABGR_EXT format to upload to GL_RGBA_FLOAT32_ATI as vertex texture
    vector<point>::size_type N = vector<point>::size_type(width*width);
    // upload principal curvature 1
    for (vector<point>::size_type i=0; i<mesh3d->vertices.size(); i++)
    {
        pd[i*4 + 0] = mesh3d->curv1[i];
        pd[i*4 + 1] = mesh3d->pdir1[i][2];
        pd[i*4 + 2] = mesh3d->pdir1[i][1];
        pd[i*4 + 3] = mesh3d->pdir1[i][0];
    }
    for (vector<point>::size_type i=mesh3d->vertices.size(); i<N; i++)
    {
        pd[i*4 + 0] =
            pd[i*4 + 1] =
            pd[i*4 + 2] =
            pd[i*4 + 3] = 0;
    }
    glGenTextures(1,&texture_primary_curvature_1_direction_and_value);
    glBindTexture(GL_TEXTURE_2D,texture_primary_curvature_1_direction_and_value);
    glTexImage2D(
        GL_TEXTURE_2D,0,
        GL_RGBA_FLOAT32_ATI,
        width,width,0,
        GL_ABGR_EXT,
        GL_FLOAT,
        pd);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D,0);

    // upload principal curvature 2
    for (vector<point>::size_type i=0; i<mesh3d->vertices.size(); i++)
    {
        pd[i*4 + 0] = mesh3d->curv2[i];
        pd[i*4 + 1] = mesh3d->pdir2[i][2];
        pd[i*4 + 2] = mesh3d->pdir2[i][1];
        pd[i*4 + 3] = mesh3d->pdir2[i][0];
    }
    for (vector<point>::size_type i=mesh3d->vertices.size(); i<N; i++)
    {
        pd[i*4 + 0] =
            pd[i*4 + 1] =
            pd[i*4 + 2] =
            pd[i*4 + 3] = 0;
    }
    glGenTextures(1,&texture_primary_curvature_2_direction_and_value);
    glBindTexture(GL_TEXTURE_2D,texture_primary_curvature_2_direction_and_value);
    glTexImage2D(
        GL_TEXTURE_2D,0,
        GL_RGBA_FLOAT32_ATI,
        width,width,0,
        GL_ABGR_EXT,
        GL_FLOAT,
        pd);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D,0);

    // release pd and return
    delete []pd;
    return true;
    return true;
}
void MyRender_OUR::init_framebuffer_object()
{
    // LOG("Init FrameBuffer and Check Status\n");
    /************************************************************************/
    /*       FBO AND BUFFER FOR VIEW DEPENDENT CURVATURE COMPUTING          */
    /************************************************************************/

    // generate fbo and textures for view dependent curvature
    glGenFramebuffersEXT(1, &fbo_view_dependent_curvature);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo_view_dependent_curvature);

    // attach 3 channel float point buffer to this fbo
    glGenTextures(1,&texture_view_dependent_curvature);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB,texture_view_dependent_curvature);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(
        GL_TEXTURE_RECTANGLE_ARB,
        0,
        GL_FLOAT_RGB32_NV,
        fbo_size,
        fbo_size,
        0,
        GL_RGB,
        GL_FLOAT,
        NULL);
    glFramebufferTexture2DEXT(
        GL_FRAMEBUFFER_EXT,
        GL_COLOR_ATTACHMENT0_EXT,
        GL_TEXTURE_RECTANGLE_ARB,
        texture_view_dependent_curvature,
        0);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB,0);

    // attach depth-buffer rbo
    glGenRenderbuffersEXT(1,&rbo_view_dependent_curvature_depth_buffer);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, rbo_view_dependent_curvature_depth_buffer);
    glRenderbufferStorageEXT(
        GL_RENDERBUFFER_EXT,
        GL_DEPTH_COMPONENT24,
        fbo_size,
        fbo_size);
    glFramebufferRenderbufferEXT(
        GL_FRAMEBUFFER_EXT,
        GL_DEPTH_ATTACHMENT_EXT,
        GL_RENDERBUFFER_EXT,
        rbo_view_dependent_curvature_depth_buffer);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT,0);

    //check_frame_buffer_object_status("View Dependent Curvature");
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);

   // LOG("End.\n\n");
}
void MyRender_OUR::init_cg_context_and_programs()
{
    // LOG("Load Cg programs and Check status\n");
    /************************************************************************/
    /*           Create Cg Context                                          */
    /************************************************************************/
    cg_context = cgCreateContext();

    /************************************************************************/
    /*           Get latest Cg Profile                                      */
    /************************************************************************/
    cg_vprofile = cgGLGetLatestProfile(CG_GL_VERTEX);
    assert(cg_vprofile != CG_PROFILE_UNKNOWN);
    cgGLSetOptimalOptions(cg_vprofile);

    cg_fprofile = cgGLGetLatestProfile(CG_GL_FRAGMENT);
    assert(cg_fprofile != CG_PROFILE_UNKNOWN);
    cgGLSetOptimalOptions(cg_fprofile);

    /************************************************************************/
    /*           Load Cg Programs For View Dependent Curvature Computing    */
    /************************************************************************/

    // vertex program for view dependent curvature
    cg_vprogram_view_dependent_curvature = cgCreateProgramFromFile(
        cg_context,
        CG_SOURCE,
        "viewdependentcurvature_v.cg",
        cg_vprofile,
        NULL,NULL);
    std::cout
        <<"Compile and Load View Dependent Curvature Computing Vertex Program\nLAST LISTING----"
        <<(cgGetLastListing(cg_context)==NULL?"None":cgGetLastListing(cg_context))
        <<"----\n";
    /*
    std::cout
        <<"---- PROGRAM BEGIN ----\n"
        <<string(cgGetProgramString(cg_vprogram_view_dependent_curvature, CG_COMPILED_PROGRAM))
        <<"---- PROGRAM END ----\n\n";
    //*/
    cgGLLoadProgram(cg_vprogram_view_dependent_curvature);
    cgparam_model_view_project_matrix_view_dependent_curvature =
        cgGetNamedParameter(cg_vprogram_view_dependent_curvature, "model_view_project");

    // fragment program for view dependent curvature
    cg_fprogram_view_dependent_curvature =
        cgCreateProgramFromFile(
        cg_context,
        CG_SOURCE,
        "viewdependentcurvature_f.cg",
        cg_fprofile,
        NULL,NULL);
    std::cout
        <<"Compile and Load View Dependent Curvature Computing Fragment Program\nLAST LISTING----"
        <<(cgGetLastListing(cg_context)==NULL?"None":cgGetLastListing(cg_context))
        <<"----\n";
    /*
    std::cout
        <<"---- PROGRAM BEGIN ----\n"
        <<cgGetProgramString(cg_fprogram_view_dependent_curvature, CG_COMPILED_PROGRAM)
        <<"---- PROGRAM END ----\n\n";
    //*/
    cgGLLoadProgram(cg_fprogram_view_dependent_curvature);

    // LOG("End\n\n");

}

void MyRender_OUR::build_segmentation(std::vector<int> &parts)
{
    if (mesh3d == 0)
    {
        int test = 1;// LOG("load model first...");
        return;
    }
    if (mesh3d->segmentation.size() == 0)
    {
        int test = 1;// LOG("load segmentation first...");
        return;
    }

    vector<Mesh::Face> segmentation_parts_faces; // considered faces in view selection, if 0, the whole model
    vector<int> segmentation_parts; // considered parts.
    std::vector<unsigned int> segmentation_parts_segno;

    segmentation_parts = parts;
    for (unsigned int i=0; i<mesh3d->faces.size(); i++)
    {
        if (find(parts.begin(),parts.end(),mesh3d->segmentation[i]) != parts.end())
        {
            segmentation_parts_faces.push_back(mesh3d->faces[i]);
            segmentation_parts_segno.push_back(mesh3d->segmentation[i]);
        }
    }

    // build new mesh according to segmentations;
    // no parts selected, reset to whole model
    if (parts.size() == 0)
    {
        mesh3d = whole_mesh_buff;
    }
    // segments the model
    else
    {
        mesh3d = boost::shared_ptr<Mesh>(new Mesh());
        // copy data from whole_mesh_buff to mesh3d;
        mesh3d->vertices = whole_mesh_buff->vertices;
        mesh3d->normals = whole_mesh_buff->normals;
        mesh3d->faces = segmentation_parts_faces;
        mesh3d->segmentation = segmentation_parts_segno;
        mesh3d->curv1 = whole_mesh_buff->curv1;
        mesh3d->curv2 = whole_mesh_buff->curv2;

        // removed dumb points, set they position to the first vertex of faces[0]
        std::vector<int> sts;
        sts.resize(mesh3d->vertices.size(),0);
        for (unsigned int i=0; i<mesh3d->faces.size(); i++)
        {
            sts[mesh3d->faces[i][0]]++;
            sts[mesh3d->faces[i][1]]++;
            sts[mesh3d->faces[i][2]]++;
        }
        unsigned int dumb = mesh3d->faces[0][0];
        for (unsigned int i=0; i<sts.size(); i++)
        {
            if (sts[i] == 0)
            {
                mesh3d->vertices[i] = mesh3d->vertices[dumb];
            }
        }
    }
    // recompute the data and initiation
    mesh3d->compute_bsphere();
    init_mesh_vertex_array_pointer();
    eye_pos_to_bsphere_center = point(0,0,3*mesh3d->bsphere.r);
    reset_trackball();
}

typedef struct  
{
	Mesh *viewPlane;
	float *standard_deviation;
	float *shannon_entropy;
	int number_of_histogram_intervals;
} ThreadViewPlane;

float compute_entropy_han(Mesh *viewplane,int number_of_histogram_intervals, float &current_standard_deviation)
{
	if (viewplane == NULL)
	{
		return -1;
	}
    // curvatures for computing entropy
    std::vector<float> curvatures;
/*
    // use 'as' for abbreviation
    int as = adaptive_box_size;
    // statistic of curvature
    for (int j=0; j<fbo_size; j+=as)
    {
        for (int i=0; i<fbo_size; i+=as)
        {
            // at box (i,j)->(i+as,j+as)
            float mc=0;
            for (int m=0; m<as; m++)
            {
                for (int n=0; n<as; n++)
                {
                    // add mean curvature at (i+n, j+m)
                    mc  +=  mesh_viewplane->curv1[(j+m)*fbo_size+i+n]
             +mesh_viewplane->curv2[(j+m)*fbo_size+i+n];
                }
            }
            mc/=2; //get real mean curvature
            mc/=as*as; // get average mean curvature over an adaptive box
            curvatures.push_back(mc);
        }
    }
//*/
    for (unsigned int i=0; i<viewplane->vertices.size(); i++)
    {
        curvatures.push_back((viewplane->curv1[i]+viewplane->curv2[i])/2.0);
    }
#if 0
	FILE *fp =fopen("parts_hist.txt","w");
	fprintf(fp,"%d\n",fbo_size);
	for (int i=0;i<fbo_size;i++)
	{
		for (int j=0;j<fbo_size;j++)
		{
			int idx = mesh_viewplane_idx[j*fbo_size+i];
			if (idx!=-1) fprintf(fp,"%d %d %d %f\n",i,j,idx,curvatures[idx]);
			else fprintf(fp,"%d %d %d\n",i,j,idx);
		}
	}
	fclose(fp);
#endif
    float maxc=0,minc=0;
    for (unsigned int i=0; i<curvatures.size(); i++)
    {
        if (curvatures[i] > maxc)
        {
            maxc = curvatures[i];
        }
        if (curvatures[i] < minc)
        {
            minc = curvatures[i];
        }
    }
    //int test = 1;// LOG("max mc: "<<maxc<<std::endl;
    //int test = 1;// LOG("min mc: "<<minc<<std::endl;

    const unsigned int HISTO_WIDTH = number_of_histogram_intervals;
    std::vector<unsigned int> histogram;
    float hs = (maxc - minc)/HISTO_WIDTH; //histogram step width
    histogram.resize(HISTO_WIDTH,0);
    for (unsigned int i=0; i<curvatures.size(); i++)
    {
        unsigned int idx = int((curvatures[i] - minc)/hs);
        if (idx < histogram.size())
        {
            histogram[idx] += 1;
        }
    }

    // consider all pixels in the image
    histogram[int(-minc)/hs] += MyRender_OUR::fbo_size*MyRender_OUR::fbo_size - viewplane->vertices.size();

#if 1 // for output curvature kinds
        std::ofstream of("hist.txt");
        for(unsigned int i=0; i<HISTO_WIDTH; i++)
        {
            of<<histogram[i]<<'\n';
        }
        of.close();
#endif
    // computing standard deviation
    float avg_mc=0; // average mean curvature,
    for (vector<int>::size_type i=0; i<histogram.size(); i++)
    {
        avg_mc += histogram[i];
    }
    avg_mc/=histogram.size();
    float xigma_mc = 0;
    for (vector<int>::size_type i=0; i<histogram.size(); i++)
    {
        xigma_mc += (avg_mc - histogram[i])*(avg_mc - histogram[i]);
    }
    xigma_mc /= histogram.size() - 1;
    current_standard_deviation = sqrt(xigma_mc);

    // computing entropy 
    float N=0;
    for (vector<int>::size_type i=0; i<histogram.size(); i++)
    {
        N += histogram[i];
    }
    float pi,E=0;
    for (vector<int>::size_type i=0; i<histogram.size(); i++)
    {
        if (histogram[i] != 0)
        {
            pi = histogram[i]/N;
            E += -pi*log(pi)/log(2.0f);
        }
    }

    // LOG("Shannon Entropy: "+QString::number(E));
    // LOG("Standard Deviation: "+QString::number(current_standard_deviation));
    return E;
}
float MyRender_OUR::compute_entropy()
{
    if (!mesh3d)
    {
        int test = 1;// LOG("plz load a mesh first....\n");
        return -1;
    }
    if (showing_type != SHOWING_VIEW_DEPENDENT_CURVATURE)
    {
        int test = 1;// LOG("plz select a correct type first\n");
        return -1;
    }

  return compute_entropy_han(mesh_viewplane,  number_of_histogram_intervals, current_standard_deviation);
}

float MyRender_OUR::compute_revised_entropy()
{
    // LOG("You'd try this after a viewpoint sampling procedure, or the result will be undefined.\n");

    float e = compute_shannon_entropy_II();
    float sd = current_standard_deviation;
    float avg_e = current_avg_shannon_entropy;
    float avg_sd = current_avg_stand_deviation;

    float re = e - 3*abs(sd-avg_sd)*avg_e/avg_sd;
    // LOG("Revised entropy: "+QString::number(re));
    return re;
}

// 利用多线程来加速计算
/*DWORD WINAPI*/unsigned __stdcall ThreadComputeCurvatureEntropy(LPVOID param)
{
	ThreadViewPlane *vp = (ThreadViewPlane*) param;
	vp->viewPlane->compute_curvatures();

	*(vp->shannon_entropy) = compute_entropy_han(vp->viewPlane, vp->number_of_histogram_intervals , *(vp->standard_deviation));
	//delete vp->viewPlane;
	//nCurrThreadNum --;
	return 0;
}

void MyRender_OUR::viewpoint_selection(bool bMultiThread)
{
    if (!mesh3d)
    {
        int test = 1;// LOG("plz load a mesh first....\n");
        return;
    }
    if (showing_type!= SHOWING_VIEW_DEPENDENT_CURVATURE)
    {
        int test = 1;// LOG("plz select a correct image type first\n");
        return;
    }

    //reset trackball quat, time accumulator
    sgi_trackball_space::trackball(trackball_quat,0,0,0,0);

    //time cost
    clock_t start_t = clock();


	if (bMultiThread)
	{
		ThreadViewPlane *vp[VIEWPOINT_SAMPLES_NUM];
		/*DWORD*/unsigned dwThreadId[VIEWPOINT_SAMPLES_NUM];
		HANDLE hThread[VIEWPOINT_SAMPLES_NUM]; 
		Mesh* vpMesh[VIEWPOINT_SAMPLES_NUM];
		for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
		{
			for (int j = 0; j < i; j++)
			{
				DWORD exitCode = 0;
				GetExitCodeThread(vp[j], &exitCode); 
				if (exitCode !=  STILL_ACTIVE)
				{
					delete vpMesh[j];
					vpMesh[j] = NULL;
					HeapFree(GetProcessHeap(), 0, vp[j]);
				}
			}
			//set auto_sampling_veiwpoint and compute the entropy
			eye_pos_to_bsphere_center = viewpoint_candidates[i]->pos;
			vp[i] = (ThreadViewPlane*)HeapAlloc(GetProcessHeap(), HEAP_ZERO_MEMORY,
				sizeof(ThreadViewPlane));
			vp[i]->viewPlane = vpMesh[i] = GetViewplane();
			//vp[i]->m_mesh_viewplane_idx = new int[fbo_size*fbo_size];
			//memcpy(vp[i]->m_mesh_viewplane_idx, mesh_viewplane_idx, sizeof(int)* fbo_size*fbo_size);
			vp[i]->shannon_entropy = &shannon_entropy[i];
			vp[i]->standard_deviation = &standard_deviations[i];
			vp[i]->number_of_histogram_intervals = number_of_histogram_intervals;
			hThread[i] = (HANDLE)_beginthreadex( NULL, 0, ThreadComputeCurvatureEntropy, vp[i], 0, &dwThreadId[i] );
		}
		// Close all thread handles and free memory allocation.
		for(int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
		{
			WaitForSingleObject(hThread[i],INFINITE);
			revised_entropy[i] = shannon_entropy[i];
			if (vpMesh[i] != NULL)
			{
				delete vpMesh[i];
				HeapFree(GetProcessHeap(), 0, vp[i]);
			}
			CloseHandle(hThread[i]);
		}
	}
	else
	{
		for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
		{
			//set auto_sampling_veiwpoint and compute the entropy
			eye_pos_to_bsphere_center = viewpoint_candidates[i]->pos;
			//paintGL();
			if (mesh_viewplane != NULL)
				delete mesh_viewplane;
			mesh_viewplane = GetViewplane();
			mesh_viewplane->compute_curvatures();

			revised_entropy[i] = shannon_entropy[i] = compute_entropy();
			standard_deviations[i] = current_standard_deviation; //computed in compute_shannon_entropy_II()
		}
	}

    // get mean std_deviation
    current_avg_stand_deviation=0;

    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        current_avg_stand_deviation += standard_deviations[i];
    }
    current_avg_stand_deviation/=VIEWPOINT_SAMPLES_NUM;

    // get mean shannon entropy
    current_avg_shannon_entropy=0;
    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        current_avg_shannon_entropy += shannon_entropy[i];
    }
    current_avg_shannon_entropy /= VIEWPOINT_SAMPLES_NUM;
	/*
	FILE *fp = fopen("around.txt","w");
    // get rebalanced E
    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        standard_deviations[i] -= current_avg_stand_deviation;
        standard_deviations[i] = abs(standard_deviations[i]);
		float tmp = standard_deviations[i]/(current_avg_stand_deviation/current_avg_shannon_entropy);
        //revised_entropy[i] -= 3*standard_deviations[i]/(current_avg_stand_deviation/current_avg_shannon_entropy);
		fprintf(fp,"%f\t%f\t%f\t%f\t%f\n",revised_entropy[i]-1*tmp,revised_entropy[i]-2*tmp,revised_entropy[i]-3*tmp,revised_entropy[i]-4*tmp,revised_entropy[i]-5*tmp);
		revised_entropy[i]-=3*tmp;
    }
	fclose(fp);*/
    // find max revised entropy
    float max_E = revised_entropy[0];
    int max_E_pos = 0;
    for (int i=1; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        if (max_E < revised_entropy[i])
        {
            max_E = revised_entropy[i];
            max_E_pos = i;
        }                                    
    }

	best_viewpoint_idx = max_E_pos;

    clock_t finish_t = clock();
    double sec = double(finish_t - start_t) / CLOCKS_PER_SEC;
    int test = 1;// LOG("Total time elapsed: "+QString::number(sec)+" seconds\n");
    eye_pos_to_bsphere_center = viewpoint_candidates[max_E_pos]->pos;

    // set has been sampled, and the data can be used to compute a certain view's revised entropy
    // sampled_using_revised_entropy_II = true;
}

Mesh* MyRender_OUR::GetViewplane()
{
	if (!model_loaded) return NULL;
	//ust change the packing to ensure no overruns!
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	
	glGetFloatv(GL_MODELVIEW_MATRIX, model_view_project_matrix);
	// 测试使用gl自己的方法来完成矩阵设置 [3/8/2012 Han]
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(/*vertical field of view*/ 45.,
		/*aspect ratio*/ /*(double) viewport.width/viewport.height,*/1,
		/*znear*/ 0.1f, /*zfar*/ 4000.0f);
	glGetFloatv(GL_MODELVIEW_MATRIX, model_view_project_matrix);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
		glGetFloatv(GL_MODELVIEW_MATRIX, model_view_project_matrix);
	glTranslatef(-mesh3d->bsphere.center[0],
		-mesh3d->bsphere.center[1],
		-mesh3d->bsphere.center[2]);
	glGetFloatv(GL_MODELVIEW_MATRIX, model_view_project_matrix);

	gluLookAt( eye_pos_to_bsphere_center[0],
		eye_pos_to_bsphere_center[1],
		eye_pos_to_bsphere_center[2],
		0,0,0,
		0,1,0);
	glGetFloatv(GL_MODELVIEW_MATRIX, model_view_matrix);
	glGetFloatv(GL_PROJECTION_MATRIX, project_matrix);
	multMatrix(model_view_project_matrix, model_view_matrix, project_matrix);
	transposeMatrix(model_view_project_matrix);
	//float rot1[16]; 
	//sgi_trackball_space::build_row_major_rotmatrix(rot1,trackball_quat);
	//glMultMatrixf(rot1);

		//// build the matrices
		//float trans1[16],rot1[16];  // translate and rotate matrix
		//buildTranslateMatrix(
		//	-mesh3d->bsphere.center[0],
		//	-mesh3d->bsphere.center[1],
		//	-mesh3d->bsphere.center[2],
		//	trans1);
		//sgi_trackball_space::build_row_major_rotmatrix(rot1,trackball_quat);
		////set viewpoint according to auto sampling parameter
		//buildLookAtMatrix(
		//	eye_pos_to_bsphere_center[0],
		//	eye_pos_to_bsphere_center[1],
		//	eye_pos_to_bsphere_center[2],
		//	0,0,0,
		//	0,1,0,
		//	view_matrix);
		//buildOrthoMatrix(
		//	-mesh3d->bsphere.r, mesh3d->bsphere.r,
		//	-mesh3d->bsphere.r, mesh3d->bsphere.r,
		//	2*mesh3d->bsphere.r, 4*mesh3d->bsphere.r,
		//	project_matrix);
		////buildPerspectiveMatrix(40,1.0,
		////    2*mesh3d->bsphere.r,
		////    4*mesh3d->bsphere.r,project_matrix);

		//multMatrix(model_matrix, rot1,trans1);
		//multMatrix(model_view_matrix, view_matrix, model_matrix);
		//multMatrix(model_view_project_matrix,project_matrix,model_view_matrix);

		// bind to the fbo
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,fbo_view_dependent_curvature);
		glClearColor(0,0,0,0);
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

		cgGLBindProgram(cg_vprogram_view_dependent_curvature);
		cgGLEnableProfile(cg_vprofile);
		cgGLBindProgram(cg_fprogram_view_dependent_curvature);
		cgGLEnableProfile(cg_fprofile);
		cgSetMatrixParameterfr(
			cgparam_model_view_project_matrix_view_dependent_curvature, 
			model_view_project_matrix);

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_NORMAL_ARRAY);

		glDrawElements(GL_TRIANGLES,(GLsizei)mesh3d->faces.size()*3,GL_UNSIGNED_INT,&mesh3d->faces[0]);
		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);

		//// disable the profiles
		cgGLDisableProfile(cg_vprofile);
		cgGLDisableProfile(cg_fprofile);

		//        // read back radial curvatures, and normalization, and upload to GPU for rendering.
		glReadPixels(0,0,fbo_size,fbo_size,GL_RGB,GL_FLOAT,view_dependent_normal);

#if USE_PROJECTING_COVERED_VIEWPLANE 
		// construct mesh_viewplane, just in the projected area, and a round empty pixel
		// considering the projected area and an empty around pixel adjacent to this area.
		// set all to -1, then set vertex index from 0 to ...
#define IDX(paramx,paramy) ((paramy)*fbo_size+paramx)
		Mesh *viewplane = new Mesh();
		memset(mesh_viewplane_idx,-1,sizeof(int)*fbo_size*fbo_size);
		unsigned int counter=0;
		for (int i=1; i<fbo_size-1; i++) // not consider the edge
		{
			for (int j=1; j<fbo_size-1; j++)
			{
				// if there is any pixel near current is covered, then count it
				bool should_be_a_vertex = false;
				for (int ii=-1; ii<2; ii++)
				{
					for (int jj=-1; jj<2; jj++)
					{
						if( view_dependent_normal[3*IDX(i+ii,j+jj) + 0] != 0 ||
							view_dependent_normal[3*IDX(i+ii,j+jj) + 1] != 0 ||
							view_dependent_normal[3*IDX(i+ii,j+jj) + 2] != 0)
						{
							should_be_a_vertex = true;
							goto IS_A_VERTEX;
						}
					}
				}
				// covered, add to vertex map.
IS_A_VERTEX:
				if (should_be_a_vertex)
				{
					mesh_viewplane_idx[IDX(i,j)] = counter++;
					viewplane->vertices.push_back(point((i+0.5),(j+0.5),0.0));
					viewplane->normals.push_back(point(
						view_dependent_normal[3*IDX(i,j) + 0],
						view_dependent_normal[3*IDX(i,j) + 1],
						view_dependent_normal[3*IDX(i,j) + 2]));
				}
			}
		}
		//		QString o;

		int test = 1;// LOG("vertice:" + QString::number( viewplane->vertices.size()));
		// add to mesh3d data
		for (int i=1; i<fbo_size-1; i++) // not consider the edge
		{
			for (int j=1; j<fbo_size-1; j++)
			{
				if (mesh_viewplane_idx[IDX(i,j)] != -1 &&
					mesh_viewplane_idx[IDX(i+1,j)] != -1 &&
					mesh_viewplane_idx[IDX(i,j+1)] != -1 &&
					mesh_viewplane_idx[IDX(i+1,j+1)] != -1)
				{
					viewplane->faces.push_back(Mesh::Face(
						mesh_viewplane_idx[IDX(i+1,j)],
						mesh_viewplane_idx[IDX(i+1,j+1)],
						mesh_viewplane_idx[IDX(i,j)]));
					viewplane->faces.push_back(Mesh::Face(
						mesh_viewplane_idx[IDX(i+1,j+1)],
						mesh_viewplane_idx[IDX(i,j+1)],
						mesh_viewplane_idx[IDX(i,j)]));
				}
			}
		}
#undef IDX

#else
		// init view dependent curvature view plane mesh,
		static bool init_mesh_viewplane = false;
		if (!init_mesh_viewplane)
		{
			if(mesh_viewplane != 0) delete mesh_viewplane;
			mesh_viewplane = new Mesh();
			// set viewplane vertices
			mesh_viewplane->vertices.resize(fbo_size*fbo_size);
			float mesh_viewplane_step = 1;
			for (int i=0; i<fbo_size; i++)
			{
				for (int j=0; j<fbo_size; j++)
				{
					mesh_viewplane->vertices[j*fbo_size+i][0] = (i+0.5f)*mesh_viewplane_step;
					mesh_viewplane->vertices[j*fbo_size+i][1] = (j+0.5f)*mesh_viewplane_step;
					mesh_viewplane->vertices[j*fbo_size+i][2] = 0;
				}
			}
			// set viewplane faces
			for (int j=0; j<fbo_size-1; j++)
			{
				for (int i=0; i<fbo_size-1; i++)
				{
					mesh_viewplane->faces.push_back(
						Mesh::Face(j*fbo_size+i+1,(j+1)*fbo_size+i+1,j*fbo_size+i));
					mesh_viewplane->faces.push_back(
						Mesh::Face((j+1)*fbo_size+i+1,(j+1)*fbo_size+i,j*fbo_size+i));
				}
			}
			// resize normal vector, and for used by render_scene procedure
			mesh_viewplane->normals.resize(fbo_size*fbo_size);
			init_mesh_viewplane = true;
		}

		for (int i=0; i<fbo_size*fbo_size; i++)
		{
			mesh_viewplane->normals[i][0] = view_dependent_normal[i*3 + 0];
			mesh_viewplane->normals[i][1] = view_dependent_normal[i*3 + 1];
			mesh_viewplane->normals[i][2] = view_dependent_normal[i*3 + 2];
		}
		mesh_viewplane->compute_bsphere();
		mesh_viewplane->curv1.clear();
		mesh_viewplane->curv2.clear();
		mesh_viewplane->compute_curvatures();

#if 0 // for output curvature kinds
		std::ofstream of("a.txt");
		std::map<int,int> cnt;
		for(unsigned int i=0; i<mesh_viewplane->vertices.size(); i++)
		{
			int c = (mesh_viewplane->curv1[i] + mesh_viewplane->curv2[i])/2.0;
			cnt[c]+=1;
		}
		for(std::map<int,int>::iterator it=cnt.begin(); it!=cnt.end(); it++)
		{
			of<<it->first<<' '<<it->second<<"\n";
		}
		of.close();
#endif

		//reset the canvas background
		memset(RGB_mapped_from_curvature,96,
			3*fbo_size*fbo_size*sizeof(unsigned char));
		for (int i=0; i<fbo_size*fbo_size; i++)
		{   // map the curvatures to color space for displaying
			// just for view, so we only handle curvature < 255
			//int idx = mesh_viewplane_idx[i];
			//if (idx == -1) continue; // not covered

			float mc = (mesh_viewplane->curv1[i] + mesh_viewplane->curv2[i])/2.0f;
			mc*=1000;

			if (mc > 255)
			{
				RGB_mapped_from_curvature[i*3+0] = 255;
				RGB_mapped_from_curvature[i*3+1] = 0;
				RGB_mapped_from_curvature[i*3+2] = 0;
			}
			else if (mc > 0 && mc <= 255)
			{
				RGB_mapped_from_curvature[i*3+0] = 255;
				RGB_mapped_from_curvature[i*3+1] = 255 - int(mc);
				RGB_mapped_from_curvature[i*3+2] = 255 - int(mc);
			}
			else if (mc < 0 && mc >= -255)
			{
				RGB_mapped_from_curvature[i*3+0] = 255 + int(mc);
				RGB_mapped_from_curvature[i*3+1] = 255;
				RGB_mapped_from_curvature[i*3+2] = 255 + int(mc);
			}
			else if (mc < -255)
			{
				RGB_mapped_from_curvature[i*3+0] = 0;
				RGB_mapped_from_curvature[i*3+1] = 255;
				RGB_mapped_from_curvature[i*3+2] = 0;
			}
		}
#endif
	// restore default alignment
	glPixelStorei(GL_PACK_ALIGNMENT, 4);
		return viewplane;
}

float MyRender_OUR::viewpoint_selection2(bool bMultiThread)
{
	if (!mesh3d)
	{
		int test = 1;// LOG("plz load a mesh first....\n");
		return 0 ;
	}
	if (showing_type!= SHOWING_VIEW_DEPENDENT_CURVATURE)
	{
		int test = 1;// LOG("plz select a correct image type first\n");
		return 0;
	}

	//reset trackball quat, time accumulator
	sgi_trackball_space::trackball(trackball_quat,0,0,0,0);

	//time cost
	//clock_t start_t = clock();

	//int nCPU = GetCpuNum();
	// 测试，设定视口
	glViewport(0,0,fbo_size,fbo_size);
	glClearColor(0,0,0,0);

	int nMaxThread =/* 2*nCPU*/8;
	int nCurrThreadNum = 0;
	DWORD idx = 0;
	if (bMultiThread)
	{
		ThreadViewPlane **vp = new ThreadViewPlane *[nMaxThread];
		/*DWORD*/unsigned *dwThreadId = new unsigned [nMaxThread];
		HANDLE *hThread = new HANDLE [nMaxThread]; 
		Mesh** vpMesh = new Mesh* [nMaxThread];
		for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
		{
			if (nCurrThreadNum < nMaxThread)
				idx = nCurrThreadNum;
			else
			{
				idx = WaitForMultipleObjects(nMaxThread, hThread,false, INFINITE );
				idx -= WAIT_OBJECT_0;
				delete vpMesh[idx];
				vpMesh[idx] = NULL;
				delete vp[idx];
				vp[idx] = NULL;
				CloseHandle(hThread[idx]);
			}
			//set auto_sampling_veiwpoint and compute the entropy
			eye_pos_to_bsphere_center = viewpoint_candidates[i]->pos;
			vp[idx] = new ThreadViewPlane;/*(ThreadViewPlane*)HeapAlloc(GetProcessHeap(), HEAP_ZERO_MEMORY,
				sizeof(ThreadViewPlane));*/
			vp[idx]->viewPlane = vpMesh[idx] = GetViewplane();
			vp[idx]->shannon_entropy = &shannon_entropy[i];
			vp[idx]->standard_deviation = &standard_deviations[i];
			vp[idx]->number_of_histogram_intervals = number_of_histogram_intervals;
			hThread[idx] = (HANDLE)_beginthreadex( NULL, 0, ThreadComputeCurvatureEntropy, vp[idx], 0, &dwThreadId[idx] );
			nCurrThreadNum++;
		}

		// Close all thread handles and free memory allocation.
		for (int i = 0; i < nMaxThread; i++)
		{
			if (WAIT_FAILED == WaitForSingleObject(hThread[i],INFINITE))
				continue;
			if (vpMesh[i] != NULL)
				delete vpMesh[i];
			if (vp[i] != NULL)
				delete vp[i]/*HeapFree(GetProcessHeap(), 0, vp[i])*/;
			CloseHandle(hThread[i]);
		}
		delete []vp ;
		delete []dwThreadId;
		delete []hThread; 
		delete []vpMesh;


		for(int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
		{
			revised_entropy[i] = shannon_entropy[i];
		}
	}
	else
	{
		for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
		{
			//set auto_sampling_veiwpoint and compute the entropy
			eye_pos_to_bsphere_center = viewpoint_candidates[i]->pos;
			//paintGL();
			if (mesh_viewplane != NULL)
				delete mesh_viewplane;
			mesh_viewplane = GetViewplane();
			mesh_viewplane->compute_curvatures();

			revised_entropy[i] = shannon_entropy[i] = compute_entropy();
			standard_deviations[i] = current_standard_deviation; //computed in compute_shannon_entropy_II()
		}
	}


	// get mean std_deviation
	current_avg_stand_deviation=0;

	for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
	{
		current_avg_stand_deviation += standard_deviations[i];
	}
	current_avg_stand_deviation/=VIEWPOINT_SAMPLES_NUM;

	// get mean shannon entropy
	//current_avg_shannon_entropy=0;
	//for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
	//{
	//	current_avg_shannon_entropy += shannon_entropy[i];
	//}
	//current_avg_shannon_entropy /= VIEWPOINT_SAMPLES_NUM;
	//FILE *fp = fopen("around.txt","w");
	// get rebalanced E

	int alpha = 3;
	

	
	for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
	{
		standard_deviations[i] -= current_avg_stand_deviation;
		standard_deviations[i] = fabs(standard_deviations[i]);
		stdDeviations[i].idx = i;
		stdDeviations[i].standard_deviations = standard_deviations[i];
		//int test = 1;// LOG(QString::number(fabs(standard_deviations[i])));
		
		revised_entropy[i] = -1.0;
		//standard_deviations[i] = abs(standard_deviations[i]);
		//float tmp = standard_deviations[i]/(current_avg_stand_deviation/current_avg_shannon_entropy);
		//revised_entropy[i] -= 3*standard_deviations[i]/(current_avg_stand_deviation/current_avg_shannon_entropy);
		//fprintf(fp,"%f\t%f\t%f\t%f\t%f\n",revised_entropy[i]-1*tmp,revised_entropy[i]-2*tmp,revised_entropy[i]-3*tmp,revised_entropy[i]-4*tmp,revised_entropy[i]-5*tmp);
		//revised_entropy[i]-=3*tmp;
	}
	std::sort(stdDeviations,stdDeviations+VIEWPOINT_SAMPLES_NUM,stdDeviationsLess);

	for (int i=0;i<VIEWPOINT_SAMPLES_NUM;i++)
	{
		revised_entropy[i] = shannon_entropy[i];
	}

	for (int i=0;i<alpha;i++)
	{
		int test = 1;// LOG(QString::number(shannon_entropy[i])+" "+QString::number(stdDeviations[i].idx)+" "+QString::number(stdDeviations[i].standard_deviations));
	}

	//fclose(fp);
	// find max revised entropy
	float max_E = revised_entropy[stdDeviations[0].idx];
	int max_E_pos = stdDeviations[0].idx;
	for (int i=1; i<alpha; i++)
	{
		if (max_E < revised_entropy[stdDeviations[i].idx])
		{
			max_E = revised_entropy[stdDeviations[i].idx];
			max_E_pos = stdDeviations[i].idx;
		}                                    
	}

	best_viewpoint_idx = max_E_pos;
	//clock_t finish_t = clock();
	//double sec = double(finish_t - start_t) / CLOCKS_PER_SEC;
	int test = 1;// LOG("Total time elapsed: "+QString::number(sec)+" seconds\n");
	eye_pos_to_bsphere_center = viewpoint_candidates[best_viewpoint_idx]->pos;
	return revised_entropy[best_viewpoint_idx];

	// set has been sampled, and the data can be used to compute a certain view's revised entropy
	// sampled_using_revised_entropy_II = true;
}


void MyRender_OUR::mouseDrag(int startx, int starty, int endx, int endy, int state)
{
    if (startx > endx)
        std::swap(startx,endx);
    if (starty > endy)
        std::swap(starty,endy);
    if (startx < 0) startx = 0;
    if (starty < 0) starty = 0;
    if (endx > fbo_size-1) endx = fbo_size-1;
    if (endy > fbo_size-1) endy = fbo_size-1;

    if (state == 0) // show the box
    {
        memcpy(parts_color_for_rendering, parts_color_buffer,
            sizeof(unsigned char)*fbo_size*fbo_size*3);
        for (int i=startx; i<endx; i++)
        {
            parts_color_for_rendering[3*((fbo_size-starty)*fbo_size + i) + 0] = 0;
            parts_color_for_rendering[3*((fbo_size-starty)*fbo_size + i) + 1] = 0;
            parts_color_for_rendering[3*((fbo_size-starty)*fbo_size + i) + 2] = 0;

            parts_color_for_rendering[3*((fbo_size-endy)*fbo_size + i) + 0] = 0;
            parts_color_for_rendering[3*((fbo_size-endy)*fbo_size + i) + 1] = 0;
            parts_color_for_rendering[3*((fbo_size-endy)*fbo_size + i) + 2] = 0;
        }
        for (int i=starty; i<endy; i++)
        {
            parts_color_for_rendering[3*((fbo_size - i)*fbo_size + startx) + 0] = 0;
            parts_color_for_rendering[3*((fbo_size - i)*fbo_size + startx) + 1] = 0;
            parts_color_for_rendering[3*((fbo_size - i)*fbo_size + startx) + 2] = 0;

            parts_color_for_rendering[3*((fbo_size - i)*fbo_size + endx) + 0] = 0;
            parts_color_for_rendering[3*((fbo_size - i)*fbo_size + endx) + 1] = 0;
            parts_color_for_rendering[3*((fbo_size - i)*fbo_size + endx) + 2] = 0;
        }
    }

    if (state == 1) // select parts
    {
        // LOG("Selecting: "+QString::number(startx)+","+QString::number(starty)+","+QString::number(endx)+","+QString::number(endy));
        int wmd[sizeof(seg_colors)/sizeof(seg_colors[0])]; //whole model color distribution
        int spd[sizeof(seg_colors)/sizeof(seg_colors[0])]; //selected parts color distribution
        memset(wmd,0,sizeof(int)*sizeof(seg_colors)/sizeof(seg_colors[0]));
        memset(spd,0,sizeof(int)*sizeof(seg_colors)/sizeof(seg_colors[0]));
        unsigned char color_p[3] = {255,255,255}; // previous color
        unsigned char color_c[3] = {0,0,0}; // current color
        bool previous_in_model = false;
        int previous_in_which_part = -1;
        for (int i=0; i<fbo_size*fbo_size; ++i)
        {
            color_c[0] = parts_color_buffer[3*i + 0];
            color_c[1] = parts_color_buffer[3*i + 1];
            color_c[2] = parts_color_buffer[3*i + 2];

            // previous pixel is not in model area, and the next one is the same
            if (!previous_in_model && color_p[0] == color_c[0] &&
                color_p[1] == color_c[1] && color_p[2] == color_c[2])
            {
                continue;
            }
            // previous pixel is not in model area, but the next one is not the same
            else if (!previous_in_model && (color_p[0] != color_c[0] ||
                color_p[1] != color_c[1] || color_p[2] != color_c[2]))
            {
                // find covered by which part
                for (previous_in_which_part = 0; 
                    previous_in_which_part < sizeof(seg_colors)/sizeof(seg_colors[0]); 
                    previous_in_which_part++)
                {
                    if( seg_colors[previous_in_which_part][0] == color_c[0] &&
                        seg_colors[previous_in_which_part][1] == color_c[1] &&
                        seg_colors[previous_in_which_part][2] == color_c[2])
                        break;
                }
                if (previous_in_which_part < sizeof(seg_colors)/sizeof(seg_colors[0]))
                {
                    ++wmd[previous_in_which_part];
                    if (i%fbo_size >= startx && i%fbo_size <= endx &&
                        i/fbo_size >= fbo_size - endy && i/fbo_size <= fbo_size - starty)
                    {
                        ++spd[previous_in_which_part];
                    }
                    previous_in_model = true;
                }
                else
                {
                    previous_in_model = false;
                }
                color_p[0] = color_c[0];
                color_p[1] = color_c[1];
                color_p[2] = color_c[2];
            }
            // previous pixel is in model area, and the next one is the same
            else if (previous_in_model && color_p[0] == color_c[0] &&
                color_p[1] == color_c[1] && color_p[2] == color_c[2])
            {
                ++wmd[previous_in_which_part];
                if (i%fbo_size >= startx && i%fbo_size <= endx &&
                    i/fbo_size >= fbo_size - endy && i/fbo_size <= fbo_size - starty)
                {
                    ++spd[previous_in_which_part];
                }
            }
            // previous pixel is in model area, but the next one is not the same
            else if (previous_in_model && (color_p[0] != color_c[0] ||
                color_p[1] != color_c[1] || color_p[2] != color_c[2]))
            {
                // find covered by which part
                for (previous_in_which_part = 0; 
                    previous_in_which_part < sizeof(seg_colors)/sizeof(seg_colors[0]); 
                    previous_in_which_part++)
                {
                    if( seg_colors[previous_in_which_part][0] == color_c[0] &&
                        seg_colors[previous_in_which_part][1] == color_c[1] &&
                        seg_colors[previous_in_which_part][2] == color_c[2])
                        break;
                }
                if (previous_in_which_part < sizeof(seg_colors)/sizeof(seg_colors[0]))
                {
                    ++wmd[previous_in_which_part];
                    if (i%fbo_size >= startx && i%fbo_size <= endx &&
                        i/fbo_size >= fbo_size - endy && i/fbo_size <= fbo_size - starty)
                    {
                        ++spd[previous_in_which_part];
                    }
                    previous_in_model = true;
                }
                else
                {
                    previous_in_model = false;
                }
                color_p[0] = color_c[0];
                color_p[1] = color_c[1];
                color_p[2] = color_c[2];
            }
        }

        float const ratio = 0.9f;
        std::vector<int> parts;
        // LOG("part");
        for (int i=0; i<sizeof(seg_colors)/sizeof(seg_colors[0]); ++i)
        {
            if (float(spd[i])/wmd[i] > ratio)
            {
                parts.push_back(i);
                int test = 1;// LOG(" "+QString::number(i)+" ");
            }
        }
        // LOG("selected.");
        build_segmentation(parts);

        // update the saved color buffer
        change_mouse_mode(0);
        paintGL();
        change_mouse_mode(1);
        paintGL();
    }
}

// mode = 0, simulate trackball
// mode = 1, select parts by mouse
void MyRender_OUR::change_mouse_mode(int mode)
{
    if (mode == 0)
    {
        selecting_parts_by_mouse = false;
    }
    if (mode == 1)
    {
		//ust change the packing to ensure no overruns!
		glPixelStorei(GL_PACK_ALIGNMENT, 1);

        selecting_parts_by_mouse = true;
        glReadPixels(0,0,fbo_size,fbo_size,GL_RGB,GL_UNSIGNED_BYTE,parts_color_buffer);
        memcpy(parts_color_for_rendering, parts_color_buffer,
            sizeof(unsigned char)*fbo_size*fbo_size*3);
 		// restore default alignment
		glPixelStorei(GL_PACK_ALIGNMENT, 4);
   }
}

void MyRender_OUR::sample_using_revised_entropy_II()
{
    // used for compare
    if (!mesh3d)
    {
        int test = 1;// LOG("plz load a mesh first....\n");
        return;
    }
    if (showing_type != SHOWING_VIEW_DEPENDENT_CURVATURE)
    {
        int test = 1;// LOG("plz select a correct image type first\n");
        return;
    }

    //reset trackball quat, time accumulator
    sgi_trackball_space::trackball(trackball_quat,0,0,0,0);

    //time cost
    clock_t start_t = clock();

    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        //set auto_sampling_veiwpoint and compute the entropy
        eye_pos_to_bsphere_center = viewpoint_candidates[i]->pos;
        paintGL();
        shannon_entropy[i] = compute_shannon_entropy_II();
        standard_deviations[i] = current_standard_deviation; //computed in compute_shannon_entropy_II()
    }

    // get mean std_deviation
    current_avg_stand_deviation=0;

    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        current_avg_stand_deviation += standard_deviations[i];
    }
    current_avg_stand_deviation/=VIEWPOINT_SAMPLES_NUM;

    // get mean shannon entropy
    current_avg_shannon_entropy=0;
    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        current_avg_shannon_entropy += shannon_entropy[i];
    }
    current_avg_shannon_entropy /= VIEWPOINT_SAMPLES_NUM;

    // get rebalanced E
    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        standard_deviations[i] -= current_avg_stand_deviation;
        standard_deviations[i] = abs(standard_deviations[i]);
        revised_entropy[i] -= 3*standard_deviations[i]/(current_avg_stand_deviation/current_avg_shannon_entropy);
    }

    // find max revised entropy
    float max_E = revised_entropy[0];
    int max_E_pos = 0;
    for (int i=1; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        if (max_E < revised_entropy[i])
        {
            max_E = revised_entropy[i];
            max_E_pos = i;
        }                                    
    }

    clock_t finish_t = clock();
    secs = double(finish_t - start_t) / CLOCKS_PER_SEC;
    int test = 1;// LOG("Total time elapsed: "+QString::number(secs)+" seconds\n");
    eye_pos_to_bsphere_center = viewpoint_candidates[max_E_pos]->pos;
}

float MyRender_OUR::compute_shannon_entropy_II()
{
    if (!mesh3d)
    {
        int test = 1;// LOG("plz load a mesh first....\n");
        return -1;
    }
    if (showing_type != SHOWING_VIEW_DEPENDENT_CURVATURE)
    {
        int test = 1;// LOG("plz select a correct image type first\n");
        return -1;
    }

    // curvatures for computing entropy
    std::vector<float> curvatures;
    // use 'as' for abbreviation
    //int as = adaptive_box_size;
    int as = 1;
    // statistic of curvature
    for (int j=0; j<fbo_size; j+=as)
    {
        for (int i=0; i<fbo_size; i+=as)
        {
            // at box (i,j)->(i+as,j+as)
            float mc=0;
            for (int m=0; m<as; m++)
            {
                for (int n=0; n<as; n++)
                {
                    // add mean curvature at (i+n, j+m)
                    mc  +=  mesh_viewplane->curv1[(j+m)*fbo_size+i+n]
                    +mesh_viewplane->curv2[(j+m)*fbo_size+i+n];
                }
            }
            mc/=2; //get real mean curvature
            mc/=as*as; // get average mean curvature over an adaptive box
            curvatures.push_back(mc);
        }
    }

    float maxc=0,minc=0;
    for (unsigned int i=0; i<curvatures.size(); i++)
    {
        if (curvatures[i] > maxc)
        {
            maxc = curvatures[i];
        }
        if (curvatures[i] < minc)
        {
            minc = curvatures[i];
        }
    }
    //int test = 1;// LOG("max mc: "<<maxc<<std::endl;
    //int test = 1;// LOG("min mc: "<<minc<<std::endl;

    const unsigned int HISTO_WIDTH = number_of_histogram_intervals;
    std::vector<unsigned int> histogram;
    float hs = (maxc - minc)/HISTO_WIDTH; //histogram step width
    histogram.resize(HISTO_WIDTH,0);
    for (unsigned int i=0; i<curvatures.size(); i++)
    {
        unsigned int idx = int((curvatures[i] - minc)/hs);
        if (idx < histogram.size())
        {
            histogram[idx] += 1;
        }
    }

    // computing standard deviation
    float avg_mc=0; // average mean curvature,
    for (vector<int>::size_type i=0; i<histogram.size(); i++)
    {
        avg_mc += histogram[i];
    }
    avg_mc/=histogram.size();
    float xigma_mc = 0;
    for (vector<int>::size_type i=0; i<histogram.size(); i++)
    {
        xigma_mc += (avg_mc - histogram[i])*(avg_mc - histogram[i]);
    }
    xigma_mc /= histogram.size() - 1;
    current_standard_deviation = sqrt(xigma_mc);

    // computing entropy 
    float N=0;
    for (vector<int>::size_type i=0; i<histogram.size(); i++)
    {
        N += histogram[i];
    }
    float pi,E=0;
    for (vector<int>::size_type i=0; i<histogram.size(); i++)
    {
        if (histogram[i] != 0)
        {
            pi = histogram[i]/N;
            E += -pi*log(pi)/log(2.0f);
        }
    }

    // LOG("Shannon Entropy II: "+QString::number(E));
    // LOG("Standard Deviation II: "+QString::number(current_standard_deviation));
    return E;
}

void MyRender_OUR::save_viewpoint()
{
    // save result to file:
    std::ofstream msresult("ourresult",std::ios::out | std::ios::app);
    msresult<<secs<<eye_pos_to_bsphere_center<<std::endl;
    msresult.close();
}


void MyRender_OUR::testing_a_round()
{
#define VIEWPOINT_SAMPLES_NUM1 360
	if (!mesh3d)
	{
		std::cout<<"plz load a mesh first....\n";
		return;
	}

	//generating viewpoint candidates
	vector<point> viewpoint_candidates;
	float r = mesh3d->bsphere.r;   //radius of bounding sphere
	// sampling a round 
	for (int i=-90; i<VIEWPOINT_SAMPLES_NUM1-90; i++)
	{
		float alpha = i/360.0f*2*PI;        //angle with positive-x-axis
		point p;
		p[0] = r*cos(alpha);
		p[1] = 0;
		p[2] = r*sin(alpha);
		viewpoint_candidates.push_back(p);
	}

	float entropy[VIEWPOINT_SAMPLES_NUM1];
	float stddeviation[VIEWPOINT_SAMPLES_NUM1];
	float entropy_r[VIEWPOINT_SAMPLES_NUM1];

	std::ofstream of;
	of.open("around.txt",std::ios::out);
	for (int i=0; i<VIEWPOINT_SAMPLES_NUM1; i++)
	{
		//set auto_sampling_veiwpoint and compute the entropy
		eye_pos_to_bsphere_center = viewpoint_candidates[i];

		//showing_image_type = mesh_saliency;
		//paintGL();
		//saliency[i] = compute_mesh_saliency();

		showing_type = SHOWING_VIEW_DEPENDENT_CURVATURE;
		paintGL();
		entropy[i] = compute_entropy();
		stddeviation[i] = current_standard_deviation;
	}

	double esum = 0;
	double ssum = 0;
	for (int i=0; i<VIEWPOINT_SAMPLES_NUM1; i++)
	{
		esum += entropy[i];
		ssum += stddeviation[i];
	}

	for (int i=0; i<VIEWPOINT_SAMPLES_NUM1; i++)
	{
		entropy_r[i] = entropy[i] - 3*abs(stddeviation[i] - ssum/VIEWPOINT_SAMPLES_NUM1)*esum/ssum;
	}


	for (int i=0; i<VIEWPOINT_SAMPLES_NUM1; i++)
	{
		of<<i;
		of<<' '<<entropy[i];
		of<<' '<<stddeviation[i];
		of<<' '<<entropy_r [i]<<std::endl;
	}
	of.close();

	//static int ep = 0;
	//std::cout<<ep<<std::endl;
	//eye_pos_to_bsphere_center = viewpoint_candidates[ep++];
#undef VIEWPOINT_SAMPLES_NUM1
}


void MyRender_OUR::show_A()
{
	float maxE = -10;
	int maxEk = 0;
	for (int i=0;i<VIEWPOINT_SAMPLES_NUM;i++)
	{
		if (maxE<shannon_entropy[i]) 
		{
			maxE=shannon_entropy[i];
			maxEk = i;
		}
	}
	eye_pos_to_bsphere_center = viewpoint_candidates[maxEk]->pos;
	int test = 1;// LOG(QString::number(shannon_entropy[maxEk])+" " + QString::number(standard_deviations[maxEk]));
}

void MyRender_OUR::show_B1()
{
	eye_pos_to_bsphere_center = viewpoint_candidates[stdDeviations[0].idx]->pos;
	int test = 1;// LOG(QString::number(shannon_entropy[stdDeviations[0].idx])+" " + QString::number(stdDeviations[0].standard_deviations));
}

void MyRender_OUR::show_B2()
{
	eye_pos_to_bsphere_center = viewpoint_candidates[stdDeviations[1].idx]->pos;
	int test = 1;// LOG(QString::number(shannon_entropy[stdDeviations[1].idx])+" " + QString::number(stdDeviations[1].standard_deviations));
}

void MyRender_OUR::show_B3()
{
	eye_pos_to_bsphere_center = viewpoint_candidates[stdDeviations[2].idx]->pos;
	int test = 1;// LOG(QString::number(shannon_entropy[stdDeviations[2].idx])+" " + QString::number(stdDeviations[2].standard_deviations));
}


int MyRender_OUR::GetMeshTriNum()
{
	return mesh3d->faces.size();
}

point MyRender_OUR::GetViewPoint()
{
	if (best_viewpoint_idx >= viewpoint_candidates_spherical.size())
		return point(0, 0, 0);
	return viewpoint_candidates_spherical[best_viewpoint_idx];
}

void MyRender_OUR::RandomDistributeViews(int numViews, float r)
{
	viewpoint_candidates.clear();
	viewpoint_candidates_spherical.clear();
	time_t tm;
	srand((unsigned)time(&tm));  //set seed;

	for (int i=0; i<numViews; i++)
	{
		float alpha = (rand()%360)/360.0f*2*PI;        //angle with positive-x-axis
		float beta = (rand()%360)/360.0f*2*PI - PI;  //angle with xy-plane
		ViewPoint *vp = new ViewPoint;
		point p;
		p[0] = r*cos(beta)*sin(alpha);
		p[1] = r*cos(beta)*cos(alpha);
		p[2] = r*sin(beta);
		vp->pos = p;
		viewpoint_candidates.push_back(vp);

		p[0] = alpha;
		p[1] = beta;
		p[2] = r;
		viewpoint_candidates_spherical.push_back(p);		
	}
}

/*
   Return the midpoint between two vectors
*/
point MidPoint(point p1,point p2)
{
	point p;

   p[0] = (p1[0] + p2[0]) / 2;
   p[1] = (p1[1] + p2[1]) / 2;
   p[2] = (p1[2] + p2[2]) / 2;

   return(p);
}

/*
   Normalise a vector
*/
void Normalise(point *p, double r)
{
	double length;

	length = sqrt((*p)[0] * (*p)[0] + (*p)[1] * (*p)[1] + (*p)[2] * (*p)[2]);
	length /= r;
	if (length != 0) {
		(*p)[0] /= length;
		(*p)[1] /= length;
		(*p)[2] /= length;
	} else {
		(*p)[0] = 0;
		(*p)[1] = 0;
		(*p)[2] = 0;
	}
}
/*
http://paulbourke.net/miscellaneous/sphere_cylinder/
   Create a triangular facet approximation to a sphere
   Return the number of facets created.
   The number of facets will be (4^iterations) * 8
*/
// 对球进行平均划分，划分后的顶点作为视点 [5/15/2012 Han]
void MyRender_OUR::EqualDistributeViews(int numViewsH, int numViewsV, float r)
{
	vector<ViewPoint*>::iterator iter;
	for (iter = viewpoint_candidates.begin(); iter != viewpoint_candidates.end(); ++iter)
		delete *iter;
	viewpoint_candidates.clear();
	viewpoint_candidates_spherical.clear();

		//////////////////////////////////////////////////////////////////////////
		// 采用球面平均分布点的方式进行视点初始化
	int iterations = 3;
	/*r = 2.7;*/
	//f = (FACET3 *)calloc(pow(4.0,iterations) * 8, sizeof(FACET3));
	vf.clear();
	vf.resize(pow(4.0,iterations) * 8);
	int i,it;
	double a;
	point p[6] = {point(0,0,r),  point(0,0,-r),  point(-r,-r,0),  point(r,-r,0),  point(r,r,0), point(-r,r,0)};
	point pa,pb,pc;
	int nt = 0,ntold;

	/* Create the level 0 object */
	a = 1 / sqrt(2.0);
	for (i=0;i<6;i++) {
		p[i][0] *= a;
		p[i][1] *= a;

		ViewPoint *pV = new ViewPoint;
		pV->pos[0] = p[i][0];
		pV->pos[1] = p[i][1];
		pV->pos[2] = p[i][2];
		viewpoint_candidates.push_back(pV);

	}
	vf[0].p1 = 0; vf[0].p2 = 3; vf[0].p3 = 4;
	vf[1].p1 = 0; vf[1].p2 = 4; vf[1].p3 = 5;
	vf[2].p1 = 0; vf[2].p2 = 5; vf[2].p3 = 2;
	vf[3].p1 = 0; vf[3].p2 = 2; vf[3].p3 = 3;
	vf[4].p1 = 1; vf[4].p2 = 4; vf[4].p3 = 3;
	vf[5].p1 = 1; vf[5].p2 = 5; vf[5].p3 = 4;
	vf[6].p1 = 1; vf[6].p2 = 2; vf[6].p3 = 5;
	vf[7].p1 = 1; vf[7].p2 = 3; vf[7].p3 = 2;
	nt = 8;

	if (iterations < 1)
		return;

	/* Bisect each edge and move to the surface of a unit sphere */
	for (it=0;it<iterations;it++) {
		ntold = nt;
		for (i=0;i<ntold;i++) {
			pa = (viewpoint_candidates.at(vf[i].p1)->pos + viewpoint_candidates.at(vf[i].p2)->pos)/ 2.0f;
			//pa[0] = (vf[i].p1.pos[0] + vf[i].p2.pos[0]) / 2;
			//pa[1] = (vf[i].p1.pos[1] + vf[i].p2.pos[1]) / 2;
			//pa[2] = (vf[i].p1.pos[2] + vf[i].p2.pos[2]) / 2;
			pb = (viewpoint_candidates.at(vf[i].p2)->pos + viewpoint_candidates.at(vf[i].p3)->pos)/ 2.0f;
			//pb[0] = (vf[i].p2.pos[0] + vf[i].p3.pos[0]) / 2;
			//pb[1] = (vf[i].p2.pos[1] + vf[i].p3.pos[1]) / 2;
			//pb[2] = (vf[i].p2.pos[2] + vf[i].p3.pos[2]) / 2;
			pc = (viewpoint_candidates.at(vf[i].p3)->pos + viewpoint_candidates.at(vf[i].p1)->pos)/ 2.0f;
			//pc[0] = (vf[i].p3.pos[0] + vf[i].p1.pos[0]) / 2;
			//pc[1] = (vf[i].p3.pos[1] + vf[i].p1.pos[1]) / 2;
			//pc[2] = (vf[i].p3.pos[2] + vf[i].p1.pos[2]) / 2;
			Normalise(&pa, r);
			Normalise(&pb, r);
			Normalise(&pc, r);
			// 将新生成的点插入视点当中 [5/15/2012 Han]
			ViewPoint *pV = new ViewPoint;
			unsigned int paI, pbI, pcI;
			paI = viewpoint_candidates.size();
			pbI = paI+1, pcI = paI+2;

			pV->pos = pa;
			viewpoint_candidates.push_back(pV);

			pV = new ViewPoint;
			pV->pos = pb;
			viewpoint_candidates.push_back(pV);

			pV = new ViewPoint;
			pV->pos = pc;
			viewpoint_candidates.push_back(pV);

			vf[nt].p1 = vf[i].p1; vf[nt].p2 = paI/*pa*/; vf[nt].p3 = pcI/*pc*/; nt++;
			vf[nt].p1 = paI/*pa*/; vf[nt].p2 = vf[i].p2; vf[nt].p3 = pbI/*pb*/; nt++;
			vf[nt].p1 = pbI/*pb*/; vf[nt].p2 = vf[i].p3; vf[nt].p3 = pcI/*pc*/; nt++;
			vf[i].p1 = paI/*pa*/;
			vf[i].p2 = pbI/*pb*/;
			vf[i].p3 = pcI/*pc*/;
		}
	}
}

void MyRender_OUR::GetBestViewDist(float dis, Image_Type type, FILE *outPutFile)
{
	Image_Type ba = showing_type;
	showing_type =  type;
	time_t tm;
	srand((unsigned)time(&tm));  //set seed;

	float backUp = mesh3d->bsphere.r;
	//mesh3d->bsphere.r = dis;
	EqualDistributeViews(10, 10, dis/*mesh3d->bsphere.r*/);
	wchar_t fn[128];
	swprintf( fn,   L"%d_%f_%d.ppm ", mesh3d->faces.size(),dis,rand());
	wchar_t fn_original[128];
	swprintf( fn_original, L"%d_%f_%d_o.ppm ", mesh3d->faces.size(),dis,rand());
	double viewSlecTime; 
	float fInfo = 0;
	double sT = 0;
	switch (type)
	{
	case SHOWING_MESH_SALIENCY:
		// 利用重要度进行模型视点选择
		// 使用文章的算法:找到可见顶点saliency之和最大的视点即为最优视点
		if (mesh3d->saliency.size() == 0)
			TIMING(sT, compute_mesh_saliency_3d());
		TIMING(viewSlecTime, fInfo = viewpoint_selection_saliency());
		if (outPutFile != NULL)
			fprintf(outPutFile, "使用saliency进行视点选择,saliency计算时间：%lf秒\n", sT);
		break;
		//case model_space_curvature_image:
		//	break;
	default:
		TIMING(viewSlecTime, fInfo = viewpoint_selection2(true));
		break;
	}
	// 保存该视点的渲染信息
	if (outPutFile != NULL)
	{
		fprintf(outPutFile, 
			"最优视点编号：%d,位置：α：%f，β：%f，距离：%f\n最优视点的熵：%f\n视点选择用时：%lf秒\n"
			,best_viewpoint_idx, GetViewPoint()[0], GetViewPoint()[1], GetViewPoint()[2],fInfo, viewSlecTime);
		//fprintf(outPutFile, "所有候选视点revised_entropy：\n");			
		//for (int j = 0; j < VIEWPOINT_SAMPLES_NUM; j++)
		//{
		//	//for (int k = 0; k < 3; k++)
		//		fprintf(outPutFile, "%d:%f\t", j, revised_entropy[j]);
		//}
		//fprintf(outPutFile, "\n");
	}
	// 保存得到视点对应的渲染图像 [2/24/2012 Han]
	GLint viewport[4];											//space for viewport data
	glGetIntegerv (GL_VIEWPORT, viewport);						//get actual viewport
	int w = viewport[2], h = viewport[3];

	//ust change the packing to ensure no overruns!
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	PPM *ppmWriter = new PPM;
	ppmWriter->width = w;
	ppmWriter->height = h;
	ppmWriter->version = "P6";
	ppmWriter->data = new unsigned char [w*h*3];
	paintGL();
	glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, ppmWriter->data);
	ppmWriter->save(fn);

	// 保存真实着色图像 [3/7/2012 Han]
	showing_type =  SHOWING_MODEL;
	paintGL();
	glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, ppmWriter->data);
	ppmWriter->save(fn_original);
	showing_type = ba;

	delete ppmWriter;
	// restore default alignment
	glPixelStorei(GL_PACK_ALIGNMENT, 4);
	mesh3d->bsphere.r = backUp;

	ResetEyePos();
	///////////////////////////////////////////////////////////////////////////
	// test
	ViewSphereColor();
}

void MyRender_OUR::ViewSphereColor()
{
	double maxI = 0;
	double minI = FLT_MAX;
	for (int i = 0; i < viewpoint_candidates.size(); i++)
	{
		if (maxI < viewpoint_candidates.at(i)->importance)
			maxI = viewpoint_candidates.at(i)->importance;
		if (minI > viewpoint_candidates.at(i)->importance)
			minI = viewpoint_candidates.at(i)->importance;
	}


	double R;
	int indiceLut;
	for (int i = 0; i < viewpoint_candidates.size(); i++)
	{
		R=((viewpoint_candidates.at(i)->importance-minI)/(maxI-minI))*255;

		if(R>255)
			R=255;
		indiceLut=floor(R);
		viewpoint_candidates.at(i)->color[0] = LUT_CourbureClust[3*indiceLut];
		viewpoint_candidates.at(i)->color[1] = LUT_CourbureClust[3*indiceLut+1];
		viewpoint_candidates.at(i)->color[2] = LUT_CourbureClust[3*indiceLut+2];
	}
}

bool MyRender_OUR::IsOccluded(float p[3], GLfloat *bufferZ)
{
	GLint viewport[4];											//space for viewport data
	GLdouble mvmatrix[16], projmatrix[16];  //space for transform matricex
	GLdouble winx, winy, winz;							//space for returned projected coords
	GLdouble flareZ;												//here we will store the transformed flare Z

	// Now we will ask OGL to project some geometry for us using the gluProject function.
	// Practically we ask OGL to guess where a point in space will be projected in our current viewport,
	// using arbitrary viewport and transform matrices we pass to the function.
	// If we pass to the function the current matrices  (retrievede with the glGet funcs)
	// we will have the real position on screen where the dot will be drawn.
	// The interesting part is that we also get a Z value back, this means that 
	// reading the REAL buffer for Z values we can discover if the flare is in front or
	// if it's occluded by some objects.


	glGetIntegerv (GL_VIEWPORT, viewport);						//get actual viewport
	glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);			//get actual model view matrix
	glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);	//get actual projiection matrix

	// this asks OGL to guess the 2d position of a 3d point inside the viewport
	gluProject(p[0], p[1], p[2], mvmatrix, projmatrix, viewport, &winx, &winy, &winz);
	flareZ = winz;

	// test [5/7/2012 Han]
	flareZ -= 0.001f;

	// we read back one pixel from th depth buffer (exactly where our flare should be drawn)
	//glReadPixels(winx, winy,1,1,GL_DEPTH_COMPONENT, GL_FLOAT, &bufferZ);

	// if the buffer Z is lower than our flare guessed Z then don't draw 
	// this means there is something in front of our flare
	if (winx < 0 || winx > viewport[2]-1 || winy < 0 || winy > viewport[3]-1)
		return true;
	int l, r, u, d;
	l = int(winx);
	l = l < 0 ? 0:l;
	r = int(winx+1);
	r = r > viewport[2] - 1 ? viewport[2] - 1 : r;
	u = int(winy+1);
	u = u > viewport[3] - 1 ? viewport[3] - 1 : u;
	d = int(winy);
	d = d < 0 ? 0 : d;
	if (bufferZ[u*viewport[2] + l] < flareZ
		&&bufferZ[u*viewport[2] + r] < flareZ
		&&bufferZ[d*viewport[2] + l] < flareZ
		&&bufferZ[d*viewport[2] + r] < flareZ)
		return true;
	else
		return false;
}

void MyRender_OUR::GetVisibleTris(vector<int> &tris)
{
	//ust change the packing to ensure no overruns!
	glPixelStorei(GL_PACK_ALIGNMENT, 1);

	Image_Type ba = showing_type;
	showing_type =  VISIBLE_TEST;

	//////////////////////////////////////////////////////////////////////////
	// 利用三角形三个顶点以及中点的深度和当前像素深度比较来决定三角形是否可见
	GLint viewport[4];											//space for viewport data
	glGetIntegerv (GL_VIEWPORT, viewport);						//get actual viewport
	int w = viewport[2], h = viewport[3];

	GLfloat *zBuffer = new GLfloat [w*h]; 
	paintGL();		// render model

	//glReadPixels(0,0,w,h,GL_RGB,GL_FLOAT,zBuffer);
	//glReadPixels(0,0,w,h,GL_COLOR_INDEX,GL_UNSIGNED_INT,zBuffer);
	// we read back one pixel from th depth buffer (exactly where our flare should be drawn)
	glReadPixels(0, 0,w,h,GL_DEPTH_COMPONENT, GL_FLOAT, zBuffer);
	
	float center_T[3];
	for (unsigned int i=0; i<mesh3d->faces.size(); i++)
	{
		// 三角形中点
		for(int n = 0 ; n < 3; n++)
		{
			center_T[n] = 0;
			for (int m = 0; m < 3; m++)
			{
				center_T[n] += 	mesh3d->vertices[mesh3d->faces[i][m]][n];
			}
			center_T[n] /= 3.0f;
		}
		if (/*!IsOccluded(center_T, zBuffer)||!IsOccluded(mesh3d->vertices[mesh3d->faces[i][0]], zBuffer)||!IsOccluded(mesh3d->vertices[mesh3d->faces[i][1]], zBuffer)
			||!IsOccluded(mesh3d->vertices[mesh3d->faces[i][2]], zBuffer)*/
			!IsOccluded(center_T, zBuffer)||(!IsOccluded(mesh3d->vertices[mesh3d->faces[i][0]], zBuffer)&&!IsOccluded(mesh3d->vertices[mesh3d->faces[i][1]], zBuffer)
			&&!IsOccluded(mesh3d->vertices[mesh3d->faces[i][2]], zBuffer)))
		{
			tris.push_back(i);
		}
	}
	delete []zBuffer;
	// restore default alignment
	glPixelStorei(GL_PACK_ALIGNMENT, 4);
	//// 三个颜色通道都使用，每种颜色均分，得到对应的三角形编号
	//// 颜色1不使用
	//int num_t = mesh3d->faces.size()/3 + 1;
	//for (int i = 0; i < w*h; i++)
	//{
	//	// test delete
	//	tris.push_back(zBuffer[i]);
	//	continue;


	//	for(int n = 0; n < 3; n++)
	//	{
	//		if (zBuffer[i*3+n] != 1.0f)
	//		{
	//			int tId = int(zBuffer[i*3+n] * (num_t) + 0.5f);	// 避免数值误差
	//			tId += n * num_t;
	//			tris.push_back(tId);
	//			break;
	//		}
	//	}
	//}
	//delete []zBuffer;

	//// test delete [5/1/2012 Han]
	//// 将可见的面片绘制成红色，其他面片时白色，观察是否正确
	//glClearColor(1,1,1,1);
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	////glBegin(GL_TRIANGLES);
	////for (unsigned int i=0; i<mesh3d->faces.size(); i++)
	////{
	////	glColor3f(0, 0, 1);
	////	//glPassThrough(i);
	////	//glLoadName(i);										// Assign Object A Name (ID)
	////	glVertex3f(
	////		mesh3d->vertices[mesh3d->faces[i][0]][0],
	////		mesh3d->vertices[mesh3d->faces[i][0]][1],
	////		mesh3d->vertices[mesh3d->faces[i][0]][2]);
	////	glVertex3f(
	////		mesh3d->vertices[mesh3d->faces[i][1]][0],
	////		mesh3d->vertices[mesh3d->faces[i][1]][1],
	////		mesh3d->vertices[mesh3d->faces[i][1]][2]);
	////	glVertex3f(
	////		mesh3d->vertices[mesh3d->faces[i][2]][0],
	////		mesh3d->vertices[mesh3d->faces[i][2]][1],
	////		mesh3d->vertices[mesh3d->faces[i][2]][2]);

	////}
	////glEnd();

	//glBegin(GL_TRIANGLES);
	//for (unsigned int i=0; i<tris.size(); i++)
	//{
	//	if (tris[i] >= mesh3d->faces.size())
	//		continue;

	//	glColor3f(1, 0, 0);
	//	//glPassThrough(i);
	//	//glLoadName(i);										// Assign Object A Name (ID)
	//	glVertex3f(
	//		mesh3d->vertices[mesh3d->faces[tris[i]][0]][0],
	//		mesh3d->vertices[mesh3d->faces[tris[i]][0]][1],
	//		mesh3d->vertices[mesh3d->faces[tris[i]][0]][2]);
	//	glVertex3f(
	//		mesh3d->vertices[mesh3d->faces[tris[i]][1]][0],
	//		mesh3d->vertices[mesh3d->faces[tris[i]][1]][1],
	//		mesh3d->vertices[mesh3d->faces[tris[i]][1]][2]);
	//	glVertex3f(
	//		mesh3d->vertices[mesh3d->faces[tris[i]][2]][0],
	//		mesh3d->vertices[mesh3d->faces[tris[i]][2]][1],
	//		mesh3d->vertices[mesh3d->faces[tris[i]][2]][2]);
	//}
	//glEnd();
	//
	//wchar_t fn[128];
	//swprintf( fn,   L"test.ppm");
	//PPM *ppmWriter = new PPM;
	//ppmWriter->width = w;
	//ppmWriter->height = h;
	//ppmWriter->version = "P6";
	//ppmWriter->data = new unsigned char [w*h*3];
	////	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	//glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, ppmWriter->data);
	//ppmWriter->save(fn);
	//delete ppmWriter;
	
	////////////////////////////////////////////////////////////////////////////
	//// 使用feedback
	//int step = 2;
	//int bnum = mesh3d->faces.size()*step;
	//GLfloat	*buffer = new GLfloat [bnum];		// Set Up A Selection Buffer
	//GLint size;
	//glFeedbackBuffer (bnum, GL_3D, buffer);
	//(void) glRenderMode (GL_FEEDBACK);

	//paintGL();

	//size = glRenderMode (GL_RENDER);
	//GLfloat token;
	//CString s;
	//for (int loop = 0; loop < size; loop+=step)					// Loop Through All The Detected Hits
	//{
	//	token = buffer[loop];
	//	if (token == GL_PASS_THROUGH_TOKEN) {
	//		float t = buffer[loop + 1];
	//		int test = 1;
	//		sprintf(s.GetBuffer(), "%f\n", t);
	//		OutputDebugString(s);
	//	}
	//}
	//delete []buffer;


	//////////////////////////////////////////////////////////////////////////
	// 使用selection
	// This Sets The Array <viewport> To The Size And Location Of The Screen Relative To The Window
	//glGetIntegerv(GL_VIEWPORT, viewport);
	//glSelectBuffer(bnum, buffer);								// Tell OpenGL To Use Our Array For Selection

	//// Puts OpenGL In Selection Mode. Nothing Will Be Drawn.  Object ID's and Extents Are Stored In The Buffer.
	//(void) glRenderMode(GL_SELECT);

	//glInitNames();												// Initializes The Name Stack
	//glPushName(0);												// Push 0 (At Least One Entry) Onto The Stack

	//paintGL();

	//hits=glRenderMode(GL_RENDER);								// Switch To Render Mode, Find Out How Many
	//// Objects Were Drawn Where The Mouse Was
	//if (hits > 0)												// If There Were More Than 0 Hits
	//{
	//	CMeshSimpDoc* pDoc = GetDocument();
	//	ASSERT_VALID(pDoc);

	//	int	choose = buffer[3];									// Make Our Selection The First Object
	//	int depth = buffer[1];									// Store How Far Away It Is 

	//	for (int loop = 1; loop < hits; loop++)					// Loop Through All The Detected Hits
	//	{
	//		// If This Object Is Closer To Us Than The One We Have Selected
	//		if (buffer[loop*4+1] < GLuint(depth))
	//		{
	//			choose = buffer[loop*4+3];						// Select The Closer Object
	//			depth = buffer[loop*4+1];						// Store How Far Away It Is
	//		}       
	//	}
	//	triangle &t = pDoc->m_pProgMesh->getTri(choose);
	//	t.setSel(true);
	//	m_selIndexList.remove(choose);
	//	m_selIndexList.push_back(choose);
	//	m_bSelAny = true;
	//	Invalidate();
	//}
	showing_type = ba;
}

float MyRender_OUR::viewpoint_selection_saliency(bool bMultiThread)
{
	if (!mesh3d)
	{
		int test = 1;// LOG("plz load a mesh first....\n");
		return 0;
	}	
	//reset trackball quat, time accumulator
	sgi_trackball_space::trackball(trackball_quat,0,0,0,0);
	
	int nMaxThread =/* 2*nCPU*/8;
	int nCurrThreadNum = 0;
	DWORD idx = 0;
	if (bMultiThread)
	{
		// 使用多线程 ，暂不实现 [4/28/2012 Han]
		//ThreadViewPlane **vp = new ThreadViewPlane *[nMaxThread];
		///*DWORD*/unsigned *dwThreadId = new unsigned [nMaxThread];
		//HANDLE *hThread = new HANDLE [nMaxThread]; 
		//Mesh** vpMesh = new Mesh* [nMaxThread];
		//for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
		//{
		//	if (nCurrThreadNum < nMaxThread)
		//		idx = nCurrThreadNum;
		//	else
		//	{
		//		idx = WaitForMultipleObjects(nMaxThread, hThread,false, INFINITE );
		//		idx -= WAIT_OBJECT_0;
		//		delete vpMesh[idx];
		//		vpMesh[idx] = NULL;
		//		delete vp[idx];
		//		vp[idx] = NULL;
		//		CloseHandle(hThread[idx]);
		//	}
		//	//set auto_sampling_veiwpoint and compute the entropy
		//	eye_pos_to_bsphere_center = viewpoint_candidates[i];
		//	vp[idx] = new ThreadViewPlane;/*(ThreadViewPlane*)HeapAlloc(GetProcessHeap(), HEAP_ZERO_MEMORY,
		//		sizeof(ThreadViewPlane));*/
		//	vp[idx]->viewPlane = vpMesh[idx] = GetViewplane();
		//	vp[idx]->shannon_entropy = &shannon_entropy[i];
		//	vp[idx]->standard_deviation = &standard_deviations[i];
		//	vp[idx]->number_of_histogram_intervals = number_of_histogram_intervals;
		//	hThread[idx] = (HANDLE)_beginthreadex( NULL, 0, ThreadComputeCurvatureEntropy, vp[idx], 0, &dwThreadId[idx] );
		//	nCurrThreadNum++;
		//}

		//// Close all thread handles and free memory allocation.
		//for (int i = 0; i < nMaxThread; i++)
		//{
		//	if (WAIT_FAILED == WaitForSingleObject(hThread[i],INFINITE))
		//		continue;
		//	if (vpMesh[i] != NULL)
		//		delete vpMesh[i];
		//	if (vp[i] != NULL)
		//		delete vp[i]/*HeapFree(GetProcessHeap(), 0, vp[i])*/;
		//	CloseHandle(hThread[i]);
		//}
		//delete []vp ;
		//delete []dwThreadId;
		//delete []hThread; 
		//delete []vpMesh;


		//for(int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
		//{
		//	revised_entropy[i] = shannon_entropy[i];
		//}
	}
	else
	{
		for (int i=0; i<viewpoint_candidates.size(); i++)
		{
			ViewPoint *vp = viewpoint_candidates.at(i);
			//set auto_sampling_veiwpoint and compute the entropy
			eye_pos_to_bsphere_center = viewpoint_candidates.at(i)->pos;
			vp->importance = 0;
			vector<int> visible_tris;
			GetVisibleTris(visible_tris);
			for(int t = 0; t < visible_tris.size(); t++)
			{
				for(int v = 0; v < 3; v++)
				{
					vp->importance += mesh3d->saliency[mesh3d->faces[visible_tris[t]][v]];
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// 将可见的面片绘制成红色，其他面片时白色，观察是否正确
	//eye_pos_to_bsphere_center = viewpoint_candidates[0];
	//vector<int> visible_tris;
	//GetVisibleTris(visible_tris);
	//GLint viewport[4];											//space for viewport data
	//glGetIntegerv (GL_VIEWPORT, viewport);						//get actual viewport
	//int w = viewport[2], h = viewport[3];
	//PPM *ppmWriter = new PPM;
	////ust change the packing to ensure no overruns!
	//glPixelStorei(GL_PACK_ALIGNMENT, 1);

	//ppmWriter->width = w;
	//ppmWriter->height = h;
	//ppmWriter->version = "P6";
	//ppmWriter->data = new unsigned char [w*h*3];
	//for (int i=0; i</*VIEWPOINT_SAMPLES_NUM*/1; i++)
	//{
	//	//set auto_sampling_veiwpoint and compute the entropy
	//	eye_pos_to_bsphere_center = viewpoint_candidates[i];
	//	paintGL();		// render model
	//	glClearColor(1,1,1,1);
	//	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//	glBegin(GL_TRIANGLES);
	//	for (unsigned int t=0; t<visible_tris.size(); t++)
	//	{
	//		if (visible_tris[i] >= mesh3d->faces.size())
	//			continue;

	//		glColor3f(1, 0, 0);
	//		//glPassThrough(i);
	//		//glLoadName(i);										// Assign Object A Name (ID)
	//		glVertex3f(
	//			mesh3d->vertices[mesh3d->faces[visible_tris[t]][0]][0],
	//			mesh3d->vertices[mesh3d->faces[visible_tris[t]][0]][1],
	//			mesh3d->vertices[mesh3d->faces[visible_tris[t]][0]][2]);
	//		glVertex3f(
	//			mesh3d->vertices[mesh3d->faces[visible_tris[t]][1]][0],
	//			mesh3d->vertices[mesh3d->faces[visible_tris[t]][1]][1],
	//			mesh3d->vertices[mesh3d->faces[visible_tris[t]][1]][2]);
	//		glVertex3f(
	//			mesh3d->vertices[mesh3d->faces[visible_tris[t]][2]][0],
	//			mesh3d->vertices[mesh3d->faces[visible_tris[t]][2]][1],
	//			mesh3d->vertices[mesh3d->faces[visible_tris[t]][2]][2]);
	//	}
	//	glEnd();

	//	wchar_t fn[128];
	//	swprintf( fn,   L"test%d.ppm", i);
	//	//	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	//	glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, ppmWriter->data);
	//	ppmWriter->save(fn);
	//}
	//delete ppmWriter;
	//// restore default alignment
	//glPixelStorei(GL_PACK_ALIGNMENT, 4);
	//////////////////////////////////////////////////////////////////////////


	float max_s = 0;
	best_viewpoint_idx = -1; 
	for (int i=0; i<viewpoint_candidates.size(); i++)
	{
		if (viewpoint_candidates.at(i)->importance > max_s)
		{
			max_s = viewpoint_candidates.at(i)->importance;
			best_viewpoint_idx = i;
		}
	}
	eye_pos_to_bsphere_center = viewpoint_candidates[best_viewpoint_idx]->pos;
	return viewpoint_candidates.at(best_viewpoint_idx)->importance;
}




//////////////////////////////////////////////////////////////////////////

/*
GL_LINE_RESET_TOKEN
	30.00 30.00 0.00 0.84 0.84 0.84 1.00
	50.00 60.00 0.00 0.84 0.84 0.84 1.00
	GL_LINE_TOKEN
	50.00 60.00 0.00 0.84 0.84 0.84 1.00
	70.00 40.00 0.00 0.84 0.84 0.84 1.00
	GL_PASS_THROUGH_TOKEN
	1.00
	GL_PASS_THROUGH_TOKEN
	2.00
	GL_POINT_TOKEN
	50.00 50.00 0.00 0.84 0.84 0.84 1.00
void drawGeometry (GLenum mode)
{
	glBegin (GL_LINE_STRIP);
	glNormal3f (0.0, 0.0, 1.0);
	glVertex3f (30.0, 30.0, 0.0);
	glVertex3f (50.0, 60.0, 0.0);
	glVertex3f (70.0, 40.0, 0.0);
	glEnd ();
	if (mode == GL_FEEDBACK)
		glPassThrough (1.0);

	glBegin (GL_POINTS);
	glVertex3f (-100.0, -100.0, -100.0); 
	glEnd ();
	if (mode == GL_FEEDBACK)
		glPassThrough (2.0);

	glBegin (GL_POINTS);
	glNormal3f (0.0, 0.0, 1.0);
	glVertex3f (50.0, 50.0, 0.0);
	glEnd ();
}

void print3DcolorVertex (GLint size, GLint *count, 
	GLfloat *buffer)
{
	int i;

	printf ("  ");
	for (i = 0; i < 7; i++) {
		printf ("%4.2f ", buffer[size-(*count)]);
		*count = *count - 1;
	}
	printf ("\n");
}

void printBuffer(GLint size, GLfloat *buffer)
{
	GLint count;
	GLfloat token;

	count = size;
	while (count) {
		token = buffer[size-count]; count--;
		if (token == GL_PASS_THROUGH_TOKEN) {
			printf ("GL_PASS_THROUGH_TOKEN\n");
			printf ("  %4.2f\n", buffer[size-count]);
			count--;
		}
		else if (token == GL_POINT_TOKEN) {
			printf ("GL_POINT_TOKEN\n");
			print3DcolorVertex (size, &count, buffer);
		}
		else if (token == GL_LINE_TOKEN) {
			printf ("GL_LINE_TOKEN\n");
			print3DcolorVertex (size, &count, buffer);
			print3DcolorVertex (size, &count, buffer);
		}
		else if (token == GL_LINE_RESET_TOKEN) {
			printf ("GL_LINE_RESET_TOKEN\n");
			print3DcolorVertex (size, &count, buffer);
			print3DcolorVertex (size, &count, buffer);
		}
	}
}

void display(void)
{
	GLfloat feedBuffer[1024];
	GLint size;

	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	glOrtho (0.0, 100.0, 0.0, 100.0, 0.0, 1.0);

	glClearColor (0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT);
	drawGeometry (GL_RENDER);

	glFeedbackBuffer (1024, GL_3D_COLOR, feedBuffer);
	(void) glRenderMode (GL_FEEDBACK);
	drawGeometry (GL_FEEDBACK);

	size = glRenderMode (GL_RENDER);
	printBuffer (size, feedBuffer);
}

*/