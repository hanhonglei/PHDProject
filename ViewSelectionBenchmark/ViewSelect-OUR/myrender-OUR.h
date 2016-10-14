#pragma once


#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <stack>
#include <vector>

#include <boost/smart_ptr.hpp>
#include "mesh.h"
#include "render.h"
#include "ppm.h"
#include "volumedata.h"
#include "vec.h"
#include "trackball.h"
#include <cg/cg.h>
#include <cg/cgGL.h>
#pragma comment(lib,"cg")
#pragma comment(lib,"cgGL")

// if use whole_mesh_viewplane, or just projecting covered mesh_view_plane
#define USE_PROJECTING_COVERED_VIEWPLANE 1
#define VIEWPOINT_SAMPLES_NUM 100  // viewpoint num, too slow for large num


class Pos2D
{
public:
	int x,y;
	Pos2D(){}
	Pos2D(int x,int y):x(x),y(y){}
	int operator[](int idx)
	{
		if (idx==0) return x;
		else if (idx==1) return y;
		else return x;
	}

};

bool operator==(Pos2D const &_a,Pos2D const &_b);

class Pos2DLess
{
public:
	bool operator()(Pos2D const &_a,Pos2D const &_b) const
	{
		if (_a.x!=_b.x) return _a.x-_b.x;
		else return _a.y-_b.y;
	}

};

class Pos3D
{
public:
	int x,y,z;
	Pos3D(){}
	Pos3D(int x,int y,int z):x(x),y(y),z(z){}
	int operator[](int idx)
	{
		if (idx==0) return x;
		else if (idx==1) return y;
		else return z;
	}

};

bool operator==(Pos3D const &_a,Pos3D const &_b);


class Pos3DLess
{
public:
	bool operator()(Pos3D const &_a,Pos3D const &_b) const
	{
		if (_a.x!=_b.x) return _a.x<_b.x;
		else if (_a.y!=_b.y) return _a.y<_b.y;
		else return _a.z<_b.z;
	}
};

class ViewPoint
{
public:
	ViewPoint() {importance = 0, index = 0;}
	point pos;
	double posSphere[3];
	float color[3];
	double importance;
	int index;				// 重要度排序，0为最重要
	float entropy;
};

class ViewSphereFace
{
public:
	ViewSphereFace() {p1 = 0, p2 = 0; p3 = 0;}
	unsigned int p1,p2,p3;
} ;

point MidPoint(point p1,point p2);
void Normalise(point *p, double r);

class MyRender_OUR : public Render
{
private:
	void RandomDistributeViews(int numViews, float r);
	void EqualDistributeViews(int numViewsH, int numViewsV, float r);
	// 将paintGL函数里面的绘制函数重新封装 [3/13/2012 Han]
	void DrawModel();
	void Draw2D();
	void DrawVDCurvature();
	void DrawViewSphere();
	void DrawVoxel();
	void ViewSphereColor();

public:

	vector<ViewSphereFace> vf ;
	Mesh* GetViewplane();
	int best_viewpoint_idx;								// 最优视点的位置
	vector<point> viewpoint_candidates_spherical;		// 保存每个候选视点的球面坐标
	vector<ViewPoint*> viewpoint_candidates; // candidate on the real viewing sphere
	vector<point> meta_viewpoint_candidates; // meta data, candidates on a unit sphere
	void ResetEyePos();
	bool IsOccluded(float p[3], GLfloat *bufferZ);

	double LUT_CourbureClust[3*256]; 

	enum Image_Type{
		SHOWING_WHOLE_SKELETON, 
		SHOWING_PRIMARY_SKELETON,
		SHOWING_FINAL_SKELETON,
		SHOWING_MESH_SALIENCY,
		SHOWING_MODEL,
		SHOWING_MODEL_VOXELS,
		SHOWING_SEGMENTATION,
		SHOWING_VIEW_DEPENDENT_CURVATURE,
		SHOWING_VIEW_SPHERE_MAP,
		VISIBLE_TEST
	} showing_type;

	int GetMeshTriNum();
	point GetViewPoint();
	void GetBestViewDist(float dis, Image_Type type, FILE *outPutFile = NULL);	// 距离模型中心一定距离观察模型所得到的最优视点
	float viewpoint_selection_saliency(bool bMultiThread = false);			// Get Best View using mesh saliency
	void GetVisibleTris(vector<int> &tris);									// Get all visible tris in current view pos

    struct Skl2DNode{
        Pos2D pos;
        bool flag; // for travel 
        std::map<Pos2D, unsigned int, Pos2DLess> children;
        std::map<Pos2D, Skl2DNode*, Pos2DLess> children_pos;

        Skl2DNode(){flag = false; pos = Pos2D(0,0);}
        ~Skl2DNode()
        {
            for (std::map<Pos2D, Skl2DNode*, Pos2DLess>::iterator i = children_pos.begin();
                i != children_pos.end(); i++)
            {
                delete  i->second;
            }
        }
    };

    struct Skl3DNode{
        Pos3D pos;
        bool flag; // for travel 
        std::map<Pos3D, unsigned int, Pos3DLess> children;
        std::map<Pos3D, Skl3DNode*, Pos3DLess> children_pos;

        Skl3DNode(){flag = false; pos = Pos3D(0,0,0);}
        ~Skl3DNode()
        {
            for (std::map<Pos3D, Skl3DNode*, Pos3DLess>::iterator i = children_pos.begin();
                i != children_pos.end(); i++)
            {
                delete  i->second;
            }
        }
    };

    MyRender_OUR();
    ~MyRender_OUR();

    // called by Qt GLWidget
    void initializeGL();
    void resizeGL(int width, int height);
    void paintGL();

    // handling 2d or 3d condition
    bool model_loaded;
    enum {D2, D3} dimension;


	float trackball_quat[4];

    /************************************************************************/
    /*  2d models                                                           */
    /************************************************************************/
    boost::shared_ptr<PPM> model2d;
    boost::shared_array<float> potential2d;

    bool skeleton_computed_2d;
    std::vector<std::vector<Pos2D> > primary_skls2d;
    Skl2DNode skls_tree2d;
    std::vector<std::vector<Pos2D> > final_skls2d;
    std::vector<std::vector<Pos2D> > segmentation2d;

    boost::shared_array<float> drawing2d;

    bool load_model_2d(std::wstring fn);
    bool compute_potential_2d();
    bool extract_primary_skeleton_2d();
    bool build_skeleton_tree_2d();
    bool extract_final_skeleton_2d();
    bool segmentation_by_skeleton_2d();

    unsigned int  WOS_2D;
    void set_WOS_2D(unsigned int w){ WOS_2D = w;}
    unsigned int showing_primary_skeleton_number_2d;
    void set_showing_primary_skeleton_number_2d(unsigned int sps){showing_primary_skeleton_number_2d = sps;}
    unsigned int showing_final_skeleton_number_2d;
    void set_showing_final_skeleton_number_2d(unsigned int sfs){showing_final_skeleton_number_2d = sfs;}

    /************************************************************************/
    /*  3d models                                                           */
    /************************************************************************/
    float model_matrix[16];
    float view_matrix[16];
    float model_view_matrix[16];
    float project_matrix[16];
    float model_view_project_matrix[16];

    boost::shared_ptr<Mesh> mesh3d;                 // original mesh data
    boost::shared_ptr<VolumeData> model3d;          // voxelized volume data
    boost::shared_array<float> potential3d;

    std::vector<std::vector<Pos3D> > primary_skls3d;
    Skl3DNode skls_tree3d;
    std::vector<std::vector<Pos3D> > final_skls3d;
    std::vector<std::vector<Pos3D> > segmentation3d; // segmentation of seeds, seed are inner surface of voxelized model
    std::vector<std::vector<Pos3D> > model_seg3d; // segmentation of voxelized model;
    std::vector<Mesh*> mesh_seg3d; // segmentation of mesh

    std::wstring model_filename_3d;
    bool load_model_3d(std::wstring fn);
    bool load_segmentation(std::wstring fn);
	void show_A();
	void show_B1();
	void show_B2();
	void show_B3();
    // potential compute
    bool compute_potential_3d();
    bool save_potential_3d(std::wstring fn);
    bool load_potential_3d(std::wstring fn);

    bool extract_primary_skeleton_3d();
    bool build_skeleton_tree_3d();
    bool extract_final_skeleton_3d();
    bool segmentation_by_skeleton_3d();
    bool map_to_mesh_segmentation();

    unsigned int  WOS_3D;
    void set_WOS_3D(unsigned int w){ WOS_3D = w;}
    unsigned int showing_primary_skeleton_number_3d;
    void set_showing_primary_skeleton_number_3d(unsigned int sps){showing_primary_skeleton_number_3d = sps;}
    unsigned int showing_final_skeleton_number_3d;
    void set_showing_final_skeleton_number_3d(unsigned int sfs){showing_final_skeleton_number_3d = sfs;}

    // for 3d objects view
    void add_trackball_quat(float x, float y, float xx, float yy);
    void reset_trackball();

    // mesh saliency
    bool compute_mesh_saliency_3d();
    
    // viewpoint selection
    void build_segmentation(std::vector<int> &parts);   // build view selection considered faces
    float compute_entropy();
    float compute_revised_entropy();
    void viewpoint_selection(bool bMultiThread = false);
	float viewpoint_selection2(bool bMultiThread = false);
    void testing_a_round();
    void generate_viewpoint_candidates();

    void save_viewpoint();

    void sample_using_revised_entropy_II();
    float compute_shannon_entropy_II();

	// 这些数值都可以放入view类中，不必重新申请成员变量 [5/16/2012 Han]
    float shannon_entropy[VIEWPOINT_SAMPLES_NUM];
    float standard_deviations[VIEWPOINT_SAMPLES_NUM];
    float revised_entropy[VIEWPOINT_SAMPLES_NUM];

    // for selecting multiple vies
    Mesh *mesh_viewplane;  //this mesh is used as the view plane mesh for computing view dependent curvature
    Mesh *mesh_viewing_shpere;
    static const int fbo_size = 512;
    float view_dependent_normal[fbo_size*fbo_size*3];
    int mesh_viewplane_idx[fbo_size*fbo_size];
    bool selecting_parts_by_mouse;
    void change_mouse_mode(int mode);
    unsigned char parts_color_buffer[fbo_size*fbo_size*3];
    unsigned char parts_color_for_rendering[fbo_size*fbo_size*3];
    void mouseDrag(int startx, int starty, int endx, int endy, int state=0);

    boost::shared_ptr<Mesh> whole_mesh_buff;                 // original mesh data
    
    // for computing adaptive entropy
    int adaptive_box_size;
    int number_of_histogram_intervals;
    vector<int> distribution_of_curvature;  // for entropy computing

    // standard deviation computing,
    // used by compute_shannon_entropy as a return value
    // then used by sampling procedure
    float current_standard_deviation;

    // used in sampling procedure, and in compute a given view's revised entropy
    float current_avg_shannon_entropy;
    float current_avg_stand_deviation;

    // mapped curvature for drawing, just for show the result
    unsigned char RGB_mapped_from_curvature[fbo_size*fbo_size*3];

    // for transform, rotation
    point eye_pos_to_bsphere_center;

    // for time consuming
    float secs;

    // FBO common rendering
    unsigned int fbo_common_displaying;
    unsigned int texture_common_color_buffer;
    unsigned int texture_common_depth_buffer;

    // FOB model space curvature
    unsigned int fbo_model_space_curvature;
    unsigned int rbo_model_space_curvature_depth_buffer;
    unsigned int texture_model_space_curvature;

    // FBO radial curvature
    unsigned int fbo_radial_curvature;
    unsigned int rbo_radial_curvature_depth_buffer;
    unsigned int texture_radial_curvature;

    // FBO view dependent curvature
    unsigned int fbo_view_dependent_curvature;
    unsigned int rbo_view_dependent_curvature_depth_buffer;
    unsigned int texture_view_dependent_curvature;

    // Cg context and profile
    CGcontext cg_context;
    CGprofile cg_vprofile;
    CGprofile cg_fprofile;

    // following textures are used as input for model space curvature and radial curvature
    unsigned int texture_primary_curvature_1_direction_and_value;
    unsigned int texture_primary_curvature_2_direction_and_value;
    int          width_of_primary_curvature_texture;

    // for model space curvature Cg programs
    CGprogram cg_vprogram_model_space_curvature;
    CGprogram cg_fprogram_model_space_curvature;

    CGparameter cgparam_model_view_project_matrix_model_space_curvature;

    CGparameter cgparam_texture_pdircur1_model_space_curvature;
    CGparameter cgparam_texture_pdircur2_model_space_curvature;
    CGparameter cgparam_texture_width_model_space_curvature;

    // for radial curvature Cg programs
    CGprogram cg_vprogram_radial_curvature;
    CGprogram cg_fprogram_radial_curvature;

    CGparameter cgparam_model_view_project_matrix_radial_curvature;
    CGparameter cgparam_model_view_matrix_radial_curvature;

    CGparameter  cgparam_texture_pdircur1_radial_curvature;
    CGparameter  cgparam_texture_pdircur2_radial_curvature;
    CGparameter  cgparam_texture_width_radial_curvature;

    // for view dependent curvature Cg programs
    CGprogram cg_vprogram_view_dependent_curvature;
    CGprogram cg_fprogram_view_dependent_curvature;

    CGparameter cgparam_model_view_project_matrix_view_dependent_curvature;

private:
    // initiate and data processing functions
    bool init_mesh_vertex_array_pointer();
    bool init_mesh_curvature_texture();

    void init_framebuffer_object();
    void init_cg_context_and_programs();
};