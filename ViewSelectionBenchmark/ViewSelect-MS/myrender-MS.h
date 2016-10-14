#pragma once

//#include <QtCore/QObject>

#include "gl/glew.h"
#pragma comment(lib,"glew32")

#include <cg/cg.h>
#include <cg/cgGL.h>
#pragma comment(lib,"cg")
#pragma comment(lib,"cgGL")

#include <iostream>
#include <fstream>
#include <cassert>
#include <ctime>

#include "render.h"
#include "mesh.h"
#include "errorhandle.h"
#include "trackball.h"
#include "glmatrix.h"

class MyRender_MS : /*public QObject, */public Render
{
    //Q_OBJECT
public:
    // class scope definitions
    enum Image_Type {
            original_image,
            depth_buffer_image,
            model_space_curvature_image,
            radial_curvature_image,
            view_dependent_curvature_image,
            mesh_saliency
    };

    // FBO size, should be the same with GLWidget's size to be accurate
    static const int fbo_size = 512;		// 投影图像的分辨率，使用了正方形图片，分辨率越低锯齿越大 [2/18/2012 Han]

public:
    MyRender_MS(); 
    ~MyRender_MS();

    // virtual functions inherited from Render,
    // referred by GLWidget
    void initializeGL();
    void resizeGL(int width, int height);
    void paintGL();
    // for trackball manipulation
    void add_trackball_quat(float x, float y, float xx, float yy);

public /*slots*/:
    // slots for manipulation
    bool load_model(const wstring filename);
    void set_showing_image_type(Image_Type it);
	Image_Type get_showing_image_type();
    void set_adaptive_box_size(int i);
    void set_number_of_histogram_intervals(int i);

    float compute_shannon_entropy();
    float compute_shannon_entropy_II();

    float compute_revised_entropy_II();
    float compute_mesh_saliency();

    void sample_using_shannon_entropy();
    void sample_using_shannon_entropy_II();

    void sample_using_revised_entropy();
    void sample_using_revised_entropy_II();

    void sample_using_mesh_saliency();

    void testing_a_round();
    void save_result();
	//  [2/16/2012 Han]
	int GetMeshTriNum();
	point GetViewPoint();
	void GetBestViewDist(float dis, Image_Type type);	// 距离模型中心一定距离观察模型所得到的最优视点

private:
    // initiate and data processing functions
    bool init_mesh_vertex_array_pointer();
    bool init_mesh_curvature_texture();

    void init_framebuffer_object();
    void init_cg_context_and_programs();

    void build_project_matrix();

	vector<point> viewpoint_candidates;					// 保存为成员变量，便于使用
	vector<point> viewpoint_candidates_spherical;		// 保存每个候选视点的球面坐标
	void RandomDistributeViews(int numViews, float r);
	void EqualDistributeViews(int numViewsH, int numViewsV, float r);
	int best_viewpoint_idx;								// 最优视点的位置

private:    // data of MyRender_MS
    Mesh *mesh;
    Mesh *mesh_viewplane;  //this mesh is used as the view plane mesh for computing view dependent curvature

    // showing image type
    Image_Type showing_image_type;

    // read back of model space curvature projection
    // read back of radial curvature 
    // read back of projected view dependent normal, for view dependent curvature computing
    float model_space_curvature[fbo_size*fbo_size];
    float radial_curvature[fbo_size*fbo_size];
    float view_dependent_normal[fbo_size*fbo_size*3];
    float visible_mesh_saliency[fbo_size*fbo_size*3];
    
    // for computing adaptive entropy
    int mean_curvature_amplified_times;
    int adaptive_box_size;
    int number_of_histogram_intervals;
    vector<int> distribution_of_curvature;  // for entropy computing
    int curvature_counted_low_bound;
    int curvature_counted_up_bound;

    // standard deviation computing,
    // used by compute_shannon_entropy as a return value
    // then used by sampling procedure
    float current_standard_deviation;

    // used in sampling procedure, and in compute a given view's revised entropy
    float current_avg_shannon_entropy;
    float current_avg_stand_deviation;
    bool sampled_using_revised_entropy_II;


    // time consuming
    float secs;

    // mapped curvature for drawing, just for show the result
    unsigned char RGB_mapped_from_curvature[fbo_size*fbo_size*3];

    // for transform, rotation
    float trackball_quat[4];
    point eye_pos_to_bsphere_center;

    // canvas width and height
    int gl_canvas_width,gl_canvas_height;

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
};