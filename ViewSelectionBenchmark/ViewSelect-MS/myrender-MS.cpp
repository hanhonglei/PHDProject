#include "myrender-MS.h"
#include "../Common/com.h"
const float PI = 3.141592654f;
template<class T>
void save_result(T idx)
{
	//ust change the packing to ensure no overruns!
	glPixelStorei(GL_PACK_ALIGNMENT, 1);

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
	// restore default alignment
	glPixelStorei(GL_PACK_ALIGNMENT, 4);
}

MyRender_MS::MyRender_MS()
{
    // init mesh to Null
    mesh = 0;
    // init showing image type
    showing_image_type = original_image;

    // set curvature counted bound
    curvature_counted_low_bound = 0;
    curvature_counted_up_bound = 255;

    mean_curvature_amplified_times = 1000;
    current_standard_deviation = 0;
    adaptive_box_size = 1;
    number_of_histogram_intervals = 512;

    // init view dependent curvature view plane mesh,
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
    mesh_viewplane->compute_bsphere();
    // resize normal vector, and for used by render_scene procedure
    mesh_viewplane->normals.resize(fbo_size*fbo_size);
}

MyRender_MS::~MyRender_MS()
{
    delete mesh;
    delete mesh_viewplane;

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

void MyRender_MS::initializeGL()
{
    //init GLEW lib
    glewInit();

    // check for OpenGL extensions and shader supporting
    if (glDeleteFramebuffersEXT == 0 ||
        glDeleteRenderbuffersEXT == 0 ||
        glBindFramebufferEXT == 0 ||
        glBindRenderbufferEXT == 0 ||
        glGenFramebuffersEXT == 0 ||
        glGenRenderbuffersEXT == 0 ||
        glFramebufferRenderbufferEXT == 0 ||
        glFramebufferTexture2DEXT == 0 ||
        glRenderbufferStorageEXT == 0)
    {
        std::cout<<"the graphic card can not support the features needed, or using a wrong driver.\n";
    }

    static const GLfloat light0_color[4] = { 0.8f, 0.8f, 0.8f, 1.0f };
    static const GLfloat light1_color[4] = { 0.4f, 0.4f, 1.0f, 1.0f };
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

    // init frame buffer and Cg programs
    init_framebuffer_object();
    init_cg_context_and_programs();
}

void MyRender_MS::resizeGL(int width, int height)
{
    gl_canvas_width = width;
    gl_canvas_height = height;

    if (!mesh)return;
    build_project_matrix();
}

void MyRender_MS::paintGL()
{
    // check if any error occured before rendering process.
    check_OpenGL_errors("Before rendering");

    // if no mesh, then clear buffers and return
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    if(!mesh) return;

    // set to fbo viewport
    //glViewport(0,0,fbo_size,fbo_size);

    // build the matrices
    float trans1[16],rot1[16];  // translate and rotate matrix
    buildTranslateMatrix(
        -mesh->bsphere.center[0],
        -mesh->bsphere.center[1],
        -mesh->bsphere.center[2],
        trans1);
    sgi_trackball_space::build_row_major_rotmatrix(rot1,trackball_quat);

    //set viewpoint according to auto sampling parameter
    //eye_pos_to_bsphere_center = point(mesh->bsphere.r,0,-mesh->bsphere.r);
    buildLookAtMatrix(
        0,0,0,
        eye_pos_to_bsphere_center[0],
        eye_pos_to_bsphere_center[1],
        eye_pos_to_bsphere_center[2],
        0,1,0,
        view_matrix);
    multMatrix(model_matrix, rot1,trans1);
    multMatrix(model_view_matrix, view_matrix, model_matrix);
    multMatrix(model_view_project_matrix,project_matrix,model_view_matrix);

    // showing original image or depth buffer
    if (showing_image_type == original_image ||
        showing_image_type == depth_buffer_image)
    {
        // drawing the model and depth
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,fbo_common_displaying);
        glClearColor(1,1,1,1);
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        glMatrixMode(GL_PROJECTION);
        transposeMatrix(project_matrix);
        glLoadMatrixf(project_matrix);

        glMatrixMode(GL_MODELVIEW);
        transposeMatrix(model_view_matrix);
        glLoadMatrixf(model_view_matrix);

        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);

        // draw displaying textures
        glEnable(GL_LIGHTING);
        glColor3f(0.1f,0.1f,0.5f);  //set blue color
        glDrawElements(GL_TRIANGLES,(GLsizei)mesh->faces.size()*3,GL_UNSIGNED_INT,&mesh->faces[0]);

		// test 查看候选视点的位置 [2/19/2012 Han]
		//generating viewpoint candidates
		glColor3f(1.0f, 0, 0);
		glPointSize(5.0f);
		glBegin(GL_POINTS);
		for (unsigned int i = 0; i < viewpoint_candidates.size(); i++)
		{
			glVertex3f(viewpoint_candidates[i][0], viewpoint_candidates[i][1], viewpoint_candidates[i][2]);
		}
		glEnd();
		//  [2/19/2012 Han]

        glDisable(GL_LIGHTING);

        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);

        // set to system frame buffer viewport
        glViewport(0,0,gl_canvas_width,gl_canvas_height);
        // render to system frame buffer
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
        glClearColor(1,1,1,1);
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(-1,1,-1,1,-1,1);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        if (showing_image_type == original_image)
        {
            glBindTexture(GL_TEXTURE_RECTANGLE_ARB,texture_common_color_buffer);
            glColor3f(1,1,1);   // set color to original color
            glBegin(GL_QUADS);  // draw color buffer 0
            {
                glTexCoord2i(0, 0); glVertex3f(-1, -1, 0);
                glTexCoord2i(fbo_size, 0); glVertex3f(1, -1, 0);
                glTexCoord2i(fbo_size, fbo_size); glVertex3f(1, 1, 0);
                glTexCoord2i(0, fbo_size); glVertex3f(-1, 1, 0);
            }
            glEnd();
            glBindTexture(GL_TEXTURE_RECTANGLE_ARB,0);
        }
        else if (showing_image_type == depth_buffer_image)
        {

            glBindTexture(GL_TEXTURE_RECTANGLE_ARB,texture_common_depth_buffer);
            glColor3f(1,1,1);   // set color to original color
            glBegin(GL_QUADS);  // draw depth buffer
            {
                glTexCoord2i(0, 0); glVertex3f(-1, -1, 0);
                glTexCoord2i(fbo_size, 0); glVertex3f(1, -1, 0);
                glTexCoord2i(fbo_size, fbo_size); glVertex3f(1, 1, 0);
                glTexCoord2i(0, fbo_size); glVertex3f(-1, 1, 0);
            }
            glEnd();
            glBindTexture(GL_TEXTURE_RECTANGLE_ARB,0);
        }
    }
    // showing model space curvature image
    else if (showing_image_type == model_space_curvature_image)
    {
		//ust change the packing to ensure no overruns!
		glPixelStorei(GL_PACK_ALIGNMENT, 1);

        // compute radial curvature
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,fbo_model_space_curvature);
        glClearColor(0,0,0,0);
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

        cgGLBindProgram(cg_vprogram_model_space_curvature);
        cgGLEnableProfile(cg_vprofile);
        cgGLBindProgram(cg_fprogram_model_space_curvature);
        cgGLEnableProfile(cg_fprofile);
        
        cgSetMatrixParameterfr(
            cgparam_model_view_project_matrix_model_space_curvature, 
            model_view_project_matrix);
        cgGLSetTextureParameter(
            cgparam_texture_pdircur1_model_space_curvature,
            texture_primary_curvature_1_direction_and_value);
        cgGLEnableTextureParameter(cgparam_texture_pdircur1_model_space_curvature);
        
        cgGLSetTextureParameter(cgparam_texture_pdircur2_model_space_curvature,
            texture_primary_curvature_2_direction_and_value);
        cgGLEnableTextureParameter(cgparam_texture_pdircur2_model_space_curvature);
        
        cgSetParameter1i(cgparam_texture_width_model_space_curvature,
            width_of_primary_curvature_texture);


        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);

        glDrawElements(GL_TRIANGLES,(GLsizei)mesh->faces.size()*3,GL_UNSIGNED_INT,&mesh->faces[0]);

        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);

        // disable the profiles
        cgGLDisableProfile(cg_vprofile);
        cgGLDisableProfile(cg_fprofile);

        //read back radial curvatures, and normalization, and upload to GPU for rendering.
        glReadPixels(0,0,fbo_size,fbo_size,GL_LUMINANCE,GL_FLOAT,model_space_curvature);
        //reset the canvas background
        memset(RGB_mapped_from_curvature,96,
            3*fbo_size*fbo_size*sizeof(unsigned char));
		// restore default alignment
		glPixelStorei(GL_PACK_ALIGNMENT, 4);

        float maxc=0,minc=0;
        for (unsigned int i=0; i<fbo_size*fbo_size; i++)
        {
            if (model_space_curvature[i] > maxc)
            {
                maxc = model_space_curvature[i];
            }
            if (model_space_curvature[i] < minc)
            {
                minc = model_space_curvature[i];
            }
        }
        std::cout<<"max mc: "<<maxc<<std::endl;
        std::cout<<"min mc: "<<minc<<std::endl;

        int mstat = int(abs(maxc)>=abs(minc)?abs(maxc):abs(minc)) + 1;
        std::vector<int> histo;
        histo.resize(mstat,0);

        for (unsigned int i=0; i<fbo_size*fbo_size; i++)
        {
            int idx = int(abs(model_space_curvature[i]));
            histo[idx] += 1;
        }

        float sum = 0;
        for (unsigned int i=1; i<histo.size(); i++)
        {
            sum += histo[i];
        }

        const float RATIO = 0.999999f;
        float counted = 0;
        int upbound = 0;
        for (unsigned int i=1; i<histo.size(); i++)
        {
            counted += histo[i];
            if (counted/sum > RATIO)
            {
                upbound = i;
                break;
            }
        }
        std::cout<<"Up bound:"<<upbound<<std::endl;
        upbound = 30000;
        for (int i=0; i<fbo_size*fbo_size; i++)
        {   // map the curvatures to color space for displaying
            // just for view, so we only handle curvature < 255
            if (int(model_space_curvature[i]) > upbound)
            {
                RGB_mapped_from_curvature[i*3+0] = 255;
                RGB_mapped_from_curvature[i*3+1] = 0;
                RGB_mapped_from_curvature[i*3+2] = 0;
            }
            else if (model_space_curvature[i] > 0 && int(model_space_curvature[i]) <= upbound)
            {
                RGB_mapped_from_curvature[i*3+0] = 255;
                RGB_mapped_from_curvature[i*3+1] = 255 - int(model_space_curvature[i]/upbound*255);
                RGB_mapped_from_curvature[i*3+2] = 255 - int(model_space_curvature[i]/upbound*255);
            }
            else if (model_space_curvature[i] < 0 && int(model_space_curvature[i]) >= -upbound)
            {
                RGB_mapped_from_curvature[i*3+0] = 255 + int(model_space_curvature[i]/upbound*255);
                RGB_mapped_from_curvature[i*3+1] = 255;
                RGB_mapped_from_curvature[i*3+2] = 255 + int(model_space_curvature[i]/upbound*255);
            }
            else if (int(model_space_curvature[i]) < -upbound)
            {
                RGB_mapped_from_curvature[i*3+0] = 0;
                RGB_mapped_from_curvature[i*3+1] = 255;
                RGB_mapped_from_curvature[i*3+2] = 0;
            }
        }

        // set to system framebuffer viewport
        glViewport(0,0,gl_canvas_width, gl_canvas_height);
        // render to system framebuffer
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
        glDrawPixels(fbo_size,fbo_size,GL_RGB,GL_UNSIGNED_BYTE,RGB_mapped_from_curvature);
    }
    // showing radial curvature image
    else if (showing_image_type == radial_curvature_image)
    {
		//ust change the packing to ensure no overruns!
		glPixelStorei(GL_PACK_ALIGNMENT, 1);

        // compute radial curvature
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,fbo_radial_curvature);
        glClearColor(0,0,0,0);
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

        cgGLBindProgram(cg_vprogram_radial_curvature);
        cgGLEnableProfile(cg_vprofile);
        cgGLBindProgram(cg_fprogram_radial_curvature);
        cgGLEnableProfile(cg_fprofile);
        
        cgSetMatrixParameterfr(
            cgparam_model_view_project_matrix_radial_curvature, 
            model_view_project_matrix);
        cgSetMatrixParameterfr(
            cgparam_model_view_matrix_radial_curvature,
            model_view_matrix);
        
        cgGLSetTextureParameter(
            cgparam_texture_pdircur1_radial_curvature,
            texture_primary_curvature_1_direction_and_value);
        cgGLEnableTextureParameter(cgparam_texture_pdircur1_radial_curvature);
        
        cgGLSetTextureParameter(
            cgparam_texture_pdircur2_radial_curvature,
            texture_primary_curvature_2_direction_and_value);
        cgGLEnableTextureParameter(cgparam_texture_pdircur2_radial_curvature);
        
        cgSetParameter1i(
            cgparam_texture_width_radial_curvature,
            width_of_primary_curvature_texture);


        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);

        glDrawElements(GL_TRIANGLES,(GLsizei)mesh->faces.size()*3,GL_UNSIGNED_INT,&mesh->faces[0]);

        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);

        // disable the profiles
        cgGLDisableProfile(cg_vprofile);
        cgGLDisableProfile(cg_fprofile);

        //read back radial curvatures, and normalization, and upload to GPU for rendering.
        glReadPixels(0,0,fbo_size,fbo_size,GL_LUMINANCE,GL_FLOAT,radial_curvature);
        //reset the canvas background
        memset(RGB_mapped_from_curvature,96,
            3*fbo_size*fbo_size*sizeof(unsigned char));
		// restore default alignment
		glPixelStorei(GL_PACK_ALIGNMENT, 4);

        float maxc=0,minc=0;
        for (unsigned int i=0; i<fbo_size*fbo_size; i++)
        {
            if (radial_curvature[i] > maxc)
            {
                maxc = radial_curvature[i];
            }
            if (radial_curvature[i] < minc)
            {
                minc = radial_curvature[i];
            }
        }
        std::cout<<"max mc: "<<maxc<<std::endl;
        std::cout<<"min mc: "<<minc<<std::endl;

        int mstat = int(abs(maxc)>=abs(minc)?abs(maxc):abs(minc)) + 1;
        std::vector<int> histo;
        histo.resize(mstat,0);

        for (unsigned int i=0; i<fbo_size*fbo_size; i++)
        {
            int idx = int(abs(radial_curvature[i]));
            histo[idx] += 1;
        }

        float sum = 0;
        for (unsigned int i=1; i<histo.size(); i++)
        {
            sum += histo[i];
        }

        const double RATIO = 1 - 1e-1;
        float counted = 0;
        int upbound = 0;
        for (unsigned int i=1; i<histo.size(); i++)
        {
            counted += histo[i];
            if (counted/sum > RATIO)
            {
                upbound = i;
                break;
            }
        }
        upbound = 255;
        for (int i=0; i<fbo_size*fbo_size; i++)
        {   // map the curvatures to color space for displaying
            // just for view, so we only handle curvature < 255
            if (int(radial_curvature[i]) > upbound)
            {
                RGB_mapped_from_curvature[i*3+0] = 255;
                RGB_mapped_from_curvature[i*3+1] = 0;
                RGB_mapped_from_curvature[i*3+2] = 0;
            }
            else if (radial_curvature[i] > 0 && int(radial_curvature[i]) <= upbound)
            {
                RGB_mapped_from_curvature[i*3+0] = 255;
                RGB_mapped_from_curvature[i*3+1] = 255 - int(radial_curvature[i]/upbound*255);
                RGB_mapped_from_curvature[i*3+2] = 255 - int(radial_curvature[i]/upbound*255);
            }
            else if (radial_curvature[i] < 0 && int(radial_curvature[i]) >= -upbound)
            {
                RGB_mapped_from_curvature[i*3+0] = 255 + int(radial_curvature[i]/upbound*255);
                RGB_mapped_from_curvature[i*3+1] = 255;
                RGB_mapped_from_curvature[i*3+2] = 255 + int(radial_curvature[i]/upbound*255);
            }
            else if (int(radial_curvature[i]) < -upbound)
            {
                RGB_mapped_from_curvature[i*3+0] = 0;
                RGB_mapped_from_curvature[i*3+1] = 255;
                RGB_mapped_from_curvature[i*3+2] = 0;
            }
        }

        //for (int i=0; i<fbo_size*fbo_size; i++)
        //{   // map the curvatures to color space for displaying
        //    // just for view, so we only handle curvature < 255
        //    if (int(radial_curvature[i]) > 255)
        //    {
        //        RGB_mapped_from_curvature[i*3+0] = 255;
        //        RGB_mapped_from_curvature[i*3+1] = 0;
        //        RGB_mapped_from_curvature[i*3+2] = 0;
        //    }
        //    else if (radial_curvature[i] > 0 && int(radial_curvature[i]) <= 255)
        //    {
        //        RGB_mapped_from_curvature[i*3+0] = 255;
        //        RGB_mapped_from_curvature[i*3+1] = 255 - int(radial_curvature[i]);
        //        RGB_mapped_from_curvature[i*3+2] = 255 - int(radial_curvature[i]);
        //    }
        //    else if (radial_curvature[i] < 0 && int(radial_curvature[i]) >= -255)
        //    {
        //        RGB_mapped_from_curvature[i*3+0] = 255 + int(radial_curvature[i]);
        //        RGB_mapped_from_curvature[i*3+1] = 255;
        //        RGB_mapped_from_curvature[i*3+2] = 255 + int(radial_curvature[i]);
        //    }
        //    else if (int(radial_curvature[i]) < -255)
        //    {
        //        RGB_mapped_from_curvature[i*3+0] = 0;
        //        RGB_mapped_from_curvature[i*3+1] = 255;
        //        RGB_mapped_from_curvature[i*3+2] = 0;
        //    }
        //}

        // set to system framebuffer viewport
        glViewport(0,0,gl_canvas_width, gl_canvas_height);
        // render to system framebuffer
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
        glDrawPixels(fbo_size,fbo_size,GL_RGB,GL_UNSIGNED_BYTE,RGB_mapped_from_curvature);
    }
    //showing view dependent curvature
    else if (showing_image_type == view_dependent_curvature_image)
    {
		//ust change the packing to ensure no overruns!
		glPixelStorei(GL_PACK_ALIGNMENT, 1);

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

        glDrawElements(GL_TRIANGLES,(GLsizei)mesh->faces.size()*3,GL_UNSIGNED_INT,&mesh->faces[0]);

        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);

        // disable the profiles
        cgGLDisableProfile(cg_vprofile);
        cgGLDisableProfile(cg_fprofile);

        //read back radial curvatures, and normalization, and upload to GPU for rendering.
        glReadPixels(0,0,fbo_size,fbo_size,GL_RGB,GL_FLOAT,view_dependent_normal);

		// restore default alignment
		glPixelStorei(GL_PACK_ALIGNMENT, 4);

        for (int i=0; i<fbo_size*fbo_size; i++)
        {
            mesh_viewplane->normals[i][0] = view_dependent_normal[i*3 + 0];
            mesh_viewplane->normals[i][1] = view_dependent_normal[i*3 + 1];
            mesh_viewplane->normals[i][2] = view_dependent_normal[i*3 + 2];
        }
        mesh_viewplane->curv1.clear();
        mesh_viewplane->curv2.clear();
        mesh_viewplane->compute_curvatures();
        mesh_viewplane->compute_averaging_mean_curvature();

        //reset the canvas background
        memset(RGB_mapped_from_curvature,96,
            3*fbo_size*fbo_size*sizeof(unsigned char));
        for (int i=0; i<fbo_size*fbo_size; i++)
        {   // map the curvatures to color space for displaying
            // just for view, so we only handle curvature < 255
            float mc = (mesh_viewplane->curv1[i] + mesh_viewplane->curv2[i])/2.0f;
            mc*=mean_curvature_amplified_times;

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

        // set to system framebuffer viewport
        glViewport(0,0,gl_canvas_width,gl_canvas_height);
        // render to system framebuffer
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
        glDrawPixels(fbo_size,fbo_size,GL_RGB,GL_UNSIGNED_BYTE,RGB_mapped_from_curvature);
    }
    else if (showing_image_type == mesh_saliency)
	{	//ust change the packing to ensure no overruns!
		glPixelStorei(GL_PACK_ALIGNMENT, 1);

        //mesh->compute_mesh_saliency();
        if (mesh->saliency.size() == 0)
        {
            return;
        }
        // drawing the model and depth
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,fbo_common_displaying);
        glClearColor(0,0,0,1.0f);
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        glMatrixMode(GL_PROJECTION);
        transposeMatrix(project_matrix);
        glLoadMatrixf(project_matrix);

        glMatrixMode(GL_MODELVIEW);
        transposeMatrix(model_view_matrix);
        glLoadMatrixf(model_view_matrix);

        //glEnable(GL_LIGHTING);
        glBegin(GL_TRIANGLES);
        glColor3f(1.0,0,0);
        for (unsigned int i=0; i<mesh->faces.size(); i++)
        {
            for (int j=0; j<3; j++)
            {
                std::vector<unsigned int>::size_type v = mesh->faces[i][j];
                glNormal3f(
                    mesh->normals[v][0],
                    mesh->normals[v][1],
                    mesh->normals[v][2]);
                glColor3f(
                    mesh->mapped_color[3*v + 0],
                    mesh->mapped_color[3*v + 1],
                    mesh->mapped_color[3*v + 2]);
                glVertex3f(
                    mesh->vertices[v][0],
                    mesh->vertices[v][1],
                    mesh->vertices[v][2]);
            }
        }
        glEnd();
        //glDisable(GL_LIGHTING);

        // read back mesh saliency
        glReadPixels(0,0,fbo_size,fbo_size,GL_RGB,GL_FLOAT,visible_mesh_saliency);

		// restore default alignment
		glPixelStorei(GL_PACK_ALIGNMENT, 4);
        // set to system frame buffer viewport
        glViewport(0,0,gl_canvas_width,gl_canvas_height);
        // render to system frame buffer
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
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
            glTexCoord2i(0, 0); glVertex3f(-1, -1, 0);
            glTexCoord2i(fbo_size, 0); glVertex3f(1, -1, 0);
            glTexCoord2i(fbo_size, fbo_size); glVertex3f(1, 1, 0);
            glTexCoord2i(0, fbo_size); glVertex3f(-1, 1, 0);
        }
        glEnd();
        glBindTexture(GL_TEXTURE_RECTANGLE_ARB,0);
    }
    check_OpenGL_errors("After Rendering");
}

bool MyRender_MS::load_model(const wstring filename)
{
    Mesh *tm;    //temp mesh;
    tm = Mesh::read(filename);
    if(tm)
    {
        delete mesh;
        mesh = tm;

        mesh->compute_bsphere();
        mesh->compute_normals();
        mesh->compute_curvatures();
        mesh->compute_averaging_mean_curvature();

        std::cout
            <<"Vertices: "<<mesh->vertices.size()<<std::endl
            <<"Faces: "<<mesh->faces.size()<<std::endl;

        init_mesh_vertex_array_pointer();
        init_mesh_curvature_texture();

        // because bsphere re-computed, so, rebuild project matrix
        build_project_matrix(); 
        //reset eye position to bounding sphere center
        eye_pos_to_bsphere_center = point(0,0,-mesh->bsphere.r);
        sgi_trackball_space::trackball(trackball_quat,0,0,0,0); // and also trackball quat

        current_avg_shannon_entropy = 0;
        current_avg_stand_deviation = 0;
        sampled_using_revised_entropy_II = false;
        return true;
    }
    else
    {
        return false;
    }
}

bool MyRender_MS::init_mesh_vertex_array_pointer()
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

bool MyRender_MS::init_mesh_curvature_texture()
{
    if (!mesh)
    {
        std::cout<<"mesh is empty, no curvature texture initialized.\n";
        return false;
    }

    mesh->compute_curvatures();
    int width = (int)floor(sqrt(float(mesh->vertices.size())))+1;
    width_of_primary_curvature_texture = width;
    float *pd = new float[4*width*width];
    if (!pd)
    {
        std::cout<<"Not enough memory, make more money, and upgrade your machine!\n";
        return false;
    }

    // use GL_ABGR_EXT format to upload to GL_RGBA_FLOAT32_ATI as vertex texture
    vector<point>::size_type N = vector<point>::size_type(width*width);
    // upload principal curvature 1
    for (vector<point>::size_type i=0; i<mesh->vertices.size(); i++)
    {
        pd[i*4 + 0] = mesh->curv1[i];
        pd[i*4 + 1] = mesh->pdir1[i][2];
        pd[i*4 + 2] = mesh->pdir1[i][1];
        pd[i*4 + 3] = mesh->pdir1[i][0];
    }
    for (vector<point>::size_type i=mesh->vertices.size(); i<N; i++)
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
    for (vector<point>::size_type i=0; i<mesh->vertices.size(); i++)
    {
        pd[i*4 + 0] = mesh->curv2[i];
        pd[i*4 + 1] = mesh->pdir2[i][2];
        pd[i*4 + 2] = mesh->pdir2[i][1];
        pd[i*4 + 3] = mesh->pdir2[i][0];
    }
    for (vector<point>::size_type i=mesh->vertices.size(); i<N; i++)
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
}

void MyRender_MS::init_framebuffer_object()
{
    std::cout<<"Init FrameBuffer and Check Status\n";
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
        fbo_size,
        fbo_size,
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
        fbo_size,
        fbo_size,
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

    check_frame_buffer_object_status("Common Display");
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);
    
    /************************************************************************/
    /*       FBO AND BUFFER FOR MODEL SPACE CURVATURE COMPUTING             */
    /************************************************************************/

    // generating fbo_computing
    glGenFramebuffersEXT(1, &fbo_model_space_curvature);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo_model_space_curvature);

    glGenTextures(1, &texture_model_space_curvature);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB,texture_model_space_curvature);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(
        GL_TEXTURE_RECTANGLE_ARB,
        0,
        GL_FLOAT_R32_NV,
        fbo_size,
        fbo_size,
        0,
        GL_LUMINANCE,
        GL_FLOAT,
        NULL);
    glFramebufferTexture2DEXT(
        GL_FRAMEBUFFER_EXT,
        GL_COLOR_ATTACHMENT0_EXT,
        GL_TEXTURE_RECTANGLE_ARB,
        texture_model_space_curvature,
        0);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB,0);

    // attach depth-buffer rbo
    glGenRenderbuffersEXT(1,&rbo_model_space_curvature_depth_buffer);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, rbo_model_space_curvature_depth_buffer);
    glRenderbufferStorageEXT(
        GL_RENDERBUFFER_EXT,
        GL_DEPTH_COMPONENT24,
        fbo_size,
        fbo_size);
    glFramebufferRenderbufferEXT(
        GL_FRAMEBUFFER_EXT,
        GL_DEPTH_ATTACHMENT_EXT,
        GL_RENDERBUFFER_EXT,
        rbo_model_space_curvature_depth_buffer);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT,0);

    check_frame_buffer_object_status("Model Space Curvature");
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);



    /************************************************************************/
    /*       FBO AND BUFFER FOR RADIAL CURVATURE COMPUTING                  */
    /************************************************************************/

    // generating fbo_computing
    glGenFramebuffersEXT(1, &fbo_radial_curvature);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo_radial_curvature);

    // attach curvature texture
    glGenTextures(1,&texture_radial_curvature);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB,texture_radial_curvature);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage2D(
        GL_TEXTURE_RECTANGLE_ARB,
        0,
        GL_FLOAT_R32_NV,
        fbo_size,
        fbo_size,
        0,
        GL_LUMINANCE,
        GL_FLOAT,
        NULL);
    glFramebufferTexture2DEXT(
        GL_FRAMEBUFFER_EXT,
        GL_COLOR_ATTACHMENT0_EXT,
        GL_TEXTURE_RECTANGLE_ARB,
        texture_radial_curvature,
        0);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB,0);

    // attach depth-buffer rbo
    glGenRenderbuffersEXT(1,&rbo_radial_curvature_depth_buffer);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, rbo_radial_curvature_depth_buffer);
    glRenderbufferStorageEXT(
        GL_RENDERBUFFER_EXT,
        GL_DEPTH_COMPONENT24,
        fbo_size,
        fbo_size);
    glFramebufferRenderbufferEXT(
        GL_FRAMEBUFFER_EXT,
        GL_DEPTH_ATTACHMENT_EXT,
        GL_RENDERBUFFER_EXT,
        rbo_radial_curvature_depth_buffer);
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT,0);

    check_frame_buffer_object_status("Radial Curvature");
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);

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

    check_frame_buffer_object_status("View Dependent Curvature");
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT,0);

    std::cout<<"End\n\n";
}

void MyRender_MS::init_cg_context_and_programs()
{
    std::cout<<"Load Cg programs and Check status\n";
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
    /*           Load Cg Programs                                           */
    /************************************************************************/

    /************************************************************************/
    /*           Load Cg Programs For Model Space Curvature Computing       */
    /************************************************************************/

    // vertex program for radial curvature
    cg_vprogram_model_space_curvature = cgCreateProgramFromFile(
        cg_context,
        CG_SOURCE,
        "modelspacecurvature_v.cg",
        cg_vprofile,
        NULL,NULL);
    std::cout
        <<"Compile and Load Model Space Curvature Computing Vertex Program\nLAST LISTING----"
        <<(cgGetLastListing(cg_context)==NULL?"None":cgGetLastListing(cg_context))
        <<"----\n";
    /*
    std::cout
        <<"---- PROGRAM BEGIN ----\n"
        <<string(cgGetProgramString(cg_vprogram_model_space_curvature, CG_COMPILED_PROGRAM))
        <<"---- PROGRAM END ----\n\n";
    //*/
    cgGLLoadProgram(cg_vprogram_model_space_curvature);
    cgparam_model_view_project_matrix_model_space_curvature =
        cgGetNamedParameter(cg_vprogram_model_space_curvature, "model_view_project");
    cgparam_texture_pdircur1_model_space_curvature = 
        cgGetNamedParameter(cg_vprogram_model_space_curvature, "texture_pdircur1");
    cgparam_texture_pdircur2_model_space_curvature = 
        cgGetNamedParameter(cg_vprogram_model_space_curvature, "texture_pdircur2");
    cgparam_texture_width_model_space_curvature = 
        cgGetNamedParameter(cg_vprogram_model_space_curvature, "texture_width");

    // fragment program for radial curvature
    cg_fprogram_model_space_curvature =
        cgCreateProgramFromFile(
        cg_context,
        CG_SOURCE,
        "modelspacecurvature_f.cg",
        cg_fprofile,
        NULL,NULL);
    std::cout
        <<"Compile and Load Model Space Curvature Computing Fragment Program\nLAST LISTING----"
        <<(cgGetLastListing(cg_context)==NULL?"None":cgGetLastListing(cg_context))
        <<"----\n";
    /*
    std::cout
        <<"---- PROGRAM BEGIN ----\n"
        <<cgGetProgramString(cg_fprogram_radial_curvature, CG_COMPILED_PROGRAM)
        <<"---- PROGRAM END ----\n\n";
    //*/
    cgGLLoadProgram(cg_fprogram_model_space_curvature);

    /************************************************************************/
    /*           Load Cg Programs For Radial Curvature Computing            */
    /************************************************************************/

    // vertex program for radial curvature
    cg_vprogram_radial_curvature = cgCreateProgramFromFile(
        cg_context,
        CG_SOURCE,
        "radialcurvature_v.cg",
        cg_vprofile,
        NULL,NULL);
    std::cout
        <<"Compile and Load Radial Curvature Computing Vertex Program\nLAST LISTING----"
        <<(cgGetLastListing(cg_context)==NULL?"None":cgGetLastListing(cg_context))
        <<"----\n";
    /*
    std::cout
        <<"---- PROGRAM BEGIN ----\n"
        <<string(cgGetProgramString(cg_vprogram_radial_curvature, CG_COMPILED_PROGRAM))
        <<"---- PROGRAM END ----\n\n";
    //*/
    cgGLLoadProgram(cg_vprogram_radial_curvature);
    cgparam_model_view_project_matrix_radial_curvature =
        cgGetNamedParameter(cg_vprogram_radial_curvature, "model_view_project");
    cgparam_model_view_matrix_radial_curvature =
        cgGetNamedParameter(cg_vprogram_radial_curvature, "model_view");
    cgparam_texture_pdircur1_radial_curvature = 
        cgGetNamedParameter(cg_vprogram_radial_curvature, "texture_pdircur1");
    cgparam_texture_pdircur2_radial_curvature = 
        cgGetNamedParameter(cg_vprogram_radial_curvature, "texture_pdircur2");
    cgparam_texture_width_radial_curvature = 
        cgGetNamedParameter(cg_vprogram_radial_curvature, "texture_width");

    // fragment program for radial curvature
    cg_fprogram_radial_curvature =
        cgCreateProgramFromFile(
        cg_context,
        CG_SOURCE,
        "radialcurvature_f.cg",
        cg_fprofile,
        NULL,NULL);
    std::cout
        <<"Compile and Load Radial Curvature Computing Fragment Program\nLAST LISTING----"
        <<(cgGetLastListing(cg_context)==NULL?"None":cgGetLastListing(cg_context))
        <<"----\n";
    /*
    std::cout
        <<"---- PROGRAM BEGIN ----\n"
        <<cgGetProgramString(cg_fprogram_radial_curvature, CG_COMPILED_PROGRAM)
        <<"---- PROGRAM END ----\n\n";
    //*/
    cgGLLoadProgram(cg_fprogram_radial_curvature);

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

    std::cout<<"End\n\n";
}

void MyRender_MS::add_trackball_quat(float x, float y, float xx, float yy)
{
    //add a new rotation operation to this trackball
    float quat[4];
    sgi_trackball_space::trackball(quat,x,y,xx,yy);
    sgi_trackball_space::add_quats(quat, trackball_quat, trackball_quat);
}

void MyRender_MS::build_project_matrix()
{
    // build the projection matrix
    buildOrthoMatrix(
        -mesh->bsphere.r,
        mesh->bsphere.r,
        -mesh->bsphere.r,
        mesh->bsphere.r,
        -mesh->bsphere.r,
        mesh->bsphere.r,
        project_matrix);
}

void MyRender_MS::set_adaptive_box_size(int i)
{
    adaptive_box_size = i;
}

void MyRender_MS::set_number_of_histogram_intervals(int i)
{
    number_of_histogram_intervals = i;
}

void MyRender_MS::set_showing_image_type(Image_Type it)
{
    showing_image_type = it;
}

MyRender_MS::Image_Type MyRender_MS::get_showing_image_type()
{
	return showing_image_type ;
}

float MyRender_MS::compute_shannon_entropy()
{
    if (!mesh)
    {
        std::cout<<"plz load a mesh first....\n";
        return -1;
    }
    if (showing_image_type != model_space_curvature_image &&
        showing_image_type != radial_curvature_image &&
        showing_image_type != view_dependent_curvature_image)
    {
        std::cout<<"plz select a correct image type first\n";
        return -1;
    }

    // initiate statistic of curvature 
    distribution_of_curvature.clear();
    if (curvature_counted_low_bound == 0)
    {
        distribution_of_curvature.resize(
            2*(curvature_counted_up_bound - curvature_counted_low_bound) + 1);
    } else{
        distribution_of_curvature.resize(
            2*(curvature_counted_up_bound - curvature_counted_low_bound) + 2);
    }

    // use 'as' for abbreviation
    int as = adaptive_box_size;

    // statistic of curvature
    if (showing_image_type == view_dependent_curvature_image)
    {
        for (int j=0; j<fbo_size; j+=as)
        {
            for (int i=0; i<fbo_size; i+=as)
            {
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
                mc*=mean_curvature_amplified_times; // amplify it for computing
                mc/=as*as; // get average mean curvature over an adaptive box

                // add to curvature statistic vector
                if (int(mc) >= -curvature_counted_up_bound &&
                    int(mc) <= -curvature_counted_low_bound)
                {
                    ++distribution_of_curvature[int(mc) + curvature_counted_up_bound];
                }
                else if (int(mc) >= curvature_counted_low_bound &&
                    int(mc) <= curvature_counted_up_bound)
                {
                    if (curvature_counted_low_bound == 0)
                    {
                        ++distribution_of_curvature[int(mc) + curvature_counted_up_bound];
                    }else
                    {
                        ++distribution_of_curvature[int(mc) - curvature_counted_low_bound
                            + curvature_counted_up_bound - curvature_counted_low_bound + 1];
                    } 
                }
            }// for (int i=0;...)
        }// for (int j=0;...)
    }
    else if (showing_image_type == model_space_curvature_image)
    {
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
                        mc  +=  model_space_curvature[(j+m)*fbo_size+i+n];
                    }
                }
                mc/=as*as; // get average mean curvature over an adaptive box

                // add to curvature statistic vector
                if (int(mc) >= -curvature_counted_up_bound &&
                    int(mc) <= -curvature_counted_low_bound)
                {
                    ++distribution_of_curvature[int(mc) + curvature_counted_up_bound];
                }
                else if (int(mc) >= curvature_counted_low_bound &&
                    int(mc) <= curvature_counted_up_bound)
                {
                    if (curvature_counted_low_bound == 0)
                    {
                        ++distribution_of_curvature[int(mc) + curvature_counted_up_bound];
                    }else
                    {
                        ++distribution_of_curvature[int(mc) - curvature_counted_low_bound
                            + curvature_counted_up_bound - curvature_counted_low_bound + 1];
                    } 
                }
            }// for (int i=0;...)
        }// for (int j=0;...)
    }
    else if (showing_image_type == radial_curvature_image)
    {
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
                        mc  += radial_curvature[(j+m)*fbo_size+i+n];
                    }
                }
                mc/=as*as; // get average mean curvature over an adaptive box

                // add to curvature statistic vector
                if (int(mc) >= -curvature_counted_up_bound &&
                    int(mc) <= -curvature_counted_low_bound)
                {
                    ++distribution_of_curvature[int(mc) + curvature_counted_up_bound];
                }
                else if (int(mc) >= curvature_counted_low_bound &&
                    int(mc) <= curvature_counted_up_bound)
                {
                    if (curvature_counted_low_bound == 0)
                    {
                        ++distribution_of_curvature[int(mc) + curvature_counted_up_bound];
                    }else
                    {
                        ++distribution_of_curvature[int(mc) - curvature_counted_low_bound
                            + curvature_counted_up_bound - curvature_counted_low_bound + 1];
                    } 
                }
            }// for (int i=0;...)
        }// for (int j=0;...)
    }

    // computing standard deviation
    float avg_mc=0; // average mean curvature,
    for (vector<int>::size_type i=0; i<distribution_of_curvature.size(); i++)
    {
        avg_mc += distribution_of_curvature[i];
    }
    avg_mc/=distribution_of_curvature.size();
    float xigma_mc = 0;
    for (vector<int>::size_type i=0; i<distribution_of_curvature.size(); i++)
    {
        xigma_mc += (avg_mc - distribution_of_curvature[i])*(avg_mc - distribution_of_curvature[i]);
    }
    xigma_mc /= distribution_of_curvature.size() - 1;
    current_standard_deviation = sqrt(xigma_mc);
    current_standard_deviation/=mean_curvature_amplified_times; 

    // computing entropy 
    float N=0;
    for (vector<int>::size_type i=0; i<distribution_of_curvature.size(); i++)
    {
        N += distribution_of_curvature[i];
    }
    float pi,E=0;
    for (vector<int>::size_type i=0; i<distribution_of_curvature.size(); i++)
    {
        if (distribution_of_curvature[i] != 0)
        {
            pi = distribution_of_curvature[i]/N;
            E += -pi*log(pi)/log(2.0f);
        }
    }

    std::cout<<"Shannon Entropy: "<<E<<std::endl;
    std::cout<<"Standard Deviation : "<<current_standard_deviation<<std::endl;
    return E;
}

float MyRender_MS::compute_shannon_entropy_II()
{
    if (!mesh)
    {
        std::cout<<"plz load a mesh first....\n";
        return -1;
    }
    if (showing_image_type != view_dependent_curvature_image)
    {
        std::cout<<"plz select a correct image type first\n";
        return -1;
    }

    // curvatures for computing entropy
    std::vector<float> curvatures;
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
    //std::cout<<"max mc: "<<maxc<<std::endl;
    //std::cout<<"min mc: "<<minc<<std::endl;

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

    std::cout<<"Shannon Entropy II: "<<E<<std::endl;
    std::cout<<"Standard Deviation II: "<<current_standard_deviation<<std::endl;
    return E;
}
float MyRender_MS::compute_revised_entropy_II()
{
    if (!sampled_using_revised_entropy_II)
    {
        std::cout<<"Sample it to get average entropy and stand deviation first...\n";
        return -1;
    }

    float e = compute_shannon_entropy_II();
    float s = current_standard_deviation;

    float re = e - 3*abs(s - current_avg_stand_deviation)
        *current_avg_shannon_entropy/current_avg_stand_deviation;

    std::cout<<"Revised Entropy II: "<<re<<std::endl;
    return re;
}

float MyRender_MS::compute_mesh_saliency()
{
    if (!mesh)
    {
        std::cout<<"plz load a mesh first....\n";
        return -1;
    }
    if (showing_image_type != mesh_saliency)
    {
        std::cout<<"plz select a correct image type first\n";
        return -1;
    }
    double sum = 0;
    for (unsigned int i=0; i<fbo_size*fbo_size*3; i++)
    {
        sum += visible_mesh_saliency[i];
    }
    std::cout<<"Mesh Saliency: "<<sum<<std::endl;
    return float(sum);
}

void MyRender_MS::sample_using_shannon_entropy()
{
#define VIEWPOINT_SAMPLES_NUM 200 
    if (!mesh)
    {
        std::cout<<"plz load a mesh first....\n";
        return;
    }
    if (showing_image_type != model_space_curvature_image &&
        showing_image_type != radial_curvature_image &&
        showing_image_type != view_dependent_curvature_image)
    {
        std::cout<<"plz select a correct image type first\n";
        return;
    }

    //generating viewpoint candidates
    viewpoint_candidates.clear();
    float r = mesh->bsphere.r;   //radius of bounding sphere
    time_t tm;
    srand((unsigned)time(&tm));  //set seed;
    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        float alpha = (rand()%360)/360.0f*2*PI;        //angle with positive-x-axis
        float beta = (rand()%360)/360.0f*2*PI - PI;  //angle with xy-plane
        point p;
        p[0] = r*cos(beta)*sin(alpha);
        p[1] = r*cos(beta)*cos(alpha);
        p[2] = r*sin(beta);
        viewpoint_candidates.push_back(p);
    }

    //reset trackball quat, time accumulator
    sgi_trackball_space::trackball(trackball_quat,0,0,0,0);

    clock_t start_t = clock();
    float shannon_entropy[VIEWPOINT_SAMPLES_NUM];

    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        //set auto_sampling_veiwpoint and compute the entropy
        eye_pos_to_bsphere_center = viewpoint_candidates[i];
        paintGL();
        shannon_entropy[i] = compute_shannon_entropy();
    }

    // find max Shannon entropy 
    float max_E = shannon_entropy[0];
    int max_E_pos = 0;
    for (int i=1; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        if (max_E < shannon_entropy[i])
        {
            max_E = shannon_entropy[i];
            max_E_pos = i;
        }
    }

    clock_t finish_t = clock();
    double sec = double(finish_t - start_t) / CLOCKS_PER_SEC;
    std::cout<<"Total time elapsed: "<<sec<<" seconds\n";
    eye_pos_to_bsphere_center = viewpoint_candidates[max_E_pos];
#undef VIEWPOINT_SAMPLES_NUM
}

void MyRender_MS::sample_using_shannon_entropy_II()
{
#define VIEWPOINT_SAMPLES_NUM 100 
    if (!mesh)
    {
        std::cout<<"plz load a mesh first....\n";
        return;
    }
    if (showing_image_type != model_space_curvature_image &&
        showing_image_type != radial_curvature_image &&
        showing_image_type != view_dependent_curvature_image)
    {
        std::cout<<"plz select a correct image type first\n";
        return;
    }

    //generating viewpoint candidates
	float r = mesh->bsphere.r;   //radius of bounding sphere

	//RandomDistributeViews(VIEWPOINT_SAMPLES_NUM, r);

	EqualDistributeViews(10, 10, r);

    //reset trackball quat, time accumulator
    sgi_trackball_space::trackball(trackball_quat,0,0,0,0);

    clock_t start_t = clock();
    float shannon_entropy[VIEWPOINT_SAMPLES_NUM];

    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        //set auto_sampling_veiwpoint and compute the entropy
        eye_pos_to_bsphere_center = viewpoint_candidates[i];
        paintGL();
        shannon_entropy[i] = compute_shannon_entropy_II();
    }

    // find max Shannon entropy 
    float max_E = shannon_entropy[0];
    int max_E_pos = 0;
    for (int i=1; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        if (max_E < shannon_entropy[i])
        {
            max_E = shannon_entropy[i];
            max_E_pos = i;
        }
    }

    clock_t finish_t = clock();
    double sec = double(finish_t - start_t) / CLOCKS_PER_SEC;
    std::cout<<"Total time elapsed: "<<sec<<" seconds\n";
    eye_pos_to_bsphere_center = viewpoint_candidates[max_E_pos];
	best_viewpoint_idx = max_E_pos;
#undef VIEWPOINT_SAMPLES_NUM
}

void MyRender_MS::sample_using_revised_entropy()
{
#define VIEWPOINT_SAMPLES_NUM 200
    if (!mesh)
    {
        std::cout<<"plz load a mesh first....\n";
        return;
    }
    if (showing_image_type != model_space_curvature_image &&
        showing_image_type != radial_curvature_image &&
        showing_image_type != view_dependent_curvature_image)
    {
        std::cout<<"plz select a correct image type first\n";
        return;
    }

    //generating viewpoint candidates
    viewpoint_candidates.clear();
    float r = mesh->bsphere.r;   //radius of bounding sphere
    time_t tm;
    srand((unsigned)time(&tm));  //set seed;
    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        float alpha = (rand()%360)/360.0f*2*PI;        //angle with positive-x-axis
        float beta = (rand()%360)/360.0f*2*PI - PI;  //angle with xy-plane
        point p;
        p[0] = r*cos(beta)*sin(alpha);
        p[1] = r*cos(beta)*cos(alpha);
        p[2] = r*sin(beta);
        viewpoint_candidates.push_back(p);
    }

    //reset trackball quat, time accumulator
    sgi_trackball_space::trackball(trackball_quat,0,0,0,0);

    clock_t start_t = clock();
    float shannon_entropy[VIEWPOINT_SAMPLES_NUM];
    float standard_deviations[VIEWPOINT_SAMPLES_NUM];

    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        //set auto_sampling_veiwpoint and compute the entropy
        eye_pos_to_bsphere_center = viewpoint_candidates[i];
        paintGL();
        shannon_entropy[i] = compute_shannon_entropy();
        standard_deviations[i] = current_standard_deviation; //computed in compute_shannon_entropy()
    }

    // get mean std_deviation
    float mean_standard_deviation=0;
    float revised_entropy[VIEWPOINT_SAMPLES_NUM];

    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        mean_standard_deviation += standard_deviations[i];
    }
    mean_standard_deviation/=VIEWPOINT_SAMPLES_NUM;
    // get rebalanced E
    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        standard_deviations[i] -= mean_standard_deviation;
        standard_deviations[i] = abs(standard_deviations[i]);
        revised_entropy[i] -= standard_deviations[i];
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
    double sec = double(finish_t - start_t) / CLOCKS_PER_SEC;
    std::cout<<"Total time elapsed: "<<sec<<" seconds\n";
    eye_pos_to_bsphere_center = viewpoint_candidates[max_E_pos];
#undef VIEWPOINT_SAMPLES_NUM
}

void MyRender_MS::sample_using_revised_entropy_II()
{
#define VIEWPOINT_SAMPLES_NUM 200
    if (!mesh)
    {
        std::cout<<"plz load a mesh first....\n";
        return;
    }
    if (showing_image_type != view_dependent_curvature_image)
    {
        std::cout<<"plz select a correct image type first\n";
        return;
    }

    //generating viewpoint candidates
    viewpoint_candidates.clear();
    float r = mesh->bsphere.r;   //radius of bounding sphere
    time_t tm;
    srand((unsigned)time(&tm));  //set seed;
    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        float alpha = (rand()%360)/360.0f*2*PI;        //angle with positive-x-axis
        float beta = (rand()%360)/360.0f*2*PI - PI;  //angle with xy-plane
        point p;
        p[0] = r*cos(beta)*sin(alpha);
        p[1] = r*cos(beta)*cos(alpha);
        p[2] = r*sin(beta);
        viewpoint_candidates.push_back(p);
    }

    //reset trackball quat, time accumulator
    sgi_trackball_space::trackball(trackball_quat,0,0,0,0);

    clock_t start_t = clock();
    float shannon_entropy[VIEWPOINT_SAMPLES_NUM];
    float standard_deviations[VIEWPOINT_SAMPLES_NUM];

    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        //set auto_sampling_veiwpoint and compute the entropy
        eye_pos_to_bsphere_center = viewpoint_candidates[i];
        paintGL();
        shannon_entropy[i] = compute_shannon_entropy_II();
        standard_deviations[i] = current_standard_deviation; //computed in compute_shannon_entropy_II()
    }

    // get mean std_deviation
    current_avg_stand_deviation=0;
    float revised_entropy[VIEWPOINT_SAMPLES_NUM];

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
    double sec = double(finish_t - start_t) / CLOCKS_PER_SEC;
    std::cout<<"Total time elapsed: "<<sec<<" seconds\n";
    eye_pos_to_bsphere_center = viewpoint_candidates[max_E_pos];

    // set has been sampled, and the data can be used to compute a certain view's revised entropy
    sampled_using_revised_entropy_II = true;

	// 将image类型设置为curvature [2/16/2012 Han]
	showing_image_type = original_image;

#undef VIEWPOINT_SAMPLES_NUM
}

void MyRender_MS::sample_using_mesh_saliency()
{
#define VIEWPOINT_SAMPLES_NUM 400 
    if (!mesh)
    {
        std::cout<<"plz load a mesh first....\n";
        return;
    }
    if (showing_image_type != mesh_saliency)
    {
        std::cout<<"plz select a correct image type first\n";
        return;
    }

    //generating viewpoint candidates
    viewpoint_candidates.clear();
    float r = mesh->bsphere.r;   //radius of bounding sphere
#if 0
    time_t tm;
    srand((unsigned)time(&tm));  //set seed;
    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        float alpha = (rand()%360)/360.0f*2*PI;        //angle with positive-x-axis
        float beta = (rand()%360)/360.0f*2*PI - PI;  //angle with xy-plane
        point p;
        p[0] = r*cos(beta)*sin(alpha);
        p[1] = r*cos(beta)*cos(alpha);
        p[2] = r*sin(beta);
        viewpoint_candidates.push_back(p);
    }
#else
    generate_viewpoint_candidates(VIEWPOINT_SAMPLES_NUM,viewpoint_candidates); // 来自OUR程序
    for (unsigned int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        viewpoint_candidates[i] *= r;
    }
#endif

    //reset trackball quat, time accumulator
    sgi_trackball_space::trackball(trackball_quat,0,0,0,0);

    clock_t start_t = clock();
    mesh->compute_mesh_saliency();
    float saliency[VIEWPOINT_SAMPLES_NUM];

    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        //set auto_sampling_veiwpoint and compute the entropy
        eye_pos_to_bsphere_center = viewpoint_candidates[i];
        paintGL();
        saliency[i] = compute_mesh_saliency();
    }

    // find max Shannon entropy 
    float max_E = saliency[0];
    int max_E_pos = 0;
    for (int i=1; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        if (max_E < saliency[i])
        {
            max_E = saliency[i];
            max_E_pos = i;
        }
    }

    clock_t finish_t = clock();
    secs = double(finish_t - start_t) / CLOCKS_PER_SEC;
    std::cout<<"Total time elapsed: "<<secs<<" seconds\n";
    eye_pos_to_bsphere_center = viewpoint_candidates[max_E_pos];
#undef VIEWPOINT_SAMPLES_NUM
}

void MyRender_MS::testing_a_round()
{
#define VIEWPOINT_SAMPLES_NUM 360
    if (!mesh)
    {
        std::cout<<"plz load a mesh first....\n";
        return;
    }

    //generating viewpoint candidates
    viewpoint_candidates.clear();
    float r = mesh->bsphere.r;   //radius of bounding sphere
    // sampling a round 
    for (int i=-90; i<VIEWPOINT_SAMPLES_NUM-90; i++)
    {
        float alpha = i/360.0f*2*PI;        //angle with positive-x-axis
        point p;
        p[0] = r*cos(alpha);
        p[1] = 0;
        p[2] = r*sin(alpha);
        viewpoint_candidates.push_back(p);
    }

    float entropy[VIEWPOINT_SAMPLES_NUM];
    float stddeviation[VIEWPOINT_SAMPLES_NUM];
    float entropy_r[VIEWPOINT_SAMPLES_NUM];

    std::ofstream of;
    of.open("around.txt",std::ios::out);
    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        //set auto_sampling_veiwpoint and compute the entropy
        eye_pos_to_bsphere_center = viewpoint_candidates[i];

        //showing_image_type = mesh_saliency;
        //paintGL();
        //saliency[i] = compute_mesh_saliency();

        showing_image_type = view_dependent_curvature_image;
        paintGL();
        entropy[i] = compute_shannon_entropy_II();
        stddeviation[i] = current_standard_deviation;
    }

    double esum = 0;
    double ssum = 0;
    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        esum += entropy[i];
        ssum += stddeviation[i];
    }

    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
    {
        entropy_r[i] = entropy[i] - 3*abs(stddeviation[i] - ssum/VIEWPOINT_SAMPLES_NUM)*esum/ssum;
    }


    for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
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
#undef VIEWPOINT_SAMPLES_NUM
}

void MyRender_MS::save_result()
{
    // save result to file:
    std::ofstream msresult("msresult",std::ios::out | std::ios::app);
    msresult<<secs<<eye_pos_to_bsphere_center<<std::endl;
    msresult.close();
}

int MyRender_MS::GetMeshTriNum()
{
	return mesh->faces.size();
}

point MyRender_MS::GetViewPoint()
{
	if (best_viewpoint_idx >= viewpoint_candidates_spherical.size())
		return point(0, 0, 0);
	return viewpoint_candidates_spherical[best_viewpoint_idx];
}

void MyRender_MS::RandomDistributeViews(int numViews, float r)
{
	viewpoint_candidates.clear();
	viewpoint_candidates_spherical.clear();
	time_t tm;
	srand((unsigned)time(&tm));  //set seed;

	for (int i=0; i<numViews; i++)
	{
	    float alpha = (rand()%360)/360.0f*2*PI;        //angle with positive-x-axis
	    float beta = (rand()%360)/360.0f*2*PI - PI;  //angle with xy-plane
	    point p;
	    p[0] = r*cos(beta)*sin(alpha);
	    p[1] = r*cos(beta)*cos(alpha);
	    p[2] = r*sin(beta);
	    viewpoint_candidates.push_back(p);

		p[0] = alpha;
		p[1] = beta;
		p[2] = r;
		viewpoint_candidates_spherical.push_back(p);		
	}
}
void MyRender_MS::EqualDistributeViews(int numViewsH, int numViewsV, float r)
{
	viewpoint_candidates.clear();
	viewpoint_candidates_spherical.clear();
	float alphaStep, betaStep;
	alphaStep = 2*PI / (numViewsH+1);
	betaStep = 2*PI / (numViewsV+1);
	for (int i=0; i<numViewsH; i++)
		for (int j = 0; j < numViewsV; j++)
		{
			float alpha = (i+1)*alphaStep;        //angle with positive-x-axis
			float beta = (j+1)*betaStep;  //angle with xy-plane
			point p;
			p[0] = r*cos(beta)*sin(alpha);
			p[1] = r*cos(beta)*cos(alpha);
			p[2] = r*sin(beta);
			viewpoint_candidates.push_back(p);

			p[0] = alpha;
			p[1] = beta;
			p[2] = r;
			viewpoint_candidates_spherical.push_back(p);
		}
}

void MyRender_MS::GetBestViewDist(float dis, Image_Type type)
{
	float r = mesh->bsphere.r;
	mesh->bsphere.r += dis*2;
	switch (type)
	{
	//case model_space_curvature_image:
	//	break;
	default:
		sample_using_shannon_entropy_II();
	}
	//mesh->bsphere.r = r;
}