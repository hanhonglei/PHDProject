// MeshSimpView.cpp : implementation of the CMeshSimpView class
//
#include "stdafx.h"
#include "MeshSimpDoc.h"
#include <gl/GL.h>
#include <gl/GLU.h>
#include "MeshSimp.h"

#include "MeshSimpView.h"
#include <math.h>
#include <gl/glut.h>




#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// arcball
// scene parameters
const vec_arcball centre( 0.0f, 0.0f, 0.0f );
const vec_arcball up( 0.0f, 1.0f, 0.0f );
float SPHERE_RADIUS = 1.1f;
const float PI = 3.141592654f;

float no_mat[] = {0.0f, 0.0f, 0.0f, 0.5f};
float mat_red[] = {0.7f, 0.1f, 0.1f, 0.5};
float mat_green[] = {0.1f, 0.7f, 0.1f, 0.5};
float mat_ambient[] = {0.7f, 0.7f, 0.7f, 0.5};
float mat_diffuse[] = {0.1f, 0.5f, 0.8f, 0.5};
float mat_specular[] = {0.1f, 0.1f, 0.1f, 0.5};
float mat_diffuse_sel[] = {0.3f, 0.2f, 0.2f, 0.5};
float mat_seg[12][4] =  {
	{0.2f, 0.0f, 0.0f, 0.5},
	{0.0f, 0.2f, 0.0f, 0.5},
	{0.0f, 0.0f, 0.2f, 0.5},
	{0.4f, 0.0f, 0.0f, 0.5},
	{0.0f, 0.4f, 0.0f, 0.5},
	{0.0f, 0.0f, 0.6f, 0.5},
	{0.8f, 0.0f, 0.0f, 0.5},
	{0.0f, 0.8f, 0.0f, 0.5},
	{0.0f, 0.0f, 0.8f, 0.5},
	{1.0f, 0.0f, 0.0f, 0.5},
	{0.0f, 0.1f, 0.0f, 0.5},
	{0.0f, 0.0f, 0.1f, 0.5}
};

float		m_lightPos[4] = {10.0f, 1.0f, 39.0f, 1.0f};			//position of the light source in the world						 


// CMeshSimpView

IMPLEMENT_DYNCREATE(CMeshSimpView, CView)

BEGIN_MESSAGE_MAP(CMeshSimpView, CView)
	//{{AFX_MSG_MAP(CTestScrollOpenglView)
	ON_WM_CREATE()
	ON_WM_SIZE()
	ON_WM_DESTROY()
	ON_WM_ERASEBKGND()
	//}}AFX_MSG_MAP

	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CView::OnFilePrintPreview)
	ON_WM_MOUSEMOVE()
//	ON_COMMAND(ID_FILL_WIREFRAME, &CMeshSimpView::OnFillWireframe)
//ON_WM_MENUSELECT()
ON_UPDATE_COMMAND_UI(ID_SHADING_SMOOTH, &CMeshSimpView::OnUpdateShadingSmooth)
ON_COMMAND(ID_SHADING_SMOOTH, &CMeshSimpView::OnShadingSmooth)
ON_COMMAND(ID_FILL_WIREFRAME, &CMeshSimpView::OnFillWireframe)
ON_UPDATE_COMMAND_UI(ID_FILL_WIREFRAME, &CMeshSimpView::OnUpdateFillWireframe)
//ON_WM_KEYDOWN()
ON_WM_KEYUP()
ON_WM_LBUTTONUP()
ON_WM_MOUSEWHEEL()
ON_WM_KEYDOWN()
ON_WM_LBUTTONDOWN()
ON_COMMAND(ID_EDIT_UNSELECTALL, &CMeshSimpView::OnEditUnselectall)
ON_UPDATE_COMMAND_UI(ID_EDIT_UNSELECTALL, &CMeshSimpView::OnUpdateEditUnselectall)
ON_WM_TIMER()
ON_UPDATE_COMMAND_UI(ID_MY_STATE, &CMeshSimpView::OnUpdateStatus)
ON_UPDATE_COMMAND_UI(ID_MY_STATE2, &CMeshSimpView::OnUpdateStatus2)
ON_UPDATE_COMMAND_UI(ID_MY_STATE3, &CMeshSimpView::OnUpdateStatus3)
ON_COMMAND(ID_NEXT_COLL, &CMeshSimpView::OnNextColl)
ON_COMMAND(ID_PREV_COLL, &CMeshSimpView::OnPrevColl)
ON_COMMAND(ID_SHADING_DRAWCOLLDIST, &CMeshSimpView::OnShadingDrawcolldist)
ON_UPDATE_COMMAND_UI(ID_SHADING_DRAWCOLLDIST, &CMeshSimpView::OnUpdateShadingDrawcolldist)
ON_COMMAND(ID_SHADING_LIGHTING, &CMeshSimpView::OnShadingLighting)
ON_UPDATE_COMMAND_UI(ID_SHADING_LIGHTING, &CMeshSimpView::OnUpdateShadingLighting)
ON_COMMAND(ID_SHADING_SKELETON, &CMeshSimpView::OnShadingSkeleton)
ON_UPDATE_COMMAND_UI(ID_SHADING_SKELETON, &CMeshSimpView::OnUpdateShadingSkeleton)
ON_COMMAND(ID_SHADING_TRANSPARANT, &CMeshSimpView::OnShadingTransparant)
ON_UPDATE_COMMAND_UI(ID_SHADING_TRANSPARANT, &CMeshSimpView::OnUpdateShadingTransparant)
ON_COMMAND(ID_SHADING_DRAWSKELMAP, &CMeshSimpView::OnShadingDrawskelmap)
ON_UPDATE_COMMAND_UI(ID_SHADING_DRAWSKELMAP, &CMeshSimpView::OnUpdateShadingDrawskelmap)
ON_UPDATE_COMMAND_UI(ID_NEXT_COLL, &CMeshSimpView::OnUpdateNextColl)
ON_UPDATE_COMMAND_UI(ID_PREV_COLL, &CMeshSimpView::OnUpdatePrevColl)
ON_COMMAND(ID_SHADING_DRAWSEGMENT, &CMeshSimpView::OnShadingDrawsegment)
ON_UPDATE_COMMAND_UI(ID_SHADING_DRAWSEGMENT, &CMeshSimpView::OnUpdateShadingDrawsegment)
ON_COMMAND(ID_SHADING_DRAWMESH, &CMeshSimpView::OnShadingDrawmesh)
ON_UPDATE_COMMAND_UI(ID_SHADING_DRAWMESH, &CMeshSimpView::OnUpdateShadingDrawmesh)
ON_COMMAND(ID_VIEW_SELECTION, &CMeshSimpView::OnViewSelection)
ON_UPDATE_COMMAND_UI(ID_VIEW_SELECTION, &CMeshSimpView::OnUpdateViewSelection)
ON_COMMAND(ID_SHADING_RESETVIEW, &CMeshSimpView::OnShadingResetview)
ON_COMMAND(ID_SHADING_DRAW_VIEW_SPHERE, &CMeshSimpView::OnShadingDrawViewSphere)
ON_COMMAND(ID_SHADING_ALL6VIEWS, &CMeshSimpView::OnShadingAll6views)
ON_COMMAND(ID_VIEWPOINTSELECTION_UPRIGHTDIR, &CMeshSimpView::OnViewpointselectionUprightdir)
ON_COMMAND(ID_SHADING_DRAWNORMAL, &CMeshSimpView::OnShadingDrawnormal)
END_MESSAGE_MAP()



CMeshSimpView::CMeshSimpView()
{
	// TODO: add construction code here
	// leo : 初始化
	m_pDC = NULL;
	//resetOrientation();
	m_bArcBallStart = false;
	m_mouseMode = ROTATE;
	m_bSelAny = false;
	bFill_ = true;
	bSmooth_ = false;

	bLighting = true;
	bDrawSkel = true;
	bDrawMesh = true;
	bTrans = false;
	drawMeshType = PMesh::NONE;
	activeTriNum = 0;

	m_fNearClip = 10.f/*0.1f*/;
	m_fFarClip = 4000.0f;
	m_fov = 45;

	eye.x=0.0f, eye.y=0.0f, eye.z=8.0f;
	showing_type = SHOWING_MODEL;

	bShowUprightDir = false;
}

CMeshSimpView::~CMeshSimpView()
{
	// 释放指针
	if (m_pDC != NULL)
	{
		delete m_pDC;
	}
}

void CMeshSimpView::InitNewDoc()
{
	m_bSelAny = false;
	m_selIndexList.clear();
}

/////////////////////////////////////////////////////////////////////////////
// CMeshSimpView diagnostics
//#ifdef _DEBUG

CMeshSimpDoc* CMeshSimpView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CMeshSimpDoc)));
	return (CMeshSimpDoc*)m_pDocument;
}
//#endif //_DEBUG

void CMeshSimpView::DisplayModel()
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (SHOWING_NORMALMAP == showing_type)
		glClearColor(0.f, 0.f, 0.f, 0.5f);
	else
		glClearColor(1.f, 1.f, 1.f, 0.5f);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	switch (pDoc->m_mType)
	{
	case  CMeshSimpDoc::CurrentModelType::VIEW_SELECT:
		if (pDoc->m_viewMyRender)
			pDoc->m_viewMyRender->paintGL();
		break;
	case  CMeshSimpDoc::CurrentModelType::VIEW_SELECT_OUR:
		if(pDoc->m_viewMyRenderOur)
			pDoc->m_viewMyRenderOur->paintGL();
		break;
	default:

		glColor3f(0.2,0.2,0.2);
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
		glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 0);
		glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, no_mat);

		if (pDoc->m_pProgMesh)
		{		
			pDoc->m_pProgMesh->DrawOriginalMesh(drawMeshType, bSmooth_);

			if (pDoc->m_pProgMesh->isCollapsing || !pDoc->m_pProgMesh->bSimplified || !bDrawSkel)
				return;
			pDoc->m_pProgMesh->DrawSimplifiedVertices();	// 绘制骨架

		}
		else if (pDoc->m_pQslimMesh)
		{
			switch (showing_type)
			{
				/************************************************************************/
				/*    Mesh Representation                                               */
				/************************************************************************/
			case SHOWING_VIEW_SPHERE_MAP:

			case SHOWING_MODEL:
			case SHOWING_SEGMENTATION:
				glEnable(GL_LIGHTING);

				// 使用phong模型进行绘制 [5/25/2012 Han]
				m_phongShader.setUniform3f("lightPos", &VECTOR3D(m_lightPos[0], m_lightPos[1]/*+50.0f*/, m_lightPos[2]));
				m_phongShader.setUniform3f("camPos", &VECTOR3D(eye.x, eye.y, eye.z));
				m_phongShader.begin();

				// 取消背向面剔除有利于绘制出非体积型面 [7/26/2012 Han]
				glDisable(GL_CULL_FACE);
				mx_render_model(*pDoc->m_pQslimMesh);

				m_phongShader.end();

				if (showing_type == SHOWING_VIEW_SPHERE_MAP)
				{
					/************************************************************************/
					/*  Viewing sphere map to color                                         */
					/************************************************************************/
					pDoc->DrawViewSphere(pDoc->GetModelRadius()*3);
				}

				break;
			case SHOWING_MESH_SALIENCY:
			case VISIBLE_TEST:
			case SHOWING_VIEW_DEPENDENT_CURVATURE:
				//DrawVDCurvature();
				glDisable(GL_LIGHTING);
				if (bSmooth_)
					mx_render_model(*pDoc->m_pQslimMesh);
				else
				{
					glBegin(GL_TRIANGLES);
					for (int i = 0; i < pDoc->m_pQslimMesh->face_count(); i++)
					{
						if (!pDoc->m_pQslimMesh->face_is_valid(i))
							continue;
						unsigned int idx = pDoc->m_pQslimMesh->face(i).colorIndex;
						idx = idx < 0?0:idx;
						idx = idx > 255?255:idx;
						glColor3dv(&pDoc->m_pQslimMesh->LUT_CourbureClust[idx*3]);
						for (int vi = 0; vi < 3; vi++)							
							glVertex3fv(pDoc->m_pQslimMesh->vertex(pDoc->m_pQslimMesh->face(i).v[vi]).as.pos);
					}
					glEnd();
				}
				break;
			case SHOWING_CURVATURE_QUANTITY:
				glDisable(GL_LIGHTING);
				glBegin(GL_TRIANGLES);
				for (int i = 0; i < pDoc->m_pQslimMesh->face_count(); i++)
				{
					if (!pDoc->m_pQslimMesh->face_is_valid(i))
						continue;
					for (int vi = 0; vi < 3; vi++)	
					{
						unsigned int idx = pDoc->m_pQslimMesh->vertex(pDoc->m_pQslimMesh->face(i).v[vi]).colorIndex;
						idx = idx < 0?0:idx;
						idx = idx > 255?255:idx;
						glColor3ub(idx, idx, idx);
						glVertex3fv(pDoc->m_pQslimMesh->vertex(pDoc->m_pQslimMesh->face(i).v[vi]).as.pos);
					}
				}
				glEnd();
				break;
			case 		SHOWING_NORMALMAP:
				//DrawVDCurvature();
				glDisable(GL_LIGHTING);
				glBegin(GL_TRIANGLES);
					for (int i = 0; i < pDoc->m_pQslimMesh->face_count(); i++)
					{
						if (!pDoc->m_pQslimMesh->face_is_valid(i))
							continue;
						if( pDoc->m_pQslimMesh->normal_binding()==MX_PERFACE )  
							glColor3sv(pDoc->m_pQslimMesh->normal(i).raw());
						for (int vi = 0; vi < 3; vi++)		
						{
							if (pDoc->m_pQslimMesh->normal_binding()== MX_PERVERTEX)
								glColor3sv(pDoc->m_pQslimMesh->normal(pDoc->m_pQslimMesh->face(i).v[vi]).raw());
							glVertex3fv(pDoc->m_pQslimMesh->vertex(pDoc->m_pQslimMesh->face(i).v[vi]).as.pos);
						}

					}
				glEnd();

				break;

			case SHOWING_MODEL_WIRE:
				glEnable(GL_POLYGON_OFFSET_FILL);
				glPolygonOffset(1.0, 1.0);
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				//glShadeModel(GL_FLAT);
				mx_render_model(*pDoc->m_pQslimMesh);

				//glBegin(GL_TRIANGLES);
				//for (int i = 0; i < pDoc->m_pQslimMesh->face_count(); i++)
				//{
				//	glColor3dv(pDoc->m_pQslimMesh->v_infor_colors[i]);
				//	for (int vi = 0; vi < 3; vi++)							
				//		glVertex3fv(pDoc->m_pQslimMesh->vertex(pDoc->m_pQslimMesh->face(i).v[vi]).as.pos);
				//}
				//glEnd();
				glDisable(GL_POLYGON_OFFSET_FILL);
				if (!bFill_)
				{
					glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

					glColor3f(1 , 1 , 1);
					glBegin(GL_TRIANGLES);
					for (int i = 0; i < pDoc->m_pQslimMesh->face_count(); i++)
					{
						if (!pDoc->m_pQslimMesh->face_is_valid(i))
							continue;

						for (int vi = 0; vi < 3; vi++)							
							glVertex3fv(pDoc->m_pQslimMesh->vertex(pDoc->m_pQslimMesh->face(i).v[vi]).as.pos);
					}
					glEnd();
				}
				break;

			case SHOWING_MODEL_VOXELS:
			case SHOWING_WHOLE_SKELETON:
			case SHOWING_PRIMARY_SKELETON: 
			case SHOWING_FINAL_SKELETON:
			// Draw bounding sphere [5/25/2012 Han]
			glEnable(GL_CULL_FACE); // backface culling
			glDisable(GL_LIGHTING);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

			glColor4f(0.2f, 0.6, 0.1, 0.3);

			glutSolidSphere(pDoc->GetModelRadius(), 20, 20);				
				break;
			default:
				break;
			}
			// 如果绘制向上方向的话，则进行绘制 [9/18/2012 Han]
			if (bShowUprightDir)
				pDoc->DrawUprightDir();
		}
		break;
	}
}
/////////////////////////////////////////////////////////////////////////////
// CMeshSimpView drawing

void CMeshSimpView::OnDraw(CDC* pDC)
{	
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	// 每次都重新确定视口 [5/7/2012 Han]
	glViewport(0, 0, width_, height_);

	// 如果载入图像成功，则进行绘制
	if (pDoc->m_isLoaded)
	{
		if (bSmooth_)
			glShadeModel(GL_SMOOTH); // already defined in initOpenGL
		else 
			glShadeModel(GL_FLAT);

		if (bFill_) // fill in triangles or just display outlines?
			glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
		else
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );


		if (bLighting)
			glEnable(GL_LIGHTING);
		else
			glDisable(GL_LIGHTING);

		glLoadIdentity();
		arcball_rotate();

		DisplayModel();

		// 交换缓存到设备显示
		SwapBuffers(wglGetCurrentDC());
	}
}

/////////////////////////////////////////////////////////////////////////////
// CMeshSimpView message handlers

int CMeshSimpView::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
	if (CView::OnCreate(lpCreateStruct) == -1)
		return -1;
	// leo : 初始化opengl
	IniOpenGL();

	// 暂时不进行计时 [7/24/2012 Han]
	//SetTimer(1, 20, NULL);
	return 0;
}

void CMeshSimpView::OnSize(UINT nType, int cx, int cy) 
{
	CView::OnSize(nType, cx, cy);
	// 调整视窗	
	if (cy > 0)
	{
		int nWinWidth  = cx;
		int nWinHeight = cy;
		if(nWinHeight == 0)
		{
			nWinHeight = 1;
		}

		width_ = nWinWidth;
		height_ = nWinHeight;

		AdjustView( nWinWidth,nWinHeight);
	}

}

void CMeshSimpView::OnDestroy() 
{
	CView::OnDestroy();

	// TODO: Add your message handler code here
	// leo:应该添加这些释放资源的操作
	HGLRC   hrc;

	hrc = ::wglGetCurrentContext();

	::wglMakeCurrent(NULL,  NULL);

	if (hrc)
		::wglDeleteContext(hrc);	

}
// leo:不进行背景清除，以免绘制两次，产生闪烁
BOOL CMeshSimpView::OnEraseBkgnd(CDC* pDC) 
{
	// TODO: Add your message handler code here and/or call default
	return TRUE;
}
// leo : 初始化openggl绘图环境
void CMeshSimpView::IniOpenGL()
{
	m_pDC = new CClientDC(this);
	PIXELFORMATDESCRIPTOR pfd;
	int n;
	//定义一个绘制上下文的句柄
	HGLRC   hrc;

	//初始化过程中主要就是初始化了一个客户区的设备环境指针
	ASSERT(m_pDC != NULL);

	//建立应用所需的像素格式，并与当前设备上下文相关连
	if (!bSetPixelFormat())
		return;

	//得到指定设备环境的象素模式索引
	n = ::GetPixelFormat(m_pDC->GetSafeHdc());

	//根据上面得到的索引值来声明一个象素模式
	::DescribePixelFormat(m_pDC->GetSafeHdc(), n, sizeof(pfd), &pfd);

	//创建一个上下文设备环境
	hrc = wglCreateContext(m_pDC->GetSafeHdc());

	//将刚生成的设备上下文指针设为当前环境
	wglMakeCurrent(m_pDC->GetSafeHdc(), hrc);		

	glShadeModel(GL_SMOOTH); // Gouraud shading
	glClearColor(1.f, 1.f, 1.f, 0.5f);
	glClearDepth(1.0f);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE); // backface culling
	glDisable(GL_TEXTURE_2D);
	// leo:使用后备缓存绘制，避免闪烁
	glDrawBuffer (GL_BACK);


	//light properties
	float ambient[] = {0.0f, 0.0f, 0.0f, 1.0f};
	float diffuse[] = {0.5f, 0.5f, 0.5f, 1.0f};
	float specular[] = {0.3f, 0.3f, 0.3f, 1.0f};

	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
	glLightfv(GL_LIGHT0, GL_POSITION, m_lightPos);

	//light model properties
	float model_ambient[] = {0.4f, 0.4f, 0.4f, 1.0f};
	int model_two_side = 1;                                //0=2sided, 1=1sided
	int viewpoint = 0;                                     //0=infiniteViewpoint, 1=localViewpoint
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, model_ambient);     //small white ambient light

	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, viewpoint);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

	//init shaders
	m_phongShader.init("Shaders/phong.cg", "Shaders/phong.cg");
}
// leo ：设置适于OpenGL使用的像素格式
BOOL CMeshSimpView::bSetPixelFormat()
{
	//定义一种像素格式
	static PIXELFORMATDESCRIPTOR pfd =
	{   
		sizeof(PIXELFORMATDESCRIPTOR),			// size of this pfd
		1,									// version number
		PFD_DRAW_TO_WINDOW |				// support window
		PFD_SUPPORT_OPENGL |				// support OpenGL  支持OpenGL
		PFD_DOUBLEBUFFER,					// double buffered 支持又缓冲
		PFD_TYPE_RGBA,						// RGBA type使用RGBA模式，不用调色板
		24,									// 24-bit color depth  使用24位真彩色
		0, 0, 0, 0, 0, 0,					// color bits ignored
		0,									// no alpha buffer
		0,									// shift bit ignored
		0,									// no accumulation buffer
		0, 0, 0, 0,							// accum bits ignored
		32,									// 32-bit z-buffer   32位Z轴缓冲
		0,									// no stencil buffer
		0,									// no auxiliary buffer
		PFD_MAIN_PLANE,						// main layer
		0,									// reserved
		0, 0, 0								// layer masks ignored
	};
	int pixelformat;

	//如果可以得到指定的像素格式
	if ( (pixelformat = ChoosePixelFormat(m_pDC->GetSafeHdc(), &pfd)) == FALSE )
	{
		AfxMessageBox(LPCSTR("ChoosePixelFormat failed"));
		return false;
	}

	//用上面取到的格式设置设备环境
	if (SetPixelFormat(m_pDC->GetSafeHdc(), pixelformat, &pfd) == FALSE)
	{
		AfxMessageBox(LPCSTR("SetPixelFormat failed"));
		return false;
	}
	return true;
}

// reset the mesh orientation 将物体回复到最初状态，消除摄像机旋转、缩放的影响
void CMeshSimpView::resetOrientation() 
{	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	eye.x = 0.0f;
	eye.y = 0.0f;
	eye.z = 8.0f;
	// calc model's distance to the viewer, when it's projection occupy full screen [3/28/2012 Han]
	// l = R * arcsin(α/2)
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	float r = pDoc->GetModelRadius();
	float rd = r;
	if ( r > 0)
	{
		// 屏幕空间和三维空间的对应关系:
		// glFrstrom中,假定eye位置为原地,最近裁剪面表示投影平面里原点的距离,投影四边形的4个顶点位置就确定了投影平面的大小,以及Fov角度
		// gluLookAt函数所指定的eye位置即调整后的视点到原点的距离,其实相当于做相反方向的glTranslate函数
		// 故而,可以利用以下公式得到当模型包围球半径为r时,模型离视点的距离为多少时,模型投影能占据整个屏幕
		r /= sin(M_PI*(m_fov/2.0)/180);
		eye.z = r;
	}
	// nearclip 修改为可能的最远距离 [7/9/2012 Han]
	m_fNearClip = (eye.z - rd)/2.0f;
	// far clip 修改为可能的最近距离 [7/24/2012 Han]
	m_fFarClip = eye.z + 7 * rd;

	arcball_reset();
	// 取消模型向上方向，避免其干扰 [3/13/2014 Han]
	if (pDoc->m_baseViewpoint != NULL)
	{
		delete pDoc->m_baseViewpoint;
		pDoc->m_baseViewpoint = NULL;
	}
	
	AdjustView(width_, height_);	
}

void CMeshSimpView::SelectTri(CPoint point, float radius)
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	// 暂时只实现了progMesh的点选功能，所以，如果progMesh没有被载入，则不点选 [8/15/2013 Han]
	if (pDoc->m_pProgMesh == NULL)
		return;

	GLuint	buffer[512];										// Set Up A Selection Buffer
	GLint	hits;												// The Number Of Objects That We Selected

	// The Size Of The Viewport. [0] Is <x>, [1] Is <y>, [2] Is <length>, [3] Is <width>
	GLint	viewport[4];

	// This Sets The Array <viewport> To The Size And Location Of The Screen Relative To The Window
	glGetIntegerv(GL_VIEWPORT, viewport);
	glSelectBuffer(512, buffer);								// Tell OpenGL To Use Our Array For Selection

	// Puts OpenGL In Selection Mode. Nothing Will Be Drawn.  Object ID's and Extents Are Stored In The Buffer.
	(void) glRenderMode(GL_SELECT);

	glInitNames();												// Initializes The Name Stack
	glPushName(0);												// Push 0 (At Least One Entry) Onto The Stack

	glMatrixMode(GL_PROJECTION);								// Selects The Projection Matrix
	glPushMatrix();												// Push The Projection Matrix
	glLoadIdentity();											// Resets The Matrix

	// This Creates A Matrix That Will Zoom Up To A Small Portion Of The Screen, Where The Mouse Is.
	gluPickMatrix((GLdouble) point.x, (GLdouble) (viewport[3]-point.y), radius, radius, viewport);

	// 修改为依靠传入的参数决定透视坐标系
	gluPerspective(/*vertical field of view*/ m_fov,
		/*aspect ratio*/ /*(double) viewport.width/viewport.height,*/(double)width_/height_,
		/*znear*/ m_fNearClip, /*zfar*/ m_fFarClip);

	gluLookAt(
		eye.x, eye.y, eye.z,
		centre.x, centre.y, centre.z,
		up.x, up.y, up.z );
	//// set up the arcball using the current projection matrix
	//arcball_setzoom( SPHERE_RADIUS, eye, up );


	glMatrixMode(GL_MODELVIEW);									// Select The Modelview Matrix
	OnDraw(NULL);
	//pDoc->DisplayMesh();												// Render The Targets To The Selection Buffer
	glMatrixMode(GL_PROJECTION);								// Select The Projection Matrix
	glPopMatrix();												// Pop The Projection Matrix
	glMatrixMode(GL_MODELVIEW);									// Select The Modelview Matrix
	hits=glRenderMode(GL_RENDER);								// Switch To Render Mode, Find Out How Many
	// Objects Were Drawn Where The Mouse Was
	if (hits > 0)												// If There Were More Than 0 Hits
	{
		CMeshSimpDoc* pDoc = GetDocument();
		ASSERT_VALID(pDoc);

		int	choose = buffer[3];									// Make Our Selection The First Object
		int depth = buffer[1];									// Store How Far Away It Is 

		for (int loop = 1; loop < hits; loop++)					// Loop Through All The Detected Hits
		{
			// If This Object Is Closer To Us Than The One We Have Selected
			if (buffer[loop*4+1] < GLuint(depth))
			{
				choose = buffer[loop*4+3];						// Select The Closer Object
				depth = buffer[loop*4+1];						// Store How Far Away It Is
			}       
		}
		triangle &t = pDoc->m_pProgMesh->getTri(choose);
		t.setSel(true);
		m_selIndexList.remove(choose);
		m_selIndexList.push_back(choose);
		m_bSelAny = true;
		Invalidate();
	}
}


void CMeshSimpView::OnMouseMove(UINT nFlags, CPoint point)
{
	CMeshSimpDoc* pDoc = GetDocument();
	if (nFlags&MK_LBUTTON )
	{
		switch(m_mouseMode)
		{
		case NORMAL:
			break;
		case ROTATE:
			arcball_move(point.x,(height_ - point.y) - 1);

			// 使用view的mouse控制方式 [2/16/2012 Han]
			ASSERT_VALID(pDoc);
			if (pDoc->m_mType == CMeshSimpDoc::CurrentModelType::VIEW_SELECT)
			{
				pDoc->m_viewMyRender->add_trackball_quat(
					(2.0*m_mousePos.x - width_) / width_,
					(height_ - 2.0*m_mousePos.y) / height_,
					(2.0*point.x - width_)    / width_,
					(height_ - 2.0*point.y)    / height_);
				m_mousePos = point;
			}
			else if (pDoc->m_mType == CMeshSimpDoc::CurrentModelType::VIEW_SELECT_OUR)
			{
				pDoc->m_viewMyRenderOur->add_trackball_quat(
					-(2.0*m_mousePos.x - width_) / width_,
					-(height_ - 2.0*m_mousePos.y) / height_,
					-(2.0*point.x - width_)    / width_,
					-(height_ - 2.0*point.y)    / height_);
				m_mousePos = point;

			}
			if (pDoc->m_pQslimMesh != NULL)
				pDoc->m_pQslimMesh->bRotated = true;
			if (pDoc->m_pOriginalQMesh != NULL)
				pDoc->m_pOriginalQMesh->bRotated = true;
			break;
		case SELECT:
			SelectTri(point, 2.0);
		    break;
		default:
		    break;
		}
		Invalidate();
	}

	CView::OnMouseMove(nFlags, point);
}


void CMeshSimpView::OnLButtonDown(UINT nFlags, CPoint point)
{
	m_mousePos = point;
	switch(m_mouseMode)
	{
	case SELECT:
		SelectTri(point, 2.0);
		Invalidate();
		break;
	case  ROTATE:
			arcball_start(point.x,(height_ - point.y) - 1);
			m_bArcBallStart = true;
			break;

	default:
		break;
	}
	CView::OnLButtonDown(nFlags, point);
}

void CMeshSimpView::OnLButtonUp(UINT nFlags, CPoint point)
{
	m_bArcBallStart = false;
	CView::OnLButtonUp(nFlags, point);
}

void CMeshSimpView::OnUpdateShadingSmooth(CCmdUI *pCmdUI)
{
	pCmdUI->Enable();
	pCmdUI->SetCheck(bSmooth_);

}

void CMeshSimpView::OnShadingSmooth()
{
	bSmooth_ = !bSmooth_;
	Invalidate();
}

void CMeshSimpView::OnUpdateFillWireframe(CCmdUI *pCmdUI)
{
	pCmdUI->Enable();
	pCmdUI->SetCheck(!bFill_);
}
void CMeshSimpView::OnFillWireframe()
{
	bFill_ = !bFill_;
	Invalidate();
}

void CMeshSimpView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags)
{
	switch(nChar)
	{
	case VK_CONTROL:
		m_mouseMode = SELECT;
		// 转换鼠标指针
		::SetCursor(AfxGetApp()->LoadCursor(MAKEINTRESOURCE(IDC_CURSOR_SEL)));
		SetClassLong(m_hWnd,
			GCL_HCURSOR,
			(LONG)LoadCursor(AfxGetInstanceHandle(),
			MAKEINTRESOURCE(IDC_CURSOR_SEL)));
		Invalidate();
		break;
	}

	CView::OnKeyDown(nChar, nRepCnt, nFlags);
}

void CMeshSimpView::OnKeyUp(UINT nChar, UINT nRepCnt, UINT nFlags)
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	switch(nChar)
	{
	case VK_PRIOR:
		pDoc->OnResampleIncreaseby5();
		break;
	case VK_NEXT:
		pDoc->OnResampleReduceby5();
		break;
	case VK_UP:
		pDoc->OnResampleAdd1tri();
		break;
	case VK_DOWN:
		pDoc->OnResampleRemove1tri();
		break;
	case VK_CONTROL:
		m_mouseMode = ROTATE;
		// 转换鼠标指针
		SetClassLong( m_hWnd,
			GCL_HCURSOR,
			(LONG)AfxGetApp()->LoadStandardCursor(IDC_ARROW) );
		::SetCursor(AfxGetApp()->LoadStandardCursor(IDC_ARROW));
		break;
	case 'I' :
		{
			if (pDoc->m_mType == CMeshSimpDoc::CurrentModelType::VIEW_SELECT)
			{
				switch (pDoc->m_viewMyRender->get_showing_image_type())
				{
				case MyRender_MS::Image_Type::original_image:
					pDoc->m_viewMyRender->set_showing_image_type(MyRender_MS::Image_Type::depth_buffer_image);
					break;
				case MyRender_MS::Image_Type::depth_buffer_image:
					pDoc->m_viewMyRender->set_showing_image_type(MyRender_MS::Image_Type::mesh_saliency);
					break;
				case MyRender_MS::Image_Type::mesh_saliency:
					pDoc->m_viewMyRender->set_showing_image_type(MyRender_MS::Image_Type::model_space_curvature_image);
					break;
				case MyRender_MS::Image_Type::model_space_curvature_image:
					pDoc->m_viewMyRender->set_showing_image_type(MyRender_MS::Image_Type::radial_curvature_image);
					break;
				case MyRender_MS::Image_Type::radial_curvature_image:
					pDoc->m_viewMyRender->set_showing_image_type(MyRender_MS::Image_Type::view_dependent_curvature_image);
					break;
				case MyRender_MS::Image_Type::view_dependent_curvature_image:
					pDoc->m_viewMyRender->set_showing_image_type(MyRender_MS::Image_Type::original_image);
					break;
				}			
			}
			else if (pDoc->m_mType == CMeshSimpDoc::CurrentModelType::VIEW_SELECT_OUR)
			{
				int tmp = pDoc->m_viewMyRenderOur->showing_type;
				tmp++;
				tmp = tmp % 9;
				pDoc->m_viewMyRenderOur->showing_type = (MyRender_OUR::Image_Type)tmp;
			}
		}
		Invalidate();
	}

	CView::OnKeyUp(nChar, nRepCnt, nFlags);
}

void CMeshSimpView::SaveCurrentModelPic(wchar_t fn[128])
{
	DisplayModel();
	glFinish();

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

	glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, ppmWriter->data);
	wchar_t fnRe[256];
	swprintf( fnRe,   L"OutPutFiles\\%s", fn);

	ppmWriter->save(fnRe);
	delete ppmWriter;
}

// 每次的delta值为120
BOOL CMeshSimpView::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt)
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	float fzDelta = zDelta/1200.0f;
	//eye.x *= (1+fzDelta);
	//eye.y *= (1+fzDelta);
	//eye.z *= (1+fzDelta);
	eye = eye*(1+fzDelta);
	AdjustView(width_, height_);
	// Adjust LOD of the model depend on distance [9/28/2011 Han Honglei]
	pDoc->AdjustLOD(sqrt(pow(eye.x-centre.x, 2)+pow(eye.y-centre.y,2)+ pow(eye.z-centre.z,2)));

	if (pDoc->m_mType == CMeshSimpDoc::CurrentModelType::VIEW_SELECT_OUR)
	{
		pDoc->m_viewMyRenderOur->eye_pos_to_bsphere_center *= (1+fzDelta);

	}
	Invalidate();
	return CView::OnMouseWheel(nFlags, zDelta, pt);
}

float CMeshSimpView::GetViewDistance()
{
	return sqrt(pow(eye.x-centre.x, 2)+pow(eye.y-centre.y,2)+ pow(eye.z-centre.z,2));
}

void perspectiveGL( GLdouble fovY, GLdouble aspect, GLdouble zNear, GLdouble zFar )
{
	GLdouble fW, fH;
	fH = tan( fovY / 360 * PI ) * zNear;
	fW = fH * aspect;
	glFrustum( -fW, fW, -fH, fH, zNear, zFar );
}
void CMeshSimpView::AdjustView(int nWidth, int nHeight)
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	glViewport ( 0, 0, nWidth, nHeight );

	/* setup perspective m_camera with OpenGL */
	glMatrixMode(GL_PROJECTION);
	glClearColor(1.f, 1.f, 1.f, 0.5f);
	glLoadIdentity();
	// 修改为依靠传入的参数决定透视坐标系
	//gluPerspective(/*vertical field of view*/ m_fov.,
	//	/*aspect ratio*/ /*(double) viewport.width/viewport.height,*/(double)nWidth/nHeight,
	//	/*znear*/ 0.1, /*zfar*/ 1000.0);
	//glFrustum(-(double)nWidth/nHeight, (double)nWidth/nHeight,-1, 1, m_fNearClip, m_fFarClip);
	//glOrtho(-nWidth/2.0,nWidth/2.0,-nHeight/2.0,nHeight/2.0, m_fNearClip, m_fFarClip) ;
	perspectiveGL(m_fov, (double)nWidth/nHeight, m_fNearClip, m_fFarClip);

	if (pDoc->m_baseViewpoint != NULL)
	{
		gluLookAt(
		eye.x, eye.y, eye.z,
			centre.x, centre.y, centre.z,
			-pDoc->m_baseViewpoint->pos[0], -pDoc->m_baseViewpoint->pos[1],-pDoc->m_baseViewpoint->pos[2]);
	}
	else
	{
		if (eye.x == 0 && eye.z == 0)
		{
			gluLookAt(
				eye.x, eye.y, eye.z,
				centre.x, centre.y, centre.z,
				0, 0, -1 );
		}
		else
			gluLookAt(
			eye.x, eye.y, eye.z,
			centre.x, centre.y, centre.z,
			up.x, up.y, up.z );
	}
	// 保存摄像机的位置 [8/17/2012 Han]
	if (pDoc->m_pQslimMesh != NULL)
	{
		pDoc->m_pQslimMesh->eye_pos[0] = eye.x;
		pDoc->m_pQslimMesh->eye_pos[1] = eye.y;
		pDoc->m_pQslimMesh->eye_pos[2] = eye.z;
		pDoc->m_pQslimMesh->bRotated = false;
	}
	if (pDoc->m_pOriginalQMesh != NULL)
	{
		pDoc->m_pOriginalQMesh->eye_pos[0] = eye.x;
		pDoc->m_pOriginalQMesh->eye_pos[1] = eye.y;
		pDoc->m_pOriginalQMesh->eye_pos[2] = eye.z;
		pDoc->m_pOriginalQMesh->bRotated = false;
	}


	// set up the arcball using the current projection matrix
	SPHERE_RADIUS =/* 5.0f; */pDoc->GetModelRadius();// abs(eye.z) / 3.0f/*1.1f*/;
	SPHERE_RADIUS = (SPHERE_RADIUS < 0.5f)?0.5f:SPHERE_RADIUS;
	arcball_setzoom( SPHERE_RADIUS, eye, up );

	//glOrtho(0.0f, nWinWidth, nWinHeight, 0.0f, -1.0f, 1.0f);			// Create Ortho 640x480 View (0,0 At Top Left)
	/* from here on we're setting modeling transformations */
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// reset light pos [5/25/2012 Han]
	m_lightPos[0] = eye.x/*pDoc->GetModelRadius()*/;
	m_lightPos[1] = eye.y/*pDoc->GetModelRadius()*/;
	m_lightPos[2] = eye.z /*25*//*pDoc->GetModelRadius()*//*eye.z*/;
	float ld = 25.f / GetViewDistance();
	for(int i = 0; i < 3; i++)
		m_lightPos[i] *= ld;
	glLightfv(GL_LIGHT0, GL_POSITION, m_lightPos);

	Invalidate();
	// 使用view的更新代码 [2/16/2012 Han]
	//if (pDoc->m_mType == CMeshSimpDoc::CurrentModelType::VIEW_SELECT)
	//	pDoc->m_viewMyRender->resizeGL(nWidth, nHeight);
	//else if (pDoc->m_mType == CMeshSimpDoc::CurrentModelType::VIEW_SELECT_OUR)
	//	pDoc->m_viewMyRenderOur->resizeGL(nWidth, nHeight);
}

void CMeshSimpView::AdjustView()
{
	AdjustView(width_, height_);
}


void CMeshSimpView::OnEditUnselectall()
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	if (pDoc->m_pProgMesh)
	{
		list<int>::iterator selIter;

		for (selIter = m_selIndexList.begin(); selIter != m_selIndexList.end(); ++selIter) 
		{
			triangle &t = pDoc->m_pProgMesh->getTri(*selIter);
			t.setSel(false);
		}
		m_bSelAny = false;
		m_selIndexList.clear();
		pDoc->changeSimplificationAlgorithm();
		pDoc->m_bSelCalc = false;
		Invalidate();
	}

}

void CMeshSimpView::OnUpdateEditUnselectall(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(m_bSelAny);
}

// 定时重绘 [7/5/2011 Han Honglei]
void CMeshSimpView::OnTimer(UINT_PTR nIDEvent)
{
	// test 暂时不进行定时重绘
	//Invalidate();

	CView::OnTimer(nIDEvent);
}

void CMeshSimpView::OnUpdateStatus(CCmdUI *pCmdUI)
{
	// TODO: Add your command update UI handler code here
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	if (pDoc->m_pProgMesh)
	{
		if (pDoc->m_pProgMesh->isCollapsing)
		{
			CString text;
			text.Format("%d iteration finished", pDoc->m_pProgMesh->collapseIterNum);
			pCmdUI->SetText(text);
		}
		else if (pDoc->m_pProgMesh->bCollapsed)
			pCmdUI->SetText("Geometry has been collapsed!");
		if (pDoc->m_pProgMesh->isSkeling)
			pCmdUI->SetText("Simpling collapsed geometry!");
		else if (pDoc->m_pProgMesh->bSimplified)
		{
			CString text;
			text.Format("Got the skeleton, there are %d nodes", pDoc->m_pProgMesh->simplifiedVertexRec.size());
			pCmdUI->SetText(text);
		}

	}
}

void CMeshSimpView::OnUpdateStatus2(CCmdUI *pCmdUI)
{
	// TODO: Add your command update UI handler code here
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	int f, v;
	if (pDoc->GetValidMeshInfo(&f, &v))
	{
		CString text;
		text.Format("Tri:%d,Ver:%d", f, v);
		pCmdUI->SetText(text);
	}
	else
		pCmdUI->SetText("No model loaded");
}

void CMeshSimpView::OnUpdateStatus3(CCmdUI *pCmdUI)
{
		CString text;
		float dis = sqrt(pow(eye.x-centre.x, 2)+pow(eye.y-centre.y,2)+ pow(eye.z-centre.z,2));
		float tilt = asin(eye.y/dis);
		float yaw = acos(eye.x/dis);
		text.Format("View Dist:%.3f;Tilt:%.2f;Yaw:%.2f", dis, tilt, yaw);
		pCmdUI->SetText(text);

}
void CMeshSimpView::OnNextColl()
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	if (pDoc->m_pProgMesh)
		pDoc->m_pProgMesh->ChangeColl(true);
	Invalidate();
}


void CMeshSimpView::OnPrevColl()
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	if (pDoc->m_pProgMesh)
		pDoc->m_pProgMesh->ChangeColl(false);
	Invalidate();
}



void CMeshSimpView::OnShadingLighting()
{
	bLighting = !bLighting;
	Invalidate();
}


void CMeshSimpView::OnUpdateShadingLighting(CCmdUI *pCmdUI)
{
	pCmdUI->Enable();
	pCmdUI->SetCheck(bLighting);
}


void CMeshSimpView::OnShadingSkeleton()
{
	bDrawSkel = !bDrawSkel;
	Invalidate();
}


void CMeshSimpView::OnUpdateShadingSkeleton(CCmdUI *pCmdUI)
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	pCmdUI->Enable(pDoc->m_pProgMesh != NULL);
	if (pDoc->m_pProgMesh)
		pCmdUI->Enable(pDoc->m_pProgMesh->bSimplified);

	pCmdUI->SetCheck(bDrawSkel);
}


void CMeshSimpView::OnShadingTransparant()
{
	bTrans = !bTrans;
	if (bTrans)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	else
		glDisable(GL_BLEND);

	Invalidate();
}


void CMeshSimpView::OnUpdateShadingTransparant(CCmdUI *pCmdUI)
{
	pCmdUI->Enable();
	pCmdUI->SetCheck(bTrans);
}


void CMeshSimpView::OnShadingDrawskelmap()
{
	drawMeshType = PMesh::SKEL_MAP;
	Invalidate();
}


void CMeshSimpView::OnUpdateShadingDrawskelmap(CCmdUI *pCmdUI)
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	pCmdUI->Enable(pDoc->m_pProgMesh != NULL);
	if (pDoc->m_pProgMesh)
		pCmdUI->Enable(pDoc->m_pProgMesh->bSimplified);

	pCmdUI->SetCheck(drawMeshType == PMesh::SKEL_MAP);
}

void CMeshSimpView::OnShadingDrawcolldist()
{
	drawMeshType = PMesh::COLL_DIST;
	Invalidate();
}

void CMeshSimpView::OnUpdateShadingDrawcolldist(CCmdUI *pCmdUI)
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	pCmdUI->Enable(pDoc->m_pProgMesh != NULL);
	if (pDoc->m_pProgMesh)
		pCmdUI->Enable(pDoc->m_pProgMesh->bSimplified);

	pCmdUI->SetCheck(drawMeshType == PMesh::COLL_DIST);
}

void CMeshSimpView::OnShadingDrawsegment()
{
	drawMeshType = PMesh::SEGMENT;
	Invalidate();
}

void CMeshSimpView::OnUpdateShadingDrawsegment(CCmdUI *pCmdUI)
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	pCmdUI->Enable(pDoc->m_pProgMesh != NULL);
	if (pDoc->m_pProgMesh)
		pCmdUI->Enable(pDoc->m_pProgMesh->HasSegment());
	pCmdUI->SetCheck(drawMeshType == PMesh::SEGMENT);
}

void CMeshSimpView::OnUpdateNextColl(CCmdUI *pCmdUI)
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	pCmdUI->Enable(pDoc->m_pProgMesh != NULL);
	if (pDoc->m_pProgMesh)
		pCmdUI->Enable(pDoc->m_pProgMesh->bCollapsed);
}


void CMeshSimpView::OnUpdatePrevColl(CCmdUI *pCmdUI)
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	pCmdUI->Enable(pDoc->m_pProgMesh != NULL);
	if (pDoc->m_pProgMesh)
		pCmdUI->Enable(pDoc->m_pProgMesh->bCollapsed);
}


void CMeshSimpView::OnShadingDrawmesh()
{
	bDrawMesh = !bDrawMesh;
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	if (pDoc->m_viewMyRenderOur)
		pDoc->m_viewMyRenderOur->showing_type = bDrawMesh?pDoc->m_viewMyRenderOur->SHOWING_MODEL:pDoc->m_viewMyRenderOur->SHOWING_VIEW_SPHERE_MAP;

	if (bDrawMesh)
		showing_type = SHOWING_MODEL;
	else
		showing_type = SHOWING_MESH_SALIENCY;

	Invalidate();
}


void CMeshSimpView::OnUpdateShadingDrawmesh(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(showing_type == SHOWING_MODEL);
}


void CMeshSimpView::OnViewSelection()
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	pDoc->ViewSelection();

}


void CMeshSimpView::OnUpdateViewSelection(CCmdUI *pCmdUI)
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	pCmdUI->Enable(pDoc->m_isLoaded);
}

// 获得和opengl渲染投影矩阵相关的参数
void CMeshSimpView::GetGLViewParam(float * dist, float * dpi, float * projPlaneDist)
{
	*dist = sqrt(pow(eye.x-centre.x, 2)+pow(eye.y-centre.y,2)+ pow(eye.z-centre.z,2));
	*projPlaneDist = m_fNearClip;
	
	//所有像素数
	int pagecx=m_pDC->GetDeviceCaps(HORZRES);
	int pagecy=m_pDC->GetDeviceCaps(VERTRES);

	//即每英寸点数
	short cxInch = m_pDC->GetDeviceCaps(LOGPIXELSX);
	short cyInch = m_pDC->GetDeviceCaps(LOGPIXELSY);

	*dpi = cxInch;
	//*dpi = height_;
}

//  [3/29/2012 Han]
double CMeshSimpView::Get1PixelInMeter(float fpPos[3])
{
	double hC = m_fNearClip *tan(PI*(m_fov/2.0)/180);
	// 将pt单位表示的长度做一定的近似，转换为像素单位为整数的值
	//1in = 2.54cm = 25.4 mm = 72pt  = ( dpi ) px
	//float fUnit = 72.0/r;
	//float heightInPt = height_ * fUnit;
	//float hR1 = 0.0254*heightInPt/(2.0);
	//float hR2 = 0.0254*height_/(2.0*r);

	//// 真实dpi换算1080宽度像素*0.254   24.9
	//float scH = 0.0254*1080/r;
	//float scW = 0.0254*1920/r;

	double pixelIn3D = 2.0*hC / height_;
	return pixelIn3D;
}



void CMeshSimpView::OnShadingResetview()
{
	resetOrientation();
}


void CMeshSimpView::OnShadingDrawViewSphere()
{
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (showing_type == SHOWING_VIEW_SPHERE_MAP)
		showing_type = SHOWING_MODEL;
	else
	{
		showing_type = SHOWING_VIEW_SPHERE_MAP;
		pDoc->UpdateSphereColorByImportance();
	}
	Invalidate();
}


void CMeshSimpView::OnShadingAll6views()
{
	resetOrientation();
	CMeshSimpDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	static int cv = 0;
	eye.x = 0;
	eye.y = 0;
	eye.z = 0;
	switch (cv)
	{
	case 0:
		eye.z = -1;
		break;
	case 1:
		eye.x = 1;
		break;
	case 2:
		eye.z = 1;
		break;
	case 3:
		eye.x = -1;
		break;
	case 4:
		eye.y = 1;
		break;
	case 5:
		eye.y = -1;
		break;
	default:
		break;
	}
	cv = (cv+1)%6;
	eye = eye*9*pDoc->GetModelRadius();
	AdjustView();
	Invalidate();
}


void CMeshSimpView::OnViewpointselectionUprightdir()
{
	bShowUprightDir = !bShowUprightDir;
	Invalidate();
}


void CMeshSimpView::OnShadingDrawnormal()
{
	showing_type = SHOWING_NORMALMAP;
	Invalidate();
}
