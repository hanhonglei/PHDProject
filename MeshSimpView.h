// MeshSimpView.h : interface of the CMeshSimpView class
//


#pragma once

// Model orientation
const float ORIG_ELEVATION = 0.0f;
const float ORIG_AZIMUTH = 0.0f;
const float ORIG_DIST = 3.0f;

const float MIN_DISTANCE = 0.1f;
const float MAX_DISTANCE = 100.0f;
#include "arcball/arcball.h"
#include "Shaders/shader.h"


class CMeshSimpView : public CView
{
protected: // create from serialization only
	CMeshSimpView();
	DECLARE_DYNCREATE(CMeshSimpView)

	// Attributes
protected:
	// 保存鼠标的位置 [2/17/2012 Han]
	CPoint m_mousePos;
	// 指向显示设备
	CClientDC *m_pDC;

	enum MOUSEMODE{NORMAL, ROTATE, SELECT};

	MOUSEMODE m_mouseMode;
	// width, height of window
	int width_;
	int height_;

	// previous window width, height
	int oldWidth_;
	int oldHeight_;


	// Full Screen display indicator
	bool bFullScreen_;

public:
	// Use Gouraud shading?
	bool bSmooth_;
protected:
	void AdjustView(int nWidth, int nHeight);
	//void DrawSimplifiedVertices();


	bool m_bArcBallStart;

	bool bDrawSkel;
	bool bLighting;
	PMesh::DrawMeshType drawMeshType;
	bool bTrans;
	bool bDrawMesh;

	int activeTriNum;

	float m_fNearClip;		// distance of near clip plane to the viewer
	float m_fFarClip;		// farest clip plane to the viewer
	float m_fov;
	//shaders for the models
	Shader		m_phongShader;			//shader for the roof (phong shaded)

public:
	void AdjustView();
	float GetViewDistance();
	bool m_bSelAny;
	list<int> m_selIndexList; 

	vec_arcball eye;
	// Fill in triangles (or use wireframe mode?)
	bool bFill_;

	enum RenderType{
		SHOWING_WHOLE_SKELETON, 
		SHOWING_PRIMARY_SKELETON,
		SHOWING_FINAL_SKELETON,
		SHOWING_MESH_SALIENCY,
		SHOWING_MODEL,
		SHOWING_MODEL_WIRE,
		SHOWING_MODEL_VOXELS,
		SHOWING_SEGMENTATION,
		SHOWING_VIEW_DEPENDENT_CURVATURE,
		SHOWING_VIEW_SPHERE_MAP,
		VISIBLE_TEST,
		SHOWING_NORMALMAP,
		SHOWING_CURVATURE_QUANTITY,
	} showing_type;

	bool bShowUprightDir;


	// Operations
public:
	CMeshSimpDoc* GetDocument();
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
	void DisplayModel();			// display models when setting everything ok

	void InitNewDoc();

	double Get1PixelInMeter(float fpPos[3] = 0);


	// Implementation
public:
	BOOL bSetPixelFormat();
	void IniOpenGL();
	virtual ~CMeshSimpView();
	void resetOrientation() ;
	void SelectTri(CPoint point, float radius);

	// Generated message map functions
protected:
	//{{AFX_MSG(CMeshSimpView)
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnDestroy();
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
//	afx_msg void OnFillWireframe();
//	afx_msg void OnMenuSelect(UINT nItemID, UINT nFlags, HMENU hSysMenu);
	afx_msg void OnUpdateShadingSmooth(CCmdUI *pCmdUI);
	afx_msg void OnShadingSmooth();
	afx_msg void OnFillWireframe();
	afx_msg void OnUpdateFillWireframe(CCmdUI *pCmdUI);
//	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnKeyUp(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnEditUnselectall();
	afx_msg void OnUpdateEditUnselectall(CCmdUI *pCmdUI);
	afx_msg void OnTimer(UINT_PTR nIDEvent);
	afx_msg void OnUpdateUIState(UINT /*nAction*/, UINT /*nUIElement*/);
	afx_msg void OnUpdateStatus(CCmdUI *pCmdUI);
	afx_msg void OnUpdateStatus2(CCmdUI *pCmdUI);
	afx_msg void OnUpdateStatus3(CCmdUI *pCmdUI);
	afx_msg void OnNextColl();
	afx_msg void OnPrevColl();
	afx_msg void OnShadingDrawcolldist();
	afx_msg void OnUpdateShadingDrawcolldist(CCmdUI *pCmdUI);
	afx_msg void OnShadingLighting();
	afx_msg void OnUpdateShadingLighting(CCmdUI *pCmdUI);
	afx_msg void OnShadingSkeleton();
	afx_msg void OnUpdateShadingSkeleton(CCmdUI *pCmdUI);
	afx_msg void OnShadingTransparant();
	afx_msg void OnUpdateShadingTransparant(CCmdUI *pCmdUI);
	afx_msg void OnShadingDrawskelmap();
	afx_msg void OnUpdateShadingDrawskelmap(CCmdUI *pCmdUI);
	afx_msg void OnUpdateNextColl(CCmdUI *pCmdUI);
	afx_msg void OnUpdatePrevColl(CCmdUI *pCmdUI);
	afx_msg void OnShadingDrawsegment();
	afx_msg void OnUpdateShadingDrawsegment(CCmdUI *pCmdUI);
	afx_msg void OnShadingDrawmesh();
	afx_msg void OnUpdateShadingDrawmesh(CCmdUI *pCmdUI);
	afx_msg void OnViewSelection();
	afx_msg void OnUpdateViewSelection(CCmdUI *pCmdUI);
	void GetGLViewParam(float * dist, float * dpi, float * projPlaneDist);
	void SaveCurrentModelPic(wchar_t fn[128]);
	afx_msg void OnShadingResetview();
	afx_msg void OnShadingDrawViewSphere();
	afx_msg void OnShadingAll6views();
	afx_msg void OnViewpointselectionUprightdir();
	afx_msg void OnShadingDrawnormal();
};

//#ifndef _DEBUG  // debug version in MeshSimpView.cpp
//inline CMeshSimpDoc* CMeshSimpView::GetDocument() const
//   { return reinterpret_cast<CMeshSimpDoc*>(m_pDocument); }
//#endif

