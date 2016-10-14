// MeshSimpDoc.h : interface of the CMeshSimpDoc class
//


#pragma once
#include "stdafx.h"
#include "ViewSelectionBenchmark/ViewSelect-MS/myrender-MS.h"
#include "ViewSelectionBenchmark/ViewSelect-OUR/myrender-OUR.h"
#include "jmspmesh/jmsmesh.h"
#include "jmspmesh/pmesh.h"
#include "stdmix.h"
#include <MxQSlim.h>
#include <MxSMF.h>
#include "MiniBall/Miniball_dynamic_d.h"

#define QSLIM_VERSION 2100
#define QSLIM_VERSION_STRING "2.1"

typedef MxDynBlock<MxPairContraction> QSlimLog;

class CMeshSimpDoc : public CDocument
{
protected: // create from serialization only
	CMeshSimpDoc();
	DECLARE_DYNCREATE(CMeshSimpDoc)

	// Attributes
public:
	bool	m_isLoaded;
	bool	m_bSelCalc;
	// Progressive jmsMesh
	PMesh* m_pProgMesh;
	void changeSimplificationAlgorithm();

protected:
	// Edge Collapse Options
	PMesh::EdgeCost m_edgemethod ;
	// file name
	char m_sFileName[256];
	char m_sModelName[256];
	// Triangle model
	jmsMesh* m_pMesh;

	// Operations
	void InitNewDoc();
	void DrawScene(CPoint upLeft, CPoint downRight);
	void loadMesh(char *fileName);

	// 将当前文件保存为ply格式的模型 [3/31/2012 Han]
	void SaveCurrMesh2PLY(CArchive& ar);
	void SaveSMF2PLY(CArchive& ar, MxStdModel* m);	// 将smf文件保存为ply
	void OutputCurrentCurvature(std::ofstream &f, double maxC[MSDM_SCALE], double minC[MSDM_SCALE], int Hnumber, double step);

	Miniball_C *m_miniBall;
	void GetMiniBall(MxStdModel *m);
	void align_center(MxStdModel *m);

	// Segmentation [7/12/2012 Han]
	bool LoadMeshSaliency(char *fn);

public:

	enum OutputFormat { SMF, PM, MMF, LOG, IV, VRML };
	enum CurrentModelType {NORMAL, QSLIM, VIEW_SELECT, VIEW_SELECT_OUR};
	enum SimpByMethod {NONE, SIMP_BY_DIST, SIMP_BY_MSDM};
	enum FilterRadius { PIXEL1, PIXEL2, CSF, CSF2, BADMETRIC};


	bool (*unparsed_hook)(char *, int, char*[], MxStdModel&); 

	////////////////////////////////////////////////////////////////////////
	//
	// Command line parsing and application initialization
	//
	void startup_and_input(int argc, char **argv);
	void output_ivrml(ostream& out, bool vrml=false);
	void setup_output();
	void cleanup_for_output();
	//void process_cmdline(int argc, char **argv);
	void slim_init();
	void slim_cleanup();
	void input_file(const char *);
	void defer_file_inclusion(char *);
	void include_deferred_files();
	void slim_history_callback(const MxPairContraction&,float);

	////////////////////////////////////////////////////////////////////////
	//
	// Output routines
	//
	bool select_output_format(const char *);
	void output_preamble();
	void output_current_model();
	void output_final_model();

	void output_iv(ostream&);
	void output_vrml(ostream&);
	void output_regressive_mmf(ostream&);
	void output_regressive_log(ostream&);
	void output_progressive_pm(ostream&);

	////////////////////////////////////////////////////////////////////////
	//
	// Other relevant things
	//
	const char *slim_copyright_notice;
	const char *slim_version_string;
	void slim_print_banner(ostream&);



	void startup_and_input(char* filename , int otherType = -1);
	static MxDynBlock<MxEdge> *target_edges;
	// qslim使用的模型
	MxStdModel *m_pQslimMesh;
	MxStdModel *m_pOriginalQMesh;
	float GetModelRadius();
	// view selection 所用的网格类
	MyRender_MS *m_viewMyRender;
	MyRender_OUR *m_viewMyRenderOur;	// 使用YLM的方法
	// 将简化以后的模型转换为 [2/16/2012 Han]
	bool UpdateCurModel();
	CurrentModelType m_mType;
	bool ViewSelection();
	FILE* outputFile;

	bool GetValidMeshInfo(int *validF, int *validV);		// 获取当前模型的面片和顶点信息

	void DrawViewSphere(float radius = 1.0f);

	FilterRadius m_filterRadius;

	double GetFilterRadiusInNearclip();
	void DrawUprightDir();
	ViewPoint *m_baseViewpoint;
	ViewPoint *m_clusterdBestView;

protected:
	// 使用最差视点决定模型竖直方向 [10/9/2013 Han]
	int* CandBaseScore(const vector<vec3*> &candBase, const vec3 &bestView, int &bestBaseIdx);
	float Subdivide(vec3 *candBase,int baseID, ViewSelectType FilterType, bool bSmallest = true);
	ViewPoint* BestBaseView(float k, vector<vec3*> &candBase, int &nBestBase,ViewSelectType FilterType);
	void Upright(ViewSelectType FilterType);
	void UpdateViewsImportance(bool bFromView2Model = true);
	void DrawArrow(vec3 posA, vec3 posB, float *color = NULL);
	// 目前找到的最差视点ID [10/23/2012 Han]
	int m_nWorstViewID;
	// 目前找到的最优视点编号 [10/23/2012 Han]
	int m_nBestViewID;

	////////////////////////////////////////////////////////////////////////
	//
	// Globally visible (configuration) variables
	//
	unsigned int face_target;
	bool will_use_fslim;
	int placement_policy;
	double boundary_weight;
	int weighting_policy;
	bool will_record_history;
	double compactness_ratio;
	double meshing_penalty;
	bool will_join_only;
	bool be_quiet;
	OutputFormat output_format;
	char *output_filename;
	ostream *output_stream;

	// segment [7/10/2012 Han]
	//////////////////////////////////////////////////////////////////////////
	// 计算各个控制参数的依据，k：最大saliency和最小saliency的1/255, min: 模型面片数或者顶点数的0.0004, [7/25/2012 Han]
	//////////////////////////////////////////////////////////////////////////
	float  sigma;
	float k ;
	float min_size;
	bool segFace;

	MxSMFReader *smf;
	MxStdModel *m_orig;
	MxQSlim *slim;
	MxEdgeQSlim *eslim;
	MxFaceQSlim *fslim;
	QSlimLog *history;
	//////////////////////////////////////////////////////////////////////////
	vector<ViewPoint*> viewpoint_candidates; // candidate on the real viewing sphere
	vector<ViewSphereFace> vf ;
	MxStdModel *m_pViewPhereMdl;
	void EqualDistributeViews(float r = 1.0f, int itera = 3/*3*/);
public:
	void UpdateSphereColorByImportance();

	// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CMeshSimpDoc)
public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
	void LoadExcuteAllFiles(CString tp, int action=1);
	//}}AFX_VIRTUAL

	// Implementation
public:
	virtual ~CMeshSimpDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
protected:
	//{{AFX_MSG(CMeshSimpDoc)
	// NOTE - the ClassWizard will add and remove member functions here.
	//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()

public:
	afx_msg void OnResampleAdd1tri();
	afx_msg void OnResampleIncreaseby5();
	afx_msg void OnResampleReduceby5();
	afx_msg void OnResampleRemove1tri();
	afx_msg void OnUpdateSimpMelax(CCmdUI *pCmdUI);
	afx_msg void OnUpdateSimpQuadric(CCmdUI *pCmdUI);
//	afx_msg void OnUpdateSimpShorttest(CCmdUI *pCmdUI);
	afx_msg void OnUpdateSimpQuadricWeighted(CCmdUI *pCmdUI);
	afx_msg void OnUpdateSimpShorttest(CCmdUI *pCmdUI);
	afx_msg void OnSimpMelax();
	afx_msg void OnSimpQuadricWeighted();
	afx_msg void OnSimpQuadric();
	afx_msg void OnSimpShorttest();
	afx_msg void OnSkeletonizerGetskeleton();
	afx_msg void OnSkelColl();
	afx_msg void OnUpdateSkelColl(CCmdUI *pCmdUI);
	afx_msg void OnUpdateSkeletonizerGetskeleton(CCmdUI *pCmdUI);
	afx_msg void OnSimplificationalgorithmQuadricweightedbysegment();
	afx_msg void OnUpdateSimplificationalgorithmQuadricweightedbysegment(CCmdUI *pCmdUI);
	afx_msg void OnNormalizeSize();
	void AdjustLOD(float fDist);
	SimpByMethod m_SimpByMethod;
	afx_msg void OnSimplificationalgorithmQuadricweightedbymotionblur();
	afx_msg void OnUpdateSimplificationalgorithmQuadricweightedbymotionblur(CCmdUI *pCmdUI);
	afx_msg void OnSimplificationalgorithmSimpbydistance();
	afx_msg void OnUpdateSimplificationalgorithmSimpbydistance(CCmdUI *pCmdUI);
	bool LoadSMF(char *filename);
	bool LoadOFF2SMF(char *filename);
	bool LoadViewSelectModel(char *filename);
	int LoadSegFiles(char *fileName);
	//bool DisplayMesh(PMesh::DrawMeshType drawType = PMesh::NONE, bool bSmooth = true);
	afx_msg void OnViewselectionMeshcurvature();
	afx_msg void OnViewselectionMeshsaliency();
	afx_msg void OnViewselectionViewdmeshcurvature();
	afx_msg void OnViewselectionViewplanecurvature();
	afx_msg void OnSimplificationalgorithmSimpbymsdm2();
	afx_msg void OnUpdateSimplificationalgorithmSimpbymsdm2(CCmdUI *pCmdUI);
	afx_msg void OnViewSlectionByDistance();

	afx_msg void OnGetVisibleFaces();
	void ViewSlectionByLOD(float dis);
	void ViewSlectionByFilter(float dis);
	void ViewSlectionByViewEntropy(float d);
	void ViewSlectionBySaliencySegment(float d);
	void ViewSlectionByViewImportance(float d, ViewSelectType FilterType, bool bSimpMesh = true);
	void SaveRenderPics(const point *bestViewPos,float d, float maxImportance,wchar_t* fnhead, int FilterType, int Rank = 1, bool bSimpMesh=true);
	void SortSaveCandidateViews(float d, ViewSelectType FilterType, FILE *outputViewSelectionFile, bool bVisibleFace, bool bSimpMesh);
	float ViewQuality(vec3 view,ViewSelectType FilterType, bool bVisibleFace = true, bool bSimpMesh = true);
	afx_msg void OnEditTestonly();
//	afx_msg void OnUpdateFilterradius1pixel(CCmdUI *pCmdUI);
	afx_msg void OnFilterradius2pixel();
	afx_msg void OnFilterradiusCsf();
	afx_msg void OnFilterradius2();
	afx_msg void OnFilterradius1pixel();
	afx_msg void OnEditTestonly2();
	afx_msg void OnEditCalcsaliency();
	afx_msg void OnEditSmoothsaliency();
	afx_msg void OnEditSegbysaliency();
	//afx_msg void OnEditDotestviewselect();
	afx_msg void OnEditSaliencyvs();
	afx_msg void OnEditViewentropyvs();
	afx_msg void OnEditSemanticdrivenvs();
	afx_msg void OnResampleReducebycsf();
	afx_msg void OnResampleReduceuserdef();
	afx_msg void OnEditVisiblefaces();
	afx_msg void OnEditVisiblevertexes();
	afx_msg void OnShadingOutputallcandviews();
	afx_msg void OnViewselectionGetcandviewscores();
	afx_msg void OnViewpointselectionCe();
	afx_msg void OnViewpointselectionMeancurvature();
	void ViewpointSelectionBaseCurvature();
	afx_msg void OnViewpointselectionSaliencyentropy();
	afx_msg void OnViewpointselectionMeancentropy();
	afx_msg void OnViewselectionMeshcurvatureGeo();
	afx_msg void OnViewpointselectionNbestviewsmce();
	afx_msg void OnViewpointselectionSaliencyatlas();
	afx_msg void OnEditLoadsegdata();
	afx_msg void OnEditSemantic2009();
	afx_msg void OnUprightCurshannon();
	afx_msg void OnUprightMcs();
	afx_msg void OnUprightMeancur();
	afx_msg void OnUprightMeshsaliency();
	afx_msg void OnUprightMiniarea();
	afx_msg void OnUprightViewshannon();
	afx_msg void OnViewpointselectionMaxarea();
	afx_msg void OnViewpointselectionViewentropy();
	afx_msg void OnUprightAll();
	afx_msg void OnUprightBasecur();
};


