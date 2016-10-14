// MeshSimpDoc.cpp : implementation of the CMeshSimpDoc class
//

#include "stdafx.h"
#include "MeshSimpDoc.h"
#include <stdmix.h>
#include <MxTimer.h>
#include "qslim.h"


#include "MeshSimp.h"

#include "MainFrm.h"
#include "MeshSimpView.h"
#include "MSDM2/MSDM2.h"
#include "ViewSelectionBenchmark/ViewSelect-OUR/ppm.h"
#include "Segment/segment-image.h"
#include "SegmentParam.h"
//#include "vld.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif
const int REDUCE_TRI_PERCENT = 5;	// when page up/page down, inc/dec # tris by this percent 
// (percent applies to # of tris in *original* mesh.)
const int NUM_PAGEUPDN_INTERVALS = 100/REDUCE_TRI_PERCENT;
// CMeshSimpDoc
IMPLEMENT_DYNCREATE(CMeshSimpDoc, CDocument)
BEGIN_MESSAGE_MAP(CMeshSimpDoc, CDocument)
	ON_COMMAND(ID_RESAMPLE_ADD1TRI, &CMeshSimpDoc::OnResampleAdd1tri)
	ON_COMMAND(ID_RESAMPLE_INCREASEBY5, &CMeshSimpDoc::OnResampleIncreaseby5)
	ON_COMMAND(ID_RESAMPLE_REDUCEBY5, &CMeshSimpDoc::OnResampleReduceby5)
	ON_COMMAND(ID_RESAMPLE_REMOVE1TRI, &CMeshSimpDoc::OnResampleRemove1tri)
	ON_UPDATE_COMMAND_UI(ID_SIMP_MELAX, &CMeshSimpDoc::OnUpdateSimpMelax)
	ON_UPDATE_COMMAND_UI(ID_SIMP_QUADRIC, &CMeshSimpDoc::OnUpdateSimpQuadric)
//	ON_UPDATE_COMMAND_UI(ID_SIMP_SHORTTEST, &CMeshSimpDoc::OnUpdateSimpShorttest)
ON_UPDATE_COMMAND_UI(ID_SIMP_QUADRIC_WEIGHTED, &CMeshSimpDoc::OnUpdateSimpQuadricWeighted)
ON_UPDATE_COMMAND_UI(ID_SIMP_SHORTTEST, &CMeshSimpDoc::OnUpdateSimpShorttest)
ON_COMMAND(ID_SIMP_MELAX, &CMeshSimpDoc::OnSimpMelax)
ON_COMMAND(ID_SIMP_QUADRIC_WEIGHTED, &CMeshSimpDoc::OnSimpQuadricWeighted)
ON_COMMAND(ID_SIMP_QUADRIC, &CMeshSimpDoc::OnSimpQuadric)
ON_COMMAND(ID_SIMP_SHORTTEST, &CMeshSimpDoc::OnSimpShorttest)
ON_COMMAND(ID_SKELETONIZER_GETSKELETON, &CMeshSimpDoc::OnSkeletonizerGetskeleton)
ON_COMMAND(ID_SKEL_COLL, &CMeshSimpDoc::OnSkelColl)
ON_UPDATE_COMMAND_UI(ID_SKEL_COLL, &CMeshSimpDoc::OnUpdateSkelColl)
ON_UPDATE_COMMAND_UI(ID_SKELETONIZER_GETSKELETON, &CMeshSimpDoc::OnUpdateSkeletonizerGetskeleton)
ON_COMMAND(ID_SIMPLIFICATIONALGORITHM_QUADRICWEIGHTEDBYSEGMENT, &CMeshSimpDoc::OnSimplificationalgorithmQuadricweightedbysegment)
ON_UPDATE_COMMAND_UI(ID_SIMPLIFICATIONALGORITHM_QUADRICWEIGHTEDBYSEGMENT, &CMeshSimpDoc::OnUpdateSimplificationalgorithmQuadricweightedbysegment)
ON_COMMAND(ID_NORMALIZE_SIZE, &CMeshSimpDoc::OnNormalizeSize)
ON_COMMAND(ID_SIMPLIFICATIONALGORITHM_QUADRICWEIGHTEDBYMOTIONBLUR, &CMeshSimpDoc::OnSimplificationalgorithmQuadricweightedbymotionblur)
ON_UPDATE_COMMAND_UI(ID_SIMPLIFICATIONALGORITHM_QUADRICWEIGHTEDBYMOTIONBLUR, &CMeshSimpDoc::OnUpdateSimplificationalgorithmQuadricweightedbymotionblur)
ON_COMMAND(ID_SIMPLIFICATIONALGORITHM_SIMPBYDISTANCE, &CMeshSimpDoc::OnSimplificationalgorithmSimpbydistance)
ON_UPDATE_COMMAND_UI(ID_SIMPLIFICATIONALGORITHM_SIMPBYDISTANCE, &CMeshSimpDoc::OnUpdateSimplificationalgorithmSimpbydistance)
ON_COMMAND(ID_VIEWSELECTION_MESHCURVATURE, &CMeshSimpDoc::OnViewselectionMeshcurvature)
ON_COMMAND(ID_VIEWSELECTION_MESHSALIENCY, &CMeshSimpDoc::OnViewselectionMeshsaliency)
ON_COMMAND(ID_VIEWSELECTION_VIEWDMESHCURVATURE, &CMeshSimpDoc::OnViewselectionViewdmeshcurvature)
ON_COMMAND(ID_VIEWSELECTION_VIEWPLANECURVATURE, &CMeshSimpDoc::OnViewselectionViewplanecurvature)
ON_COMMAND(ID_SIMPLIFICATIONALGORITHM_SIMPBYMSDM2, &CMeshSimpDoc::OnSimplificationalgorithmSimpbymsdm2)
ON_COMMAND(ID_VIEWSELECTION_VIEWDISTANCESELECTION, &CMeshSimpDoc::OnViewSlectionByDistance)
ON_UPDATE_COMMAND_UI(ID_SIMPLIFICATIONALGORITHM_SIMPBYMSDM2, &CMeshSimpDoc::OnUpdateSimplificationalgorithmSimpbymsdm2)

ON_COMMAND(ID_EDIT_TESTONLY, &CMeshSimpDoc::OnEditTestonly)
//ON_UPDATE_COMMAND_UI(ID_FILTERRADIUS_1PIXEL, &CMeshSimpDoc::OnUpdateFilterradius1pixel)
ON_COMMAND(ID_FILTERRADIUS_2PIXEL, &CMeshSimpDoc::OnFilterradius2pixel)
ON_COMMAND(ID_FILTERRADIUS_CSF, &CMeshSimpDoc::OnFilterradiusCsf)
ON_COMMAND(ID_FILTERRADIUS_2, &CMeshSimpDoc::OnFilterradius2)
ON_COMMAND(ID_FILTERRADIUS_1PIXEL, &CMeshSimpDoc::OnFilterradius1pixel)
ON_COMMAND(ID_EDIT_TESTONLY2, &CMeshSimpDoc::OnEditTestonly2)
ON_COMMAND(ID_EDIT_CALCSALIENCY, &CMeshSimpDoc::OnEditCalcsaliency)
ON_COMMAND(ID_EDIT_SMOOTHSALIENCY, &CMeshSimpDoc::OnEditSmoothsaliency)
ON_COMMAND(ID_EDIT_SEGBYSALIENCY, &CMeshSimpDoc::OnEditSegbysaliency)

//ON_COMMAND(ID_EDIT_DOTESTVIEWSELECT, &CMeshSimpDoc::OnEditDotestviewselect)
ON_COMMAND(ID_EDIT_SALIENCYVS, &CMeshSimpDoc::OnEditSaliencyvs)
ON_COMMAND(ID_EDIT_VIEWENTROPYVS, &CMeshSimpDoc::OnEditViewentropyvs)
ON_COMMAND(ID_EDIT_SEMANTICDRIVENVS, &CMeshSimpDoc::OnEditSemanticdrivenvs)
ON_COMMAND(ID_RESAMPLE_REDUCEBYCSF, &CMeshSimpDoc::OnResampleReducebycsf)
ON_COMMAND(ID_RESAMPLE_REDUCEUSERDEF, &CMeshSimpDoc::OnResampleReduceuserdef)
ON_COMMAND(ID_EDIT_VISIBLEFACES, &CMeshSimpDoc::OnEditVisiblefaces)
ON_COMMAND(ID_EDIT_VISIBLEVERTEXES, &CMeshSimpDoc::OnEditVisiblevertexes)
ON_COMMAND(ID_SHADING_OUTPUTALLCANDVIEWS, &CMeshSimpDoc::OnShadingOutputallcandviews)
ON_COMMAND(ID_VIEWSELECTION_GETCANDVIEWSCORES, &CMeshSimpDoc::OnViewselectionGetcandviewscores)
ON_COMMAND(ID_VIEWPOINTSELECTION_CE, &CMeshSimpDoc::OnViewpointselectionCe)
ON_COMMAND(ID_VIEWPOINTSELECTION_MEANCURVATURE, &CMeshSimpDoc::OnViewpointselectionMeancurvature)
ON_COMMAND(ID_VIEWPOINTSELECTION_SALIENCYENTROPY, &CMeshSimpDoc::OnViewpointselectionSaliencyentropy)
ON_COMMAND(ID_VIEWPOINTSELECTION_MEANCENTROPY, &CMeshSimpDoc::OnViewpointselectionMeancentropy)
ON_COMMAND(ID_VIEWSELECTION_MESHCURVATURE32857, &CMeshSimpDoc::OnViewselectionMeshcurvatureGeo)
ON_COMMAND(ID_VIEWPOINTSELECTION_NBESTVIEWSMCE, &CMeshSimpDoc::OnViewpointselectionNbestviewsmce)
ON_COMMAND(ID_VIEWPOINTSELECTION_SALIENCYATLAS, &CMeshSimpDoc::OnViewpointselectionSaliencyatlas)
ON_COMMAND(ID_EDIT_LOADSEGDATA, &CMeshSimpDoc::OnEditLoadsegdata)
ON_COMMAND(ID_EDIT_SEMANTIC2009, &CMeshSimpDoc::OnEditSemantic2009)
ON_COMMAND(ID_UPRIGHT_CURSHANNON, &CMeshSimpDoc::OnUprightCurshannon)
ON_COMMAND(ID_UPRIGHT_MCS, &CMeshSimpDoc::OnUprightMcs)
ON_COMMAND(ID_UPRIGHT_MEANCUR, &CMeshSimpDoc::OnUprightMeancur)
ON_COMMAND(ID_UPRIGHT_MESHSALIENCY, &CMeshSimpDoc::OnUprightMeshsaliency)
ON_COMMAND(ID_UPRIGHT_MINIAREA, &CMeshSimpDoc::OnUprightMiniarea)
ON_COMMAND(ID_UPRIGHT_VIEWSHANNON, &CMeshSimpDoc::OnUprightViewshannon)
ON_COMMAND(ID_VIEWPOINTSELECTION_MINAREA, &CMeshSimpDoc::OnViewpointselectionMaxarea)
ON_COMMAND(ID_VIEWPOINTSELECTION_VIEWENTROPY, &CMeshSimpDoc::OnViewpointselectionViewentropy)
ON_COMMAND(ID_UPRIGHT_ALL, &CMeshSimpDoc::OnUprightAll)
ON_COMMAND(ID_UPRIGHT_BASECUR, &CMeshSimpDoc::OnUprightBasecur)
END_MESSAGE_MAP()

// CMeshSimpDoc construction/destruction
MxDynBlock<MxEdge> *CMeshSimpDoc::target_edges = NULL;
CMeshSimpDoc::CMeshSimpDoc()

{
	m_isLoaded = FALSE;

	// Triangle model
	m_pMesh = NULL;

	// Progressive jmsMesh
	 m_pProgMesh = NULL;

	// Edge Collapse Options
	 m_edgemethod = PMesh::QUADRICTRI;
	 m_bSelCalc = false;

	 m_SimpByMethod = NONE;
	 /*使用Garland的Qslim方法进行处理所使用的变量*/
	 // Configuration variables and initial values
	 //
	 face_target = 0;
	 will_use_fslim = false;
	 placement_policy = /*MX_PLACE_ENDPOINTS*/MX_PLACE_OPTIMAL;
	 boundary_weight = 1000.0;
	 weighting_policy = MX_WEIGHT_UNIFORM/*MX_WEIGHT_AREA*/;
	 will_record_history = false;
	 compactness_ratio = 0.0;
	 meshing_penalty = 2.0;
	 will_join_only = false;
	 be_quiet = false;
	 output_format = SMF;
	 output_filename = NULL;
	 output_stream = NULL;

	 // Globally visible variables
	 //
	 smf = NULL;
	 m_pQslimMesh = NULL;
	 m_pOriginalQMesh = NULL;
	 m_pViewPhereMdl = NULL;
	 m_orig = NULL;
	 slim = NULL;
	 eslim = NULL;
	 fslim = NULL;
	 history = NULL;
	 target_edges = NULL;
	 m_baseViewpoint = NULL;
	 m_clusterdBestView = NULL;

	 slim_copyright_notice =
		 "Copyright (C) 1998-2002 Michael Garland.  See \"COPYING.txt\" for details.";

	 slim_version_string = "2.1";

	 //////////////////////////////////////////////////////////////////////////
	 m_viewMyRender  = NULL;
	 m_viewMyRenderOur = NULL;

	 m_mType = NORMAL;
	 //time_t tm;
	 //char fn[128];
	 //sprintf( fn, "Experiment%d.txt", time(&tm));

	 //outputFile = fopen(fn, "wt");
	 outputFile = NULL;

	 m_miniBall = NULL;

	 EqualDistributeViews();

	 m_filterRadius = CSF;

	 // segment [7/10/2012 Han]
//sigma: Used to smooth the input image before segmenting it.
//k: Value for the threshold function.
//min: Minimum component size enforced by post-processing.
	 //sigma = 0.5,k = 500, min = 20;
	 //////////////////////////////////////////////////////////////////////////
	 // 计算各个控制参数的依据，k：最大saliency和最小saliency的1/255, min: 模型面片数或者顶点数的0.0004, [7/25/2012 Han]
	 // play with these nums
	//////////////////////////////////////////////////////////////////////////
	 sigma = 0.5,k = 0.004/*1/255*/, min_size = 0.0005f, segFace = true;
}
CMeshSimpDoc::~CMeshSimpDoc()
{
	delete m_pMesh;
	delete m_pProgMesh;
	if (m_viewMyRender)
		delete m_viewMyRender;
	if (m_viewMyRenderOur)
		delete m_viewMyRenderOur;
	if (m_miniBall)
		delete m_miniBall;
	if (m_baseViewpoint)
		delete m_baseViewpoint;
	if (m_clusterdBestView)
		delete m_clusterdBestView;

	vector<ViewPoint*>::iterator iter;
	for (iter = viewpoint_candidates.begin(); iter != viewpoint_candidates.end(); ++iter)
		delete *iter;
	viewpoint_candidates.clear();
	slim_cleanup();
	delete m_pViewPhereMdl;
	if (outputFile != NULL)
		fclose(outputFile); // close the file
}
void CMeshSimpDoc::InitNewDoc()
{
	m_nWorstViewID = -1;
	m_nBestViewID = -1;
	if (m_baseViewpoint)
		delete m_baseViewpoint;
	m_baseViewpoint = NULL;
	if (m_clusterdBestView)
		delete m_clusterdBestView;
	m_clusterdBestView = NULL;
	m_edgemethod = PMesh::QUADRICTRI;
	m_bSelCalc = false;
	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 
	pView->resetOrientation(); 	
	pView->InitNewDoc();
}
BOOL CMeshSimpDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	return TRUE;
}
#ifdef _DEBUG
void CMeshSimpDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CMeshSimpDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG
// ↑ Doc的基本函数
//////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// ↓ 关于文件载入载出的函数
void CMeshSimpDoc::Serialize(CArchive& ar)
{
	CFile *file = ar.GetFile();
	int fileLen = (file->GetFilePath()).GetLength();
	char *sFileName = new char [fileLen+1];
	//char *ch = (char *)(file->GetFilePath()).GetBuffer((fileLen + 1));
	strcpy(sFileName, file->GetFilePath());

	CString szFile = ar.GetFile()->GetFileName();
	int   nPos=szFile.ReverseFind( '.'); 
	CString   szExt=szFile.Right(szFile.GetLength()-nPos-1); 
		// 保存当前的文件名,用于输出 [2/17/2012 Han]
	strcpy(m_sFileName, szFile);
	strcpy(m_sModelName, szFile);
	if (ar.IsStoring())
	{
		SaveCurrMesh2PLY(ar);
	}
	else
	{
//#define BAT_COMPUTE
#ifdef BAT_COMPUTE
		// 进行视点选择的批处理计算 [11/6/2013 Han]
		//LoadExcuteAllFiles(ar.GetFile()->GetFilePath(), 2);

		// 进行模型向上方向确定的批处理计算 [11/6/2013 Han]
		//LoadExcuteAllFiles("D:\\Projects\\obj_model\\smf\\TestDataSet\\normal", 1);
		//LoadExcuteAllFiles("D:\\Projects\\obj_model\\UprightedJiangBo", 1);
		LoadExcuteAllFiles("D:\\Projects\\UprightByWorstView\\SampleMeshes", 1);
		//LoadExcuteAllFiles("D:\\Projects\\obj_model\\smf\\2012_8_21ShapeAtNew", 1);
		//LoadExcuteAllFiles("D:\\Projects\\MeshProcessing\\PRINCETON SHAPE BENCHMARK\\psb_v1\\benchmark\\db", 1);
		//LoadExcuteAllFiles("D:\\Projects\\obj_model\\new_McGill_benchmark\\McGillDBAll", 1);
		//LoadExcuteAllFiles("D:\\Projects\\obj_model\\tsb\\tsb-0.0.1\\test", 1);
		// 以下为台式机模型存储位置
		//LoadExcuteAllFiles("C:\\3dmesh\\testdatas", 1);
		//LoadExcuteAllFiles("C:\\3dmesh\\psb_v1\\benchmark\\db", 1);
		//LoadExcuteAllFiles("C:\\3dmesh\\tsb\\tsb-0.0.1\\test", 1);
		//LoadExcuteAllFiles("C:\\3dmesh\\new_McGill_benchmark\\McGillDBAll", 1);		
		return;
#endif
		if (outputFile != NULL)
			fclose(outputFile); // close the file
		char fn[128];
		time_t tm;
		srand((unsigned)time(&tm));  //set seed;
		sprintf( fn, "Experiment_%s_%d.txt", m_sModelName, time(&tm));
		outputFile = fopen(fn, "wt");
		bool loaded = false;
		if (m_isLoaded)
		{
			if (m_pProgMesh != NULL)
			{
				if (szExt == "seg")
					loaded = m_pProgMesh->LoadSegFiles(sFileName);
				else if (szExt == "blr")
				{
					switch(m_mType){
					case NORMAL:
						loaded = m_pProgMesh->LoadBlurFiles(sFileName);
						break;
					case QSLIM:
						loaded = LoadSMF(sFileName);
						break;
					}
				}
			}
			else if (m_pQslimMesh != NULL)
			{
				loaded = LoadSegFiles(sFileName);
			}
		}
		// test 暂时删除，测试使用viewselection方法载入模型 [2/15/2012 Han]
		if (szExt == "ply" || szExt == "PLY")
		{
		double t1 = get_cpu_time();
			/*loaded =*/ loadMesh(sFileName);
			loaded = true;
		double t2 = get_cpu_time();
		fprintf(outputFile, "模型名称：%s，载入用时：%lf秒\n", sFileName,t2-t1);
		InitNewDoc();
			m_mType == NORMAL;
		}
		// test 测试使用viewSelction的方法载入模型 [2/15/2012 Han]
		//if (szExt == "ply" || szExt == "PLY" || szExt == "off" || szExt == "obj" || szExt == "OFF" || szExt == "OBJ")
		//{
		//	//if (m_viewMyRender == NULL)
		//	//	m_viewMyRender  = new MyRender_MS();
		//	//loaded = m_viewMyRender->load_model(std::wstring(sFileName, sFileName+fileLen));
		//	//m_viewMyRender->initializeGL();
		//	//m_mType = VIEW_SELECT;

		//	loaded = LoadViewSelectModel(sFileName);
		//	// 测试mesh saliency方法 [3/13/2012 Han]
		//	m_mType = VIEW_SELECT_OUR;
		//}
		else if (szExt == "smf" || szExt == "SMF" || szExt == "blr" || szExt == "BLR")
		{
			loaded = true;
			m_mType = QSLIM;
			time_t tm;
			srand((unsigned)time(&tm));  //set seed;
			loaded = LoadSMF(sFileName);
			m_mType = QSLIM;
		}
		else if (szExt == "off" || szExt == "OFF" )
		{
			loaded = true;
			m_mType = QSLIM;
			time_t tm;
			srand((unsigned)time(&tm));  //set seed;
			loaded = LoadOFF2SMF(sFileName);
			m_mType = QSLIM;
		}

		if (!loaded)
		{
				MessageBox(NULL, "请打开正确的文件类型", "ERROR", MB_OK);
				delete sFileName;
				return;
		}
		else
		{
			m_isLoaded = true;
			InitNewDoc();
		}
	}	

	// 释放指针
	delete sFileName;
}
void CMeshSimpDoc::LoadExcuteAllFiles(CString tp, int action)
{
	CFileFind fileFinder;
	CString filePath = tp + _T("//*.*");

	BOOL bFinished = fileFinder.FindFile(filePath);
	while(bFinished)  //每次循环对应一个类别目录
	{
		bFinished = fileFinder.FindNextFile();
		if(fileFinder.IsDirectory() && !fileFinder.IsDots())  //若是目录则递归调用此方法
		{
			LoadExcuteAllFiles(fileFinder.GetFilePath());
		}
		else  //再判断是否模型文件
		{
			//获取文件类型
			CString fileName = fileFinder.GetFileName();
			int dotPos=fileName.ReverseFind('.');
			CString fileExt=fileName.Right(fileName.GetLength()-dotPos);
			strcpy(m_sFileName, fileFinder.GetFilePath());
			strcpy(m_sModelName, fileName);

			if(fileExt != _T(".off")&&fileExt != _T(".OFF")
				&&fileExt != _T(".smf") && fileExt != _T(".SMF"))
				continue;

			if (outputFile != NULL)
				fclose(outputFile); // close the file
			char fn[128];
			time_t tm;
			srand((unsigned)time(&tm));  //set seed;
			sprintf( fn, "Experiment_%s_%d.txt", m_sModelName, time(&tm));
			outputFile = fopen(fn, "wt");

			if(fileExt == _T(".off")||fileExt == _T(".OFF"))  
				m_isLoaded = LoadOFF2SMF(m_sFileName);
			else if (fileExt == _T(".smf") || fileExt == _T(".SMF")/* || fileExt == _T("blr") || fileExt == _T("BLR")*/)
				m_isLoaded = LoadSMF(m_sFileName);
			else
				continue;

			if (!m_isLoaded)
				continue;
			m_mType = QSLIM;
			InitNewDoc();
			double tt = 0.0;
			switch (action)
			{
			case 1:
				OnUprightAll();

				break;
			case 2:
				// 使用MeanCurvature Entropy得到最优的N个视点 [10/23/2012 Han]
				TIMING(tt, OnViewpointselectionNbestviewsmce());//  [3/11/2013 Han]
				//TIMING(tt, OnViewpointselectionMeancentropy())

				fprintf(outputFile, "N Best viewpoints time:%lf\n", tt);
				// 使用6种方法进行投递VR的视点选择实验
				//OnViewSlectionByDistance();
				// 使用mean-curvature + entropy进行视点选择
				//OnViewpointselectionMeancentropy();
				break;
			default:
				break;
			}
		}
	} 
	if (outputFile != NULL)
		fclose(outputFile);
	fileFinder.Close();
}
// 使用YLM的方法 [2/22/2012 Han]
bool CMeshSimpDoc::LoadViewSelectModel(char *filename)
{
	bool loaded = false;
	if (m_viewMyRenderOur == NULL)
		m_viewMyRenderOur = new MyRender_OUR();
	int fileLen = strlen(filename);
	double tL;
	TIMING(tL, loaded = m_viewMyRenderOur->load_model_3d(std::wstring(filename, filename+fileLen)));
	// test
	m_viewMyRenderOur->showing_type = MyRender_OUR::Image_Type::SHOWING_MODEL;
	//m_viewMyRenderOur->showing_type = MyRender_OUR::Image_Type::SHOWING_VIEW_SPHERE_MAP;
	//Use YLM's method [2/22/2012 Han]
	m_viewMyRenderOur->initializeGL();
	fprintf(outputFile, "模型名称：%s，载入用时：%lf秒\n", filename,tL);
	// put buffer data to disk [6/8/2012 Han]
	fflush(outputFile);

	return loaded;
}
bool CMeshSimpDoc::LoadOFF2SMF(char *filename)
{
	char smfFile[256];
	CString szFilePath = filename;
	CString szFileName = m_sFileName;

	int   nPos=szFileName.ReverseFind( '.'); 
	CString   szExt=szFileName.Right(szFileName.GetLength()-nPos-1); 
	FILE* inFile;
	if (szExt == "off" || szExt == "OFF")
		strcpy(smfFile, filename);
	else
		return false;


	double input_time, init_time, slim_time, output_time;
	// 以便进行多次载入 [2/17/2012 Han]
	slim_cleanup();

	// Process command line and read input model(s)
	TIMING(input_time, startup_and_input(smfFile, 1));

	// Initial simplification process.  Collect contractions and build heap.
	TIMING(init_time, slim_init());

	// 将载入的模型复制一份，作为原始模型
	m_pOriginalQMesh = m_pQslimMesh->clone();

	fprintf(outputFile, "模型名称：%s，载入用时：%lf秒，QSlim初始化用时：%lf秒\n", m_sFileName,input_time, init_time);
	fprintf(outputFile, "包围球半径：%f,面片数：%d，顶点数：%d\n", GetModelRadius(), slim->valid_faces, slim->valid_verts);

	// put buffer data to disk [6/8/2012 Han]
	fflush(outputFile);

	// 将已经存在的辅助mesh结构删除 [8/24/2012 Han]
	if (m_viewMyRenderOur != NULL)
	{
		delete m_viewMyRenderOur;
		m_viewMyRenderOur = NULL;
	}
	return true;
}
// 载入smf格式的模型并初始化Qslim数据
bool CMeshSimpDoc::LoadSMF(char *filename)
{
	char smfFile[256];
	CString szFilePath = filename;
	CString szFileName = m_sFileName;

	int   nPos=szFileName.ReverseFind( '.'); 
	CString   szExt=szFileName.Right(szFileName.GetLength()-nPos-1); 
	FILE* inFile;
	if (szExt == "blr" || szExt == "BLR")
	{
		inFile = fopen(filename, "rt");
		char sn[256];
		fscanf(inFile, "%s", sn);					// 读取blr文件中其所对应的smf文件名
		CString smfFileName = filename;
		smfFileName.Replace(m_sFileName, sn);
		strcpy(smfFile, smfFileName);
	}
	else
		strcpy(smfFile, filename);


	double input_time, init_time, slim_time, output_time;
	// 以便进行多次载入 [2/17/2012 Han]
	slim_cleanup();

	// Process command line and read input model(s)
	TIMING(input_time, startup_and_input(smfFile));
	if (szExt == "blr" || szExt == "BLR")
	{
		m_pQslimMesh->LoadMBlur(inFile);	// 载入motion blur信息，并重新得到Q序列	
		fclose(inFile);
	}

	// Initial simplification process.  Collect contractions and build heap.
	TIMING(init_time, slim_init());

	// 将载入的模型复制一份，作为原始模型
	m_pOriginalQMesh = m_pQslimMesh->clone();

	fprintf(outputFile, "模型名称：%s，载入用时：%lf秒，QSlim初始化用时：%lf秒\n", m_sFileName,input_time, init_time);
	fprintf(outputFile, "包围球半径：%f,面片数：%d，顶点数：%d\n", GetModelRadius(), slim->valid_faces, slim->valid_verts);

	// put buffer data to disk [6/8/2012 Han]
	fflush(outputFile);

	// 将已经存在的辅助mesh结构删除 [8/24/2012 Han]
	if (m_viewMyRenderOur != NULL)
	{
		delete m_viewMyRenderOur;
		m_viewMyRenderOur = NULL;
	}

	double msdm_init_time;
	// test do not calc msdm2 Init
	//TIMING(msdm_init_time, MSDM2::MSDM2SimpInit(m_pOriginalQMesh, m_pQslimMesh));
	//fprintf(outputFile, "初始化msdm时间：%lf秒\n", msdm_init_time);

	// Decimate model until target is reached
	// 暂时不进行简化操作 [2/10/2012 Han]
	//TIMING(slim_time, slim->decimate(face_target));

	// Output the result
	//
	//TIMING(output_time, output_final_model());

	//if( !be_quiet )
	//{
	//	cerr << endl << endl;
	//	cerr << "+ Running time" << endl;
	//	cerr << "    Setup      : " << input_time << " sec" << endl;
	//	cerr << "    QSlim init : " << init_time << " sec" << endl;
	//	cerr << "    QSlim run  : " << slim_time << " sec" << endl;
	//	cerr << "    Output     : " << output_time << " sec" << endl;
	//	cerr << endl;
	//	cerr << "    Total      : "
	//		<< input_time+init_time+slim_time+output_time <<endl;
	//}
	//else
	//{
	//	cerr << slim->valid_faces << " " << init_time+slim_time << endl;
	//}

	//slim_cleanup();

	return true;
}
// 载入ply类型的模型,并使用ProgMesh数据结构来控制
void CMeshSimpDoc::loadMesh(char *fileName)
{
	char pszFileLocn[256] = {'\0'};


	//if (!m_pProgMesh) // 1st time through
	//{
	//	// load "plus sign" cursor
	//	SetClassLong(g_pWindow->getHWnd(), GCL_HCURSOR, (LONG) LoadCursor(NULL, IDC_CROSS));
	//}

	delete m_pMesh;
	m_pMesh = NULL; // not necessary, but a nice CYA habit
	delete m_pProgMesh;
	m_pProgMesh = NULL;

	m_pMesh = new jmsMesh(fileName);
	//strcpy(m_sFilePath, fileName);
	// 测试，不进行标准化，保持原来的大小
	//if (m_pMesh) m_pMesh->Normalize();// center mesh around the origin & shrink to fit

	m_pProgMesh = new PMesh(m_pMesh, m_edgemethod);

	// reset the position of the mesh
	//g_pWindow->resetOrientation();

	//g_pWindow->displayWindowTitle();

	// 载入mesh的时候需要做一些准备工作 [6/12/2011 Han Honglei]
	//public void OpenMeshFile()
	//{
	//	openFileDialog1.FileName = "";
	//	openFileDialog1.Filter = "jmsMesh files (*.obj)|*.obj";
	//	openFileDialog1.CheckFileExists = true;
	//	if (openFileDialog1.ShowDialog(this) == DialogResult.OK)
	//	{
	//		StreamReader sr = new StreamReader(openFileDialog1.FileName);
	//		jmsMesh m_pQslimMesh = new jmsMesh(sr);
	//		sr.Close();
	//		MeshRecord rec = new MeshRecord(openFileDialog1.FileName, m_pQslimMesh);
	//		meshes.Add(rec);
	//		currentMeshRecord = rec;
	//		TabPage page = new TabPage(rec.ToString()) {
	//			Tag = rec
	//		};
	//		tabControlModelList.TabPages.Add(page);
	//		tabControlModelList.SelectedTab = page;
	//		meshView1.SetModel(rec);
	//		propertyGridModel.SelectedObject = rec;
	//		PrintText("Loaded mesh " + openFileDialog1.FileName);
	//	}
	//}

}
// 将模型的当前curvature数据保存为外部文件
void CMeshSimpDoc::OutputCurrentCurvature(std::ofstream &f, double maxC[MSDM_SCALE], double minC[MSDM_SCALE], int Hnumber, double step)
{
	int *H = new int[Hnumber];

	f << "\n顶点数目：" << slim->valid_verts << "\n面片数目："<<slim->valid_faces << std::endl;
	float size[3];
	m_pQslimMesh->GetBoundSize(size);
	float max = 0;
	for (int i = 0; i < 3; i ++)
		max = max < size[i]? size[i] : max;
	double RadiusCurvature=MSDM2::mini_radius;

	// test 暂时不输出原始curvature ，只输出调整后的平均曲率
	//for (int i=0;i<MSDM_SCALE;i++)
	//{
	//	memset(H, 0, Hnumber*sizeof(int));
	//	// 计算不同半径下的模型曲率
	//	// 暂时不进行计算，因为在初始化时已经计算完毕
	//	//m_pQslimMesh->principal_curvature(true,RadiusCurvature*max, i);
	//	double step = 0.5/*(maxC[i] - minC[i]) / Hnumber*/;
	//	f << "Curvature Radius:" << RadiusCurvature*max << std::endl;
	//	f << "Max Curvature:" << maxC[i] << "\tMin Curvature:" << minC[i]<< "\tCurvature Step:" << step << std::endl;
	//	for(MxVertexID vID=0; vID<m_pQslimMesh->vert_count(); vID++)
	//	{
	//		if(!m_pQslimMesh->vertex_is_valid(vID) )
	//			continue;

	//		MxVertex &pVertex = m_pQslimMesh->vertex(vID);
	//		f << pVertex.KmaxCurv[i] << "\t";
	//		int R=(pVertex.KmaxCurv[i]/*-minC[i]*/)/step;
	//		R = R > Hnumber-1 ? Hnumber-1:R;
	//		R = R < 0 ? 0 : R;
	//		H[R]++;
	//	}
	//	f << std::endl;

	//	for (int j = 0; j < Hnumber; j++)
	//	{
	//		f << /*std::endl <<j*step << '~'<<(j+1)*step << "\t\t" <<*/ double(H[j])/slim->valid_verts << std::endl;
	//	}
	//	RadiusCurvature+=MSDM2::radius_step;
	//}

	//f << "======================================================" << std::endl;
	// 计算和模型尺寸无关的平均曲率 [5/9/2012 Han]
	RadiusCurvature=MSDM2::mini_radius;
	for (int i=0;i<MSDM_SCALE;i++)
	{
		memset(H, 0, Hnumber*sizeof(int));
		MSDM2::KmaxKmean(m_pQslimMesh,max, i);
		m_pQslimMesh->ComputeMaxMin();
		//f << "Curvature Radius:" << RadiusCurvature*max << std::endl;
		//f << "Max Mean Curvature:" << m_pQslimMesh->MaxCurv[i] << "\tMin Mean Curvature:" << m_pQslimMesh->MinCurv[i]<< "\tMean Curvature Step:" << step << std::endl;
		for(MxVertexID vID=0; vID<m_pQslimMesh->vert_count(); vID++)
		{
			if(!m_pQslimMesh->vertex_is_valid(vID) )
				continue;

			MxVertex &pVertex = m_pQslimMesh->vertex(vID);
			//f << pVertex.KmaxCurv[i] << "\t";
			int R=(pVertex.KmaxCurv[i]/*-minC[i]*/)/step;
			R = R > Hnumber-1 ? Hnumber-1:R;
			R = R < 0 ? 0 : R;
			H[R]++;
		}
		//f << std::endl;
		for (int j = 0; j < Hnumber; j++)
		{
			f << /*std::endl <<j*step << '~'<<(j+1)*step << "\t\t" <<*/ double(H[j])/slim->valid_verts << std::endl;
		}
		RadiusCurvature+=MSDM2::radius_step;
	}
	f << "Max Mean Curvature:" << m_pQslimMesh->MaxCurv[0] << "\tMin Mean Curvature:" << m_pQslimMesh->MinCurv[0]<< "\tMean Curvature Step:" << step << std::endl;
	f << "Max Dim: " << max;

	delete []H;
}
// 将当前模型保存为ply格式
void CMeshSimpDoc::SaveCurrMesh2PLY(CArchive& ar)
{
	switch (m_mType)
	{
	case NORMAL:
		if (m_pProgMesh != NULL)
			m_pProgMesh->SaveFile(ar);
		break;
	case QSLIM:
		if (m_pQslimMesh != NULL)
		{
			//CFile *file = ar.GetFile();
			//int fileLen = (file->GetFilePath()).GetLength();
			//char *sFileName = new char [fileLen+1];
			////char *ch = (char *)(file->GetFilePath()).GetBuffer((fileLen + 1));
			//strcpy(sFileName, file->GetFilePath());
			//sFileName[fileLen-1] = 'S';

			//if( output_stream != &cout && output_stream != NULL)
			//	delete output_stream;
			//output_stream = new ofstream(sFileName);
			//output_final_model();
			//delete []sFileName;

			SaveSMF2PLY(ar, m_pQslimMesh);	
		}

		break;
	}
}
// 载入模型分割信息,并赋予qmesh每个面片不同的id，按照类似atlas的方法给予不同的颜色
int CMeshSimpDoc::LoadSegFiles(char *fileName)
{
	FILE* inFile = fopen(fileName, "rt");
	int num = -1;
	if (inFile == NULL)
	{
		char pszError[_MAX_FNAME + 1];
		sprintf(pszError, "%s does not exist!\n", fileName);
		MessageBox(NULL, LPCSTR(pszError), NULL, MB_ICONEXCLAMATION);
		return num;
	}
	char tempStr[1024];
	m_pQslimMesh->MinImportanceF = FLT_MAX, m_pQslimMesh->MaxImportanceF = 0;

	for (unsigned int i = 0; i < m_pQslimMesh->face_count(); i++)
	{
		fscanf(inFile, "%s", tempStr);
		int tt = atoi(tempStr);
		m_pQslimMesh->face(i).from = tt;
		m_pQslimMesh->face(i).view_importance = tt;
		if(tt > m_pQslimMesh->MaxImportanceF)
			m_pQslimMesh->MaxImportanceF = tt;
		if (tt < m_pQslimMesh->MinImportanceF)
			m_pQslimMesh->MinImportanceF = tt;
	}
	num =  1+m_pQslimMesh->MaxImportanceF-m_pQslimMesh->MinImportanceF;
	fclose(inFile); // close the file
	m_pQslimMesh->VertexColor(3, m_pQslimMesh->MaxImportanceF, m_pQslimMesh->MinImportanceF);
	return num;
}
void CMeshSimpDoc::OnEditLoadsegdata()
{
	if (m_pQslimMesh == NULL)
	{
		MessageBox(NULL, "Please Load SMF Mesh firstly!", "ERROR", MB_OK);
		return;
	}
	CString FilePathName;
	CFileDialog dlg(TRUE);///TRUE为OPEN对话框，FALSE为SAVE AS对话框
	if(dlg.DoModal()==IDOK)
		FilePathName=dlg.GetPathName();
	else
		return ;
	int numSeg = LoadSegFiles(FilePathName.GetBuffer(FilePathName.GetLength()));
	if( numSeg!= -1)
	{
		if (!LoadMeshSaliency(NULL))
			return;
		m_pQslimMesh->saliencySegData.clear();
		m_pQslimMesh->nNumSeg = numSeg;
		m_pQslimMesh->UpdateSaliencyBySegment(true);
		m_pQslimMesh->VertexColor(3, m_pQslimMesh->MaxImportanceF, m_pQslimMesh->MinImportanceF);
	}
	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 
	pView->showing_type = pView->SHOWING_MESH_SALIENCY;
	// 进行视点选择 [8/14/2013 Han]
	ViewSlectionByViewImportance(pView->GetViewDistance(), /*SALIENCY_ENTROPY*/ /*SALIENCY*//*SEMANTIC_DRIVEN*/SEGMENT_SEMANTIC_ENTROPY, true);	
	pView->Invalidate();
}
// ↑是关于文件载入载出的函数
//////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// ↓是关于骨架提取的函数
UINT ThreadSkeletonizer(LPVOID param)
{
	PMesh *p = (PMesh*)param;
	CSkeOption* opt = new CSkeOption;
	p->isCollapsing = true;
	p->bCollapsed = false;
	//Program.displayProperty.MeshDisplayMode = DisplayProperty.EnumMeshDisplayMode.SmoothShaded;
	//Program.displayProperty.DisplaySelectedVertices = true;
	while (p->GeometryCollapse(opt));

	p->CalcCollapseDist();

	p->isCollapsing = false;
	p->bCollapsed = true;
	delete opt;
	return 0;
}
// 使用sig08算法得到模型的骨架 [6/13/2011 Han Honglei]
// 用于进行模型收缩，进而得到收缩到近似骨架的状态
void CMeshSimpDoc::OnSkeletonizerGetskeleton()
{
	// 重新开辟一个线程来处理收缩，这样可以在计算的时候实时看到结果 [7/5/2011 Han Honglei]
	if (m_pProgMesh)
		AfxBeginThread(ThreadSkeletonizer, m_pProgMesh);

	////Program.displayProperty.MeshDisplayMode = DisplayProperty.EnumMeshDisplayMode.SmoothShaded;
	////Program.displayProperty.DisplaySelectedVertices = true;


	//
	//if (opt.DisplayIntermediateMesh)
	//{
	//	Program.displayProperty.MeshDisplayMode = DisplayProperty.EnumMeshDisplayMode.None;
	//	Program.displayProperty.DisplaySelectedVertices = false;
	//}
	//if (opt.ApplyConnectivitySurgery)
	//{
	//	Simplification();
	//}
	//if (opt.ApplyEmbeddingRefinement)
	//{
	//	EmbeddingImproving();
	//}
	//if (opt.ApplyRootFinding)
	//{
	//	FindRootNode();
	//}
	//if (opt.ApplyConnectivitySurgery)
	//{
	//	DisplayOriginalMesh = true;
	//	Program.displayProperty.MeshDisplayMode = DisplayProperty.EnumMeshDisplayMode.TransparentSmoothShaded;
	//	Program.displayProperty.DisplaySelectedVertices = false;
	//}
	//else
	//{
	//	DisplayOriginalMesh = false;
	//	Program.displayProperty.MeshDisplayMode = DisplayProperty.EnumMeshDisplayMode.Wireframe;
	//	Program.displayProperty.DisplaySelectedVertices = false;
	//}
	//Program.PrintText("[Ready]");
	//Program.RefreshAllForms();
}
UINT ThreadSkeletonizerSimplification(LPVOID param)
{
	CSkeOption* opt = new CSkeOption;

	PMesh *p = (PMesh*)param;
	p->bSimplified = false;
	p->Simplification(opt);

	delete opt;

	return 0;
}
// 用于进行骨架抽取
void CMeshSimpDoc::OnSkelColl()
{
	if (m_pProgMesh)
		AfxBeginThread(ThreadSkeletonizerSimplification, m_pProgMesh);
}
// ↑是关于骨架提取的函数
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// ↓ 是得到或者修改模型信息的相关函数
// 得到当前模型的有效面片和顶点数目
bool CMeshSimpDoc::GetValidMeshInfo(int *validF, int *validV)
{
	bool succ = false;
	switch (m_mType)
	{
	case NORMAL:
		if (m_pProgMesh != NULL)
		{
			*validF = m_pProgMesh->numVisTris();
			*validV = 0;
			succ = true;
		}
		break;
	case QSLIM:
		if (m_pQslimMesh != NULL)
		{
			*validF = slim->valid_faces;
			*validV = slim->valid_verts;
			succ = true;
		}
		break;
	case VIEW_SELECT_OUR:
		if (m_viewMyRenderOur != NULL)
		{
			*validF = m_viewMyRenderOur->mesh3d->faces.size();
			*validV = m_viewMyRenderOur->mesh3d->vertices.size();
			succ = true;
		}
	}
	return succ;
}
// 得到当前模型的包围球半径
float CMeshSimpDoc::GetModelRadius()
{
	if (m_pQslimMesh)
		return sqrt(m_miniBall->squared_radius());
	else if (m_viewMyRenderOur)
		return m_viewMyRenderOur->mesh3d->bsphere.r;

	return 0.0f;
}
// 对当前模型进行单位化
void CMeshSimpDoc::OnNormalizeSize()
{
	if (m_pProgMesh) m_pProgMesh->Normalize();// center mesh around the origin & shrink to fit
}
// 近裁剪面的人眼观察极限尺寸
double CMeshSimpDoc::GetFilterRadiusInNearclip()
{
	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 
	float d, r, l; 
	pView->GetGLViewParam(&d, &r, &l);

	double fpm = pView->Get1PixelInMeter();

	switch (m_filterRadius)
	{
	case PIXEL1:
		break;
	case PIXEL2:
		fpm *= 2.0;
		break;
	case CSF:
	//////////////////////////////////////////////////////////////////////////
	// 使用严格的CSF，人类视觉极限是60c/deg，即1°最多120个黑白条 [5/23/2012 Han]
	// 故屏幕空间误差最多是1/120cm, 结果为: ppi / (120 * 2.54)
	//////////////////////////////////////////////////////////////////////////
		fpm *= r / (120 * 2.54);
		break;
	case CSF2:
		fpm *= 2.0 * r / (120 * 2.54);
		break;
	case BADMETRIC:
		fpm *= 10;
		break;
	}
	return fpm;
}
// 利用距离来调整lod级别
void CMeshSimpDoc::AdjustLOD(float fDist)
{
	double t = 0.0;
	int nOFace, nEface;
	if (m_SimpByMethod == NONE)
		return;
	fprintf(outputFile, "--Adjust mesh LOD based on view distance:\n" );

	switch(m_SimpByMethod)
	{
	case SIMP_BY_DIST:
		switch (m_mType)
		{
		case NORMAL:
			if (m_pProgMesh)
			{
				m_pProgMesh->AdjustLOD(fDist);
			}
			break;
		case QSLIM:
			if (m_pQslimMesh)
			{
				CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
				CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 
				float d, r, l; 
				pView->GetGLViewParam(&d, &r, &l);
				//double fpm = pView->Get1PixelInMeter();
				// 类型不同的话得到的近裁剪面长度不一样
				double fpm = GetFilterRadiusInNearclip();

				nOFace = slim->valid_faces;
				double err;
				//////////////////////////////////////////////////////////////////////////
				// 使用严格的CSF，人类视觉极限是60c/deg，即1°最多120个黑白条 [5/23/2012 Han]
				// 故屏幕空间误差最多是1/120cm, 结果为: ppi / (120 * 2.54)
				//fpm *= r / (120 * 2.54);
				//////////////////////////////////////////////////////////////////////////
				TIMING(t, err = slim->AdjustLOD(d, l, fpm));
				nEface = slim->valid_faces;
				fprintf(outputFile, "View distance:\t%f\t; Model space error: \t%lf\t;Decmiate from\t %d\t to \t%d\t; Using time:\t%lf seconds\n" 
					,fDist, err, nOFace, nEface, t);

				//fprintf(outputFile, "测试：vert:%d,face:%d\n", m_pQslimMesh->vert_count(), m_pQslimMesh->face_count());

				// put buffer data to disk [6/8/2012 Han]
				fflush(outputFile);

			}
			break;
		}
		break;
	case SIMP_BY_MSDM:
		switch (m_mType)
		{
		case QSLIM:
			if (m_pQslimMesh && m_pOriginalQMesh)
			{
				MSDM2::MSDM2_computation(m_pOriginalQMesh, m_pQslimMesh);
			}
			break;
		}
		break;
	}
}
// 得到当前视点下的可见面片,并赋予不同的颜色,以便观察
void CMeshSimpDoc::OnEditVisiblefaces()
{
	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 
	glDisable(GL_CULL_FACE);
	pView->Invalidate();
	pView->DisplayModel();


	int visibleFN = 0;
	bool *pVisF = m_pQslimMesh->GetVisibleFaces(visibleFN);

	glEnable(GL_CULL_FACE);

	for (int f = 0; f < m_pQslimMesh->face_count(); f++)
	{
		if (!pVisF[f])
		{
			m_pQslimMesh->face(f).colorIndex = 10;
		}
		else
			m_pQslimMesh->face(f).colorIndex = 200;
	}
	pView->bSmooth_ = false;
	pView->showing_type = CMeshSimpView::SHOWING_MESH_SALIENCY;
	pView->Invalidate();
	delete []pVisF;
	return;
}
void CMeshSimpDoc::OnEditVisiblevertexes()
{
	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 
	glDisable(GL_CULL_FACE);
	pView->Invalidate();
	pView->DisplayModel();

	int visibleVN = 0;
	bool *pVisF = m_pQslimMesh->GetVisibleVerts(visibleVN);

	glEnable(GL_CULL_FACE);

	double R = 0.0;
	int indiceLut = 0;

	for (int vID = 0; vID < m_pQslimMesh->vert_count(); vID++)
	{
		if (!pVisF[vID])
			indiceLut = 10;
		else
			indiceLut = 200;
		m_pQslimMesh->vertex(vID).colorIndex = indiceLut;
	}
	pView->bSmooth_ = true;
	pView->showing_type = CMeshSimpView::SHOWING_MESH_SALIENCY;
	pView->Invalidate();
	delete []pVisF;
	return;
}
// ↑是得到或者修改模型信息的相关函数
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// ↓和绘制相关的函数 [8/13/2013 Han]
// 绘制候选视点所在球面
void CMeshSimpDoc::DrawViewSphere(float radius)
{
	glPushMatrix();
	glPushAttrib(GL_COLOR_BUFFER_BIT|GL_POINT_BIT|GL_POLYGON_BIT);
	//////////////////////////////////////////////////////////////////////////
	// test [5/15/2012 Han]
	glEnable(GL_CULL_FACE); // backface culling
	glPointSize(8);
	glDisable(GL_LIGHTING);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glScalef(radius,radius, radius);

	glColor4f(0.2f, 0.6, 0.1, 0.3);

	// 只绘制候选视点 [10/9/2013 Han]
	/*	
	glBegin(GL_POINTS);
	vector<ViewPoint*>::iterator iter;
	for (iter = viewpoint_candidates.begin(); iter != viewpoint_candidates.end(); ++iter)
	{
		glColor3fv((*iter)->color);
		glVertex3fv((*iter)->pos);
	}
	glEnd();
	*/

	glBegin(GL_TRIANGLES);
	for(int i = 0; i < vf.size(); i++)
	{
		glColor4f(viewpoint_candidates.at(vf[i].p1)->color[0],viewpoint_candidates.at(vf[i].p1)->color[1],viewpoint_candidates.at(vf[i].p1)->color[2],0.5f);
		glVertex3fv(viewpoint_candidates.at(vf[i].p1)->pos);
		glColor4f(viewpoint_candidates.at(vf[i].p2)->color[0],viewpoint_candidates.at(vf[i].p2)->color[1],viewpoint_candidates.at(vf[i].p2)->color[2],0.5f);
		glVertex3fv(viewpoint_candidates.at(vf[i].p2)->pos);
		glColor4f(viewpoint_candidates.at(vf[i].p3)->color[0],viewpoint_candidates.at(vf[i].p3)->color[1],viewpoint_candidates.at(vf[i].p3)->color[2],0.5f);
		glVertex3fv(viewpoint_candidates.at(vf[i].p3)->pos);
	}
	glEnd();

	glPolygonMode( GL_FRONT , GL_LINE);
	glBegin(GL_TRIANGLES);
	for(int i = 0; i < vf.size(); i++)
	{
		glColor3fv(viewpoint_candidates.at(vf[i].p1)->color);
		glVertex3fv(viewpoint_candidates.at(vf[i].p1)->pos);
		glColor3fv(viewpoint_candidates.at(vf[i].p2)->color);
		glVertex3fv(viewpoint_candidates.at(vf[i].p2)->pos);
		glColor3fv(viewpoint_candidates.at(vf[i].p3)->color);
		glVertex3fv(viewpoint_candidates.at(vf[i].p3)->pos);
	}
	glEnd();

	glPolygonMode( GL_FRONT , GL_POINT);
	glColor4f(1.0f, 0.0, 0.0, 0.3);
	glBegin(GL_TRIANGLES);
	for(int i = 0; i < vf.size(); i++)
	{
		glColor4f(viewpoint_candidates.at(vf[i].p1)->color[0],viewpoint_candidates.at(vf[i].p1)->color[1],viewpoint_candidates.at(vf[i].p1)->color[2],0.5f);
		glVertex3fv(viewpoint_candidates.at(vf[i].p1)->pos);
		glColor4f(viewpoint_candidates.at(vf[i].p2)->color[0],viewpoint_candidates.at(vf[i].p2)->color[1],viewpoint_candidates.at(vf[i].p2)->color[2],0.5f);
		glVertex3fv(viewpoint_candidates.at(vf[i].p2)->pos);
		glColor4f(viewpoint_candidates.at(vf[i].p3)->color[0],viewpoint_candidates.at(vf[i].p3)->color[1],viewpoint_candidates.at(vf[i].p3)->color[2],0.5f);
		glVertex3fv(viewpoint_candidates.at(vf[i].p3)->pos);
	}
	glEnd();

	glScalef(1.0f,1.0f,1.0f);
	glPopAttrib();
	glPopMatrix();
}
// 判断当前顶点是否已经在视点集当中，如果是，则返回序号，否则返回-1 [8/21/2012 Han]
int PointIndex(vector<ViewPoint*>& pointSet, const point* p)
{
	vector<ViewPoint*>::iterator iter;
	int idx = 0;
	if (pointSet.size() > 0)
	{
		for (iter = pointSet.begin(); iter != pointSet.end(); ++iter, ++idx)
		{
			if ((*iter)->pos == *p)
			{
				return idx;
			}
		}
	}
	return -1;
}
/*
http://paulbourke.net/miscellaneous/sphere_cylinder/
   Create a triangular facet approximation to a sphere
   Return the number of facets created.
   The number of facets will be (4^iterations) * 8
*/
// 对球进行平均划分，划分后的顶点作为视点 [5/15/2012 Han]
void CMeshSimpDoc::EqualDistributeViews(float r, int itera)
{
	// 如果已经有平均分布的视点，则直接调整其半径即可
	vector<ViewPoint*>::iterator iter;
	if (viewpoint_candidates.size() > 0)
	{
		for (iter = viewpoint_candidates.begin(); iter != viewpoint_candidates.end(); ++iter)
			Normalise(&((*iter)->pos), r);
		return;
	}

		//////////////////////////////////////////////////////////////////////////
		// 采用球面平均分布点的方式进行视点初始化
	int iterations = itera;
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
			unsigned int paI, pbI, pcI;

			paI = PointIndex(viewpoint_candidates, &pa);
			if(paI == -1) 
			{
				paI = viewpoint_candidates.size();
				ViewPoint *pV = new ViewPoint;
				pV->pos = pa;
				viewpoint_candidates.push_back(pV);
			}
			
			pbI = PointIndex(viewpoint_candidates, &pb);
			if (pbI == -1)
			{
				pbI = viewpoint_candidates.size();
				ViewPoint *pV = new ViewPoint;
				pV->pos = pb;
				viewpoint_candidates.push_back(pV);
			}

			pcI = PointIndex(viewpoint_candidates, &pc);
			if (pcI == -1)
			{
				pcI = viewpoint_candidates.size();
				ViewPoint *pV = new ViewPoint;
				pV->pos = pc;
				viewpoint_candidates.push_back(pV);
			}
//每次循环得到的第一个视点和前一次循环的最后一个视点是同样的,这是由原算法造成，暂时不修改			
			vf[nt].p1 = vf[i].p1; vf[nt].p2 = paI/*pa*/; vf[nt].p3 = pcI/*pc*/; nt++;
			vf[nt].p1 = paI/*pa*/; vf[nt].p2 = vf[i].p2; vf[nt].p3 = pbI/*pb*/; nt++;
			vf[nt].p1 = pbI/*pb*/; vf[nt].p2 = vf[i].p3; vf[nt].p3 = pcI/*pc*/; nt++;
			vf[i].p1 = paI/*pa*/;
			vf[i].p2 = pbI/*pb*/;
			vf[i].p3 = pcI/*pc*/;
		}
	}

	// 载入观察球模型，用于clustering计算等 [10/9/2013 Han]
	slim_cleanup();
	// Process command line and read input model(s)
	startup_and_input("ViewSphere.smf");
	m_pViewPhereMdl = m_pQslimMesh->clone();
	slim_cleanup();

	// test [10/9/2013 Han] 经测试，模型中的顶点顺序和生成时顺序一致
	//for (int i = 0; i < m_pViewPhereMdl->vert_count(); i++)
	//{
	//	fprintf(outputFile, "v %f %f %f\n" , m_pViewPhereMdl->vertex(i)[0], m_pViewPhereMdl->vertex(i)[1], m_pViewPhereMdl->vertex(i)[2]);
	//}

	// 将只做好的观察球保存为smf文件格式的3D模型，以便进行聚类等工作 [10/9/2013 Han]
	// 如需要保存，取消注释即可
/*		FILE *tf = fopen("SphereModel.smf", "wt");
	for (iter = viewpoint_candidates.begin(); iter != viewpoint_candidates.end(); ++iter)
	{
		fprintf(tf, "v %f %f %f\n" , (*iter)->pos[0], (*iter)->pos[1],(*iter)->pos[2]);
	}

	for(int i = 0; i < vf.size(); i++)
	{
		fprintf(tf, "f %d %d %d\n" , vf[i].p1+1, vf[i].p2+1,vf[i].p3+1);

	}
	fclose(tf);
*/
}
// 按照候选视点的得分赋予候选视点颜色
void CMeshSimpDoc::UpdateSphereColorByImportance()
{
	vector<ViewPoint*>::iterator iter;
	float fmax = 0, fmin = FLT_MAX;

	for (iter = viewpoint_candidates.begin(); iter != viewpoint_candidates.end(); ++iter)
	{
		if ((*iter)->importance < fmin)
			fmin = (*iter)->importance;
		if ((*iter)->importance > fmax)
			fmax = (*iter)->importance;
	}

	double R = 0.0;
	int indiceLut = 0;
	double len = 1.0/(fmax - fmin);

	for (iter = viewpoint_candidates.begin(); iter != viewpoint_candidates.end(); ++iter)
	{
		R = ((*iter)->importance - fmin) * len * 255;

		if(R>255)
			R=255;
		if (R < 0)
			R=0;

		indiceLut=floor(R);

		(*iter)->color[0] = m_pQslimMesh->LUT_CourbureClust[3*indiceLut];
		(*iter)->color[1] = m_pQslimMesh->LUT_CourbureClust[3*indiceLut+1];
		(*iter)->color[2] = m_pQslimMesh->LUT_CourbureClust[3*indiceLut+2];
	}
}
// 输出所有的候选视点渲染的图片
void CMeshSimpDoc::OnShadingOutputallcandviews()
{
	void SaveRenderPic(const point *viewpos,float d, wchar_t* fn, CMeshSimpView *pView);
	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 
	pView->resetOrientation();
	float dv = pView->GetViewDistance();
	MxStdModel *p = m_pQslimMesh;
	m_pQslimMesh = m_pOriginalQMesh;
	wchar_t fn[128];
	vector<ViewPoint*>::iterator iter;
	int idx = 0;
	if (viewpoint_candidates.size() > 0)
	{
		for (iter = viewpoint_candidates.begin(); iter != viewpoint_candidates.end(); ++iter, ++idx)
		{
			swprintf( fn,   L"%d_%f_%f_%f_CandView.ppm ", idx, (*iter)->pos[0], (*iter)->pos[1], (*iter)->pos[2]);
			SaveRenderPic(&((*iter)->pos), dv, fn, pView);
		}
	}
	m_pQslimMesh = p;
}
// ↑和绘制相关的函数 [8/13/2013 Han]
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// ↓ 用于临时测试的函数,在edit菜单命令下
void CMeshSimpDoc::OnEditTestonly2()
{
	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 

	// 测试：利用不同尺度计算平均曲率，然后求平均，得到多尺度平均曲率之和 [11/6/2012 Han]
	double tS; 
	if (m_pQslimMesh != NULL)
	{
		float *tempC = new float[m_pQslimMesh->vert_count()];
		memset(tempC, 0.f, m_pQslimMesh->vert_count()*sizeof(float));
		// 暂时不进行直方图处理 [7/14/2012 Han]
		//float size[3];
		//m_pQslimMesh->GetBoundSize(size);
		//float max = 0;
		//for (int i = 0; i < 3; i ++)
		//	max = max < size[i]? size[i] : max;

		//double RadiusCurvature = 0.0;
		//float sss = 0.0f;
		//for (int i = 0; i < 3; i++)
		//{
		//	double RadiusCurvature=/*Play with this value*//*3**/MSDM2::mini_radius+i * 3 * MSDM2::radius_step;
		//	TIMING(tS,m_pQslimMesh->principal_curvature(true,RadiusCurvature*max, 0));
		//	fprintf(outputFile, "Radius:%lf,Face Num: %d，Calc Curvature time: %lf seconds\n", RadiusCurvature, m_pQslimMesh->face_count(), tS);
		//	MSDM2::KmaxKmean(m_pQslimMesh,max, 0);		// 求平均曲率

		//	for (int i = 0; i < m_pQslimMesh->vert_count(); i++)
		//	{
		//		sss = m_pQslimMesh->vertex(i).KmaxCurv[0];
		//		tempC[i] += sss;
		//	}
		//}

		// 再次计算一个较高尺度radius的curvature， 然后和前面小尺度curvature求平均 [11/6/2012 Han]
		//RadiusCurvature=/*Play with this value*/MSDM2::mini_radius+8 * MSDM2::radius_step;
		//TIMING(tS,m_pQslimMesh->principal_curvature(true,RadiusCurvature*max, 0));
		//fprintf(outputFile, "Radius:%lf,Face Num: %d，Calc Curvature time: %lf seconds\n", RadiusCurvature, m_pQslimMesh->face_count(), tS);
		//MSDM2::KmaxKmean(m_pQslimMesh,max, 0);		// 求平均曲率


		// copy mean curvature to importance [9/5/2012 Han]
		//m_pQslimMesh->MinImportance = FLT_MAX, m_pQslimMesh->MaxImportance = 0;
		//for (int i = 0; i < m_pQslimMesh->vert_count(); i++)
		//{
		//	sss = m_pQslimMesh->vertex(i).KmaxCurv[0] + tempC[i];
		//	m_pQslimMesh->vertex(i).view_importance = sss;
		//	if (sss > m_pQslimMesh->MaxImportance)
		//		m_pQslimMesh->MaxImportance = sss;
		//	if (sss < m_pQslimMesh->MinImportance)
		//		m_pQslimMesh->MinImportance = sss;
		//}
		// 输出为临时文件，借以统计直方图 [3/7/2013 Han]
		char fnt[128];
		//sprintf( fnt, "Histogram%d.txt", rand());
		//FILE *outputHistoFile = fopen(fnt, "wt");
		//fprintf(outputHistoFile, "Original Histogram,Ver max:%f, Ver min:%f\n", m_pQslimMesh->MaxImportanceV,m_pQslimMesh->MinImportanceV);

		sprintf( fnt, "HistogramArea%d.txt", rand());
		FILE *outputHistoArea = fopen(fnt, "wt");
		fprintf(outputHistoArea, "Original Histogram Area, Face max:%f,Face min:%f\n", m_pQslimMesh->MaxImportanceF,m_pQslimMesh->MinImportanceF);

		std::vector<float> histogram_before;
		histogram_before.resize(256,0.f);
		std::vector<float> histogram_after;
		histogram_after.resize(256,0.f);
		std::vector<float> histogram_after32;
		histogram_after32.resize(32,0.f);


		// 对结果进行smooth，并进行离散化 [11/6/2012 Han]
		//OnEditSmoothsaliency();
		// 对mean curvature 进行离散化 [9/18/2012 Han]
		OnViewselectionMeshcurvatureGeo();
		float hs = (m_pQslimMesh->MaxImportanceF - m_pQslimMesh->MinImportanceF)/256;
		for (unsigned int i=0; i<m_pQslimMesh->face_count(); i++)
		{
			float area = m_pQslimMesh->compute_face_area(i);
			unsigned int idx = int((m_pQslimMesh->face(i).view_importance - m_pQslimMesh->MinImportanceF)/hs);
			idx = idx >= 256? (256-1):idx;
			// 将面积计入原始模型面积统计当中，即q [12/29/2012 Han]
			histogram_before[idx] += area;
		}
		// 进行直方图均衡化 [3/7/2013 Han]
		m_pQslimMesh->Histoeq(256);
		for (unsigned int i=0; i<m_pQslimMesh->face_count(); i++)
		{
			float area = m_pQslimMesh->compute_face_area(i);
			unsigned int idx = int((m_pQslimMesh->face(i).view_importance - m_pQslimMesh->MinImportanceF)/hs);
			idx = idx >= 256? (256-1):idx;
			// 将面积计入原始模型面积统计当中，即q [12/29/2012 Han]
			histogram_after[idx] += area;
		}
		m_pQslimMesh->VertexColor(3);

		hs = (m_pQslimMesh->MaxImportanceF - m_pQslimMesh->MinImportanceF)/32; //histogram step width
		for (unsigned int i=0; i<m_pQslimMesh->face_count(); i++)
		{
			float area = m_pQslimMesh->compute_face_area(i);
			unsigned int idx = int((m_pQslimMesh->face(i).view_importance - m_pQslimMesh->MinImportanceF)/hs);
			idx = idx >= 32? (32-1):idx;
			// 将面积计入原始模型面积统计当中，即q [12/29/2012 Han]
			histogram_after32[idx] += area;
		}
		for (vector<float>::size_type i=0; i<256; i++)
		{
			fprintf(outputHistoArea, "%f\t", histogram_before[i]);
			fprintf(outputHistoArea, "%f\t", histogram_after[i]);
			if (i < 32)
				fprintf(outputHistoArea, "%f\t", histogram_after32[i]);
			fprintf(outputHistoArea, "\n");
		}
		//fclose(outputHistoFile);
		fclose(outputHistoArea);

		// 将经过均衡化的曲率离散化为32阶 [5/10/2013 Han]
		for (unsigned int i=0; i<m_pQslimMesh->face_count(); i++)
		{
			float area = m_pQslimMesh->compute_face_area(i);
			unsigned int idx = int((m_pQslimMesh->face(i).view_importance - m_pQslimMesh->MinImportanceF)/hs);
			idx = idx >= 32? (32-1):idx;
			// 将面积计入原始模型面积统计当中，即q [12/29/2012 Han]
			m_pQslimMesh->face(i).view_importance =  m_pQslimMesh->MinImportanceF + idx * hs;
		}
		// 对saliency进行规整化 [7/15/2012 Han]
		m_pQslimMesh->IdentityVertexImportance(m_pQslimMesh->MinImportanceV, m_pQslimMesh->MaxImportanceV);
		//m_pQslimMesh->UpdateFaceImportanceByVert();

		m_pQslimMesh->VertexColor(3/*, m_pQslimMesh->MaxImportance, m_pQslimMesh->MinImportance*/);
		delete []tempC;

	}
	pView->showing_type = pView->SHOWING_VIEW_DEPENDENT_CURVATURE;
	pView->Invalidate();
	return;
	// 测试结束 [11/6/2012 Han]


	//////////////////////////////////////////////////////////////////////////
	// 将视角重置

	//float r = GetModelRadius();
	//if ( r > 0)
	//{
	//	// 屏幕空间和三维空间的对应关系:
	//	// glFrstrom中,假定eye位置为原地,最近裁剪面表示投影平面里原点的距离,投影四边形的4个顶点位置就确定了投影平面的大小,以及Fov角度
	//	// gluLookAt函数所指定的eye位置即调整后的视点到原点的距离,其实相当于做相反方向的glTranslate函数
	//	// 故而,可以利用以下公式得到当模型包围球半径为r时,模型离视点的距离为多少时,模型投影能占据整个屏幕
	//	r /= sin(M_PI*(45/2.0)/180);
	//	pView->eye.z = r;
	//}
	//pView->AdjustView();	
	//pView->Invalidate();


	//////////////////////////////////////////////////////////////////////////
	// 首先进行saliency，然后smooth，接着使用saliency进行分割
	float dist, r, l; 
	pView->GetGLViewParam(&dist, &r, &l);
	//double fpm = pView->Get1PixelInMeter();
	// 类型不同的话得到的近裁剪面长度不一样
	double fpm = GetFilterRadiusInNearclip();

	double fUnprojRadius = fpm*dist/l;

	// using mesh saliency [7/10/2012 Han]
	//OnViewselectionMeshsaliency();
	// Using Curvature, cause it efficiency [7/10/2012 Han]
	float size[3];
	m_pQslimMesh->GetBoundSize(size);
	float max = 0;
	for (int i = 0; i < 3; i ++)
		max = max < size[i]? size[i] : max;

	m_pQslimMesh->principal_curvature(true,fUnprojRadius, 0);
	MSDM2::KmaxKmean(m_pQslimMesh,max, 0);
	for(int i = 0; i < m_pQslimMesh->vert_count(); i++)
	{
		m_pQslimMesh->vertex(i).view_importance = m_pQslimMesh->vertex(i).KmaxCurv[0];
	}
	m_pQslimMesh->UpdateFaceImportanceByVert();


	//m_pQslimMesh->LaplacianFilterVertexImportance(true, fUnprojRadius, NULL);
	// test value adjust [7/10/2012 Han]
	//int n;
	//segment_image(m_pQslimMesh,  sigma,k , min , &n, true);
	//m_pQslimMesh->UpdateSaliencyBySegment();
	//m_pQslimMesh->VertexColor(3);
	//pView->showing_type = pView->SHOWING_MESH_SALIENCY;
	//pView->Invalidate();
	OnEditSmoothsaliency();
	OnEditSegbysaliency();
}
// Test something here
void CMeshSimpDoc::OnEditTestonly()
{
	// 测试对顶点曲率进行直方图均衡化 [3/7/2013 Han]
	OnViewselectionMeshcurvatureGeo();
	m_pQslimMesh->Histoeq(256);
	m_pQslimMesh->IdentityVertexImportance();
	//m_pQslimMesh->UpdateFaceImportanceByVert();
	m_pQslimMesh->VertexColor(3);
	return;
	// 测试顶点曲率直方图均衡化结束 [3/7/2013 Han]



	// 测试某个点的周围半径顶点是否能够正确选取 [6/5/2012 Han]
	//if (m_pQslimMesh != NULL)
	//{
	//	float r = GetModelRadius();
	//	m_pQslimMesh->LaplacianFilterVertexImportance(true, r * 0.2, 0, NULL);
	//	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	//	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 
	//	pView->showing_type = pView->SHOWING_MESH_SALIENCY;
	//	pView->Invalidate();
	//}
	// 测试对每个顶点的saliency进行滤波 [6/5/2012 Han]
	// 需要首先计算完成saliency [6/5/2012 Han]
	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 

	if (m_pQslimMesh != NULL)
	{
		// Get some geometry values [5/24/2012 Han]
		float dist, r, l; 
		pView->GetGLViewParam(&dist, &r, &l);
		//double fpm = pView->Get1PixelInMeter();
		// 类型不同的话得到的近裁剪面长度不一样
		double fpm = GetFilterRadiusInNearclip();


		float fMR = GetModelRadius();

		double fUnprojRadius = fpm*dist/l;

		if (m_viewMyRenderOur == NULL )
		{
			OnViewselectionMeshsaliency();
			for(int i = 0; i < m_pQslimMesh->vert_count(); i++)
				m_pOriginalQMesh->vertex(i).view_importance = m_pQslimMesh->vertex(i).view_importance;

			pView->showing_type = pView->SHOWING_MESH_SALIENCY;
			pView->Invalidate();

			return;
		}

		for(int i = 0; i < m_pQslimMesh->vert_count(); i++)
			m_pQslimMesh->vertex(i).view_importance = m_pOriginalQMesh->vertex(i).view_importance;
		

		//vector<ViewPoint*>::iterator iter;
		//for (iter = viewpoint_candidates.begin(); iter != viewpoint_candidates.end(); ++iter)
		//{
		//	pView->eye.x = (*iter)->pos[0];
		//	pView->eye.y = (*iter)->pos[1];
		//	pView->eye.z = (*iter)->pos[2];
		//	pView->eye = pView->eye*dist;
		//	pView->AdjustView();
		//}
		// save saliency mapped mesh pic [6/6/2012 Han]

		m_pQslimMesh->LaplacianFilterVertexImportance(true, fUnprojRadius, outputFile);
		m_pQslimMesh->VertexColor(3, m_pQslimMesh->MaxImportanceV, m_pQslimMesh->MinImportanceV);
		pView->showing_type = pView->SHOWING_MESH_SALIENCY;
		pView->Invalidate();
	}
}
// ↑ 用于临时测试的函数,在edit菜单命令下
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// ↓ 计算和操作saliency值
// 备份一份顶点重要度信息
void CopyImportances(MxStdModel* pSrc, MxStdModel* pDest)
{
	// 原始模型备份一份最初的saliency值 [7/10/2012 Han]
	for (int i = 0; i < pSrc->vert_count(); i++)
		pDest->vertex(i).view_importance = pSrc->vertex(i).view_importance;
	for (int i = 0; i < pSrc->face_count(); i++)
		pDest->face(i).view_importance = pSrc->face(i).view_importance;

}
// 计算模型saliency值,如果用户选择载入文件的方式,则直接将文件中的saliency值读取
void CMeshSimpDoc::OnEditCalcsaliency()
{
	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 

	// test delete [7/14/2012 Han]
	// 测试使用心得利用测地线平均曲率的方法计算saliency [7/14/2012 Han]
	//m_pQslimMesh->compute_mesh_saliency(GetModelRadius(),NULL, true);

	//m_pQslimMesh->VertexColor(3);
	//pView->showing_type = pView->SHOWING_MESH_SALIENCY;
	//pView->Invalidate();
	//return;


	// 直接用文件载入，不进行计算 [7/12/2012 Han]
	// load saliency but calc it
	if (!LoadMeshSaliency(NULL))
		/*return;*/
		OnViewselectionMeshsaliency();

	// 计算saliency [7/12/2012 Han]
	//OnViewselectionMeshsaliency();
	// 原始模型备份一份最初的saliency值 [7/10/2012 Han]
	CopyImportances(m_pQslimMesh, m_pOriginalQMesh);

	pView->showing_type = pView->SHOWING_MESH_SALIENCY;
	pView->Invalidate();
}
// 对计算得到的模型saliency值进行平滑,平滑半径是csf计算结果
void CMeshSimpDoc::OnEditSmoothsaliency()
{
	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 

	float dist, r, l; 
	pView->GetGLViewParam(&dist, &r, &l);
	//double fpm = pView->Get1PixelInMeter();
	// 类型不同的话得到的近裁剪面长度不一样
	double fpm = GetFilterRadiusInNearclip();

	double fUnprojRadius = fpm*dist/l;

	m_pQslimMesh->LaplacianFilterVertexImportance(true, fUnprojRadius, outputFile);


	m_pQslimMesh->VertexColor(3);
	pView->showing_type = pView->SHOWING_MESH_SALIENCY;
	pView->Invalidate();
	// 原始模型备份一份最初的saliency值 [7/10/2012 Han]
	//CopyImportances(m_pQslimMesh, m_pOriginalQMesh);
}
// 利用saliency值对模型表面进行聚类分解,得到saliency atlas
void CMeshSimpDoc::OnEditSegbysaliency()
{
	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 

	// 用户输入分割参数 [8/15/2012 Han]
	//SegmentParam dlg;
	//dlg.SetValues(sigma, k, min_size,segFace);
	//if (dlg.DoModal() == IDOK)
	//{
	//	dlg.GetValues(sigma, k, min_size,segFace);
	//}
	//else
	//	return;

	float maxOrig = segFace?m_pQslimMesh->MaxImportanceF:m_pQslimMesh->MaxImportanceV;
	float minOrig = segFace?m_pQslimMesh->MinImportanceF:m_pQslimMesh->MinImportanceV;

	fprintf(outputFile, "--Segmentation by Mesh Saliency:\nMax Saliency:\t%f\t;Min Saliency:\t%f\tAfter Smooth\n", maxOrig, minOrig);

	int n;
	// 首先删除已有数据 [7/12/2012 Han]
	m_pQslimMesh->saliencySegData.clear();

	double ts =0.0;
	//////////////////////////////////////////////////////////////////////////
	// 计算各个控制参数的依据，k：最大saliency和最小saliency的1/255, min: 模型面片数或者顶点数的0.0004, [7/25/2012 Han]
	//////////////////////////////////////////////////////////////////////////
	//TIMING(ts, segment_image(m_pQslimMesh,  sigma,k , min , &n, segFace/*true*/));
	TIMING(ts, segment_image(m_pQslimMesh,  sigma
		, segFace?k*(m_pQslimMesh->MaxImportanceF - m_pQslimMesh->MinImportanceF):k*(m_pQslimMesh->MaxImportanceV - m_pQslimMesh->MinImportanceV) 
		, min_size*(segFace?slim->valid_faces:slim->valid_verts) , &n, segFace));

	m_pQslimMesh->UpdateSaliencyBySegment(segFace);

	m_pQslimMesh->VertexColor(3, maxOrig, minOrig);
	pView->showing_type = pView->SHOWING_MESH_SALIENCY;
	pView->Invalidate();

	maxOrig = segFace?m_pQslimMesh->MaxImportanceF:m_pQslimMesh->MaxImportanceV;
	minOrig = segFace?m_pQslimMesh->MinImportanceF:m_pQslimMesh->MinImportanceV;

	fprintf(outputFile, "Max Saliency:\t%f\t;Min Saliency:\t%f\tAfter Segment\n", maxOrig, minOrig);

	fprintf(outputFile, "Segment Time:%lf, Segment Params:sigma=%f,k=%f(real:%f),min=%f(real:%f)---->Segment Number:%d-Real Num:%d\n"
		,ts, sigma, k,segFace?k*(m_pQslimMesh->MaxImportanceF - m_pQslimMesh->MinImportanceF):k*(m_pQslimMesh->MaxImportanceV - m_pQslimMesh->MinImportanceV) 
		, min_size, min_size*(segFace?slim->valid_faces:slim->valid_verts) ,n,m_pQslimMesh->saliencySegData.size());

	//map<int,struct Seg>::iterator it;
	//for ( it=m_pQslimMesh->saliencySegData.begin() ; it != m_pQslimMesh->saliencySegData.end(); it++ )
	//{
	//	fprintf(outputFile, "SetID:%d, num:%d, saliency:%f\n", 		it->first, it->second.num, it->second.Saliency);
	//}

	fflush(outputFile);

	// 将saliency值进行还原 [7/10/2012 Han]
	//CopyImportances(m_pOriginalQMesh, m_pQslimMesh);
	//for (int i = 0; i < m_pQslimMesh->vert_count(); i++)
	//	m_pQslimMesh->vertex(i).view_importance = m_pOriginalQMesh->vertex(i).view_importance;

	// test delete [7/25/2012 Han]
	//char pszInfo[128];
	//sprintf(pszInfo, "Seg num:%d,Real Num:%d", n, m_pQslimMesh->saliencySegData.size());
	//MessageBox(NULL, LPCSTR(pszInfo), "Info", MB_OK);
}
// 载入已经保存的模型saliency文件
bool CMeshSimpDoc::LoadMeshSaliency(char *fn)
{
	CString FilePathName;
	if (fn == NULL)
	{
		// 不载入，直接进行saliency计算 [10/12/2013 Han]
		//CFileDialog dlg(TRUE);///TRUE为OPEN对话框，FALSE为SAVE AS对话框
		//if(dlg.DoModal()==IDOK)
		//	FilePathName=dlg.GetPathName();
		//else
		{
			OnViewselectionMeshsaliency();
			return true;
		}
		fn = FilePathName.GetBuffer(FilePathName.GetLength());
	}

	FILE *inputSaliencyFile = fopen(fn, "rt");
	int nf;
	int cv;
	double tS = 0.0;
	bool vd = true;
	switch(m_mType)
	{
	case QSLIM:
		// copy mesh saliency to qslim model [5/24/2012 Han]
		if (m_pQslimMesh != NULL)
		{

			m_pQslimMesh->MinImportanceV = FLT_MAX, m_pQslimMesh->MaxImportanceV = 0;

			fscanf(inputSaliencyFile, "%d", &nf);
			//if (nf < m_pQslimMesh->face_count())
			//{
			//	fclose(inputSaliencyFile);
			//	return false;
			//}
			// 如果是简化模型，则首先简化到指定的级别
			if (nf < slim->valid_faces)
			{
				int of = slim->valid_faces;
				double slim_time;
				TIMING(slim_time, slim->decimate(nf));
				fprintf(outputFile, "decimate from: %d faces to %d faces using:%lf seconds\n", of,slim->valid_faces, slim_time);
			}
			float sss = 0.0f;
			for (int i = 0; i < m_pQslimMesh->vert_count(); i++)
			{
				fscanf(inputSaliencyFile, "%d%f", &cv, &sss);
				m_pQslimMesh->vertex(i).view_importance = sss;
				if (sss > m_pQslimMesh->MaxImportanceV)
					m_pQslimMesh->MaxImportanceV = sss;
				if (sss < m_pQslimMesh->MinImportanceV)
					m_pQslimMesh->MinImportanceV = sss;
			}		
			//fscanf(inputSaliencyFile, "max saliency:\t%f\t;min saliency:\t%f\t\n", fmax, fmin);
			// 对saliency进行规整化 [7/15/2012 Han]
			m_pQslimMesh->IdentityVertexImportance(m_pQslimMesh->MinImportanceV, m_pQslimMesh->MaxImportanceV);
			m_pQslimMesh->UpdateFaceImportanceByVert();

			m_pQslimMesh->VertexColor(3, m_pQslimMesh->MaxImportanceV, m_pQslimMesh->MinImportanceV);
			fprintf(outputFile, "Mesh saliency loaded, Max saliency:%f\t; Min saliency:%f\n\n", m_pQslimMesh->MaxImportanceV, m_pQslimMesh->MinImportanceV);
		}
		break;
	case VIEW_SELECT_OUR:
	case VIEW_SELECT:
		// Not implement [7/12/2012 Han]
		break;
	}

	fclose(inputSaliencyFile);
	FilePathName.ReleaseBuffer();
	return true;
}
// ↑ 计算和操作saliency值
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// ↓为对菜单按钮的update函数 [8/13/2013 Han]
void CMeshSimpDoc::OnUpdateSimpMelax(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(m_pProgMesh != NULL&& m_mType == NORMAL);
	pCmdUI->SetCheck(m_edgemethod == PMesh::MELAX);
}
void CMeshSimpDoc::OnUpdateSimpQuadric(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(m_pProgMesh != NULL&& m_mType == NORMAL);
	pCmdUI->SetCheck(m_edgemethod == PMesh::QUADRIC);
}
void CMeshSimpDoc::OnUpdateSimpShorttest(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(m_pProgMesh != NULL&& m_mType == NORMAL);
	pCmdUI->SetCheck(m_edgemethod == PMesh::SHORTEST);
}
void CMeshSimpDoc::OnUpdateSimpQuadricWeighted(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(m_pProgMesh != NULL&& m_mType == NORMAL);
	pCmdUI->SetCheck(m_edgemethod == PMesh::QUADRICTRI);
}
void CMeshSimpDoc::OnUpdateSimplificationalgorithmSimpbydistance(CCmdUI *pCmdUI)
{
	//pCmdUI->SetCheck(m_bSimpByDist);
}
void CMeshSimpDoc::OnUpdateSimplificationalgorithmQuadricweightedbysegment(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_edgemethod == PMesh::QUADRIC_SEGMENT);
	pCmdUI->Enable(m_pProgMesh != NULL && m_mType == NORMAL);
	if (m_pProgMesh)
		pCmdUI->Enable(m_pProgMesh->HasSegment());
}
void CMeshSimpDoc::OnUpdateSimplificationalgorithmQuadricweightedbymotionblur(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_edgemethod == PMesh::BLUR);
	pCmdUI->Enable(m_pProgMesh != NULL && m_mType == NORMAL);
}
// 用于进行骨架抽取
void CMeshSimpDoc::OnUpdateSkelColl(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(m_pProgMesh != NULL && m_mType == NORMAL);
	if (m_pProgMesh != NULL && m_mType == NORMAL)
		pCmdUI->Enable(m_pProgMesh->bCollapsed);
}
// 用于进行模型收缩，进而得到收缩到近似骨架的状态
void CMeshSimpDoc::OnUpdateSkeletonizerGetskeleton(CCmdUI *pCmdUI)
{
	pCmdUI->Enable(m_pProgMesh != NULL&& m_mType == NORMAL);
	if (m_pProgMesh != NULL)
		pCmdUI->Enable(!m_pProgMesh->isCollapsing);
}
void CMeshSimpDoc::OnUpdateSimplificationalgorithmSimpbymsdm2(CCmdUI *pCmdUI)
{
	// TODO: Add your command update UI handler code here
}
// ↑为对菜单按钮的update函数 [8/13/2013 Han]
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// ↓为关于模型简化的函数 [8/13/2013 Han]
// 进行一次点分裂操作
void CMeshSimpDoc::OnResampleAdd1tri()
{
	if (m_pProgMesh)
	{
		CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
		CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 
		bool ret = m_pProgMesh->splitVertex();
		if (!ret) MessageBeep(0);
		pView->Invalidate();
	}
}
// 进行一次边折叠操作
void CMeshSimpDoc::OnResampleRemove1tri()
{
	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 
	double slim_time;
	switch (m_mType)
	{
	case NORMAL:
		if (m_pProgMesh)
		{
			bool ret = true;
			if (pView->m_bSelAny && !m_bSelCalc)
				if(MessageBox(pView->m_hWnd,LPCSTR("对选择的部分进行Remash？"),LPCSTR("Remash内容选择"), MB_YESNO) == IDYES)
				{
					delete m_pProgMesh;
					m_pProgMesh = new PMesh(m_pMesh, m_edgemethod, pView->m_selIndexList);
					m_bSelCalc = true;
				}

				//ret = m_pProgMesh->collapseEdge();
				ret = m_pProgMesh->collapseEdgeRealtime();	
				if (!ret) MessageBeep(0);
		}
		break;
	case QSLIM:
	case VIEW_SELECT:	// 二者都可以进行简化处理，处理完成后需要将状态设置为qslim
		if (m_pQslimMesh)
		{
			TIMING(slim_time, slim->decimate(slim->valid_faces-1));
			m_mType = QSLIM;
		}
		break;
	}
	pView->Invalidate();
}
// 细分5%
void CMeshSimpDoc::OnResampleIncreaseby5()
{
	if (m_pProgMesh)
	{
		int size = int(m_pProgMesh->numVertex()*0.05);
		if (size == 0) size = 1;
		bool ret = true;
		for (int i = 0; ret && i < size; ++i) {
			ret = m_pProgMesh->splitVertex();
		}
		if (!ret) MessageBeep(0);
		CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
		pMain-> GetActiveView()->Invalidate(); 
	}
}
// 简化5%
void CMeshSimpDoc::OnResampleReduceby5()
{
	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 
	double input_time, init_time, slim_time, output_time;

	switch (m_mType)
	{
	case NORMAL:

		if (m_pProgMesh)
		{
			int size = int(m_pProgMesh->numVertex()*0.05);
			if (size == 0) size = 1;
			bool ret = true;
			if (pView->m_bSelAny && !m_bSelCalc)
				if(MessageBox(pView->m_hWnd,LPCSTR("对选择的部分进行Remash？"),LPCSTR("Remash内容选择"), MB_YESNO) == IDYES)
				{
					//test
					size = 5;

					delete m_pProgMesh;
					m_pProgMesh = new PMesh(m_pMesh, m_edgemethod, pView->m_selIndexList);
					m_bSelCalc = true;
				}

				for (int i = 0; ret && i < size; ++i) {
					//ret = m_pProgMesh->collapseEdge();
					// test [9/20/2011 Han Honglei]
					ret = m_pProgMesh->collapseEdgeRealtime();
				}
				if (!ret) MessageBeep(0);	
		}
		break;
	case QSLIM:
	case VIEW_SELECT:	// 二者都可以进行简化处理，处理完成后需要将状态设置为qslim
		if (m_pQslimMesh)
		{
			int dF = slim->valid_faces - slim->valid_faces*0.7;
			TIMING(slim_time, slim->decimate(slim->valid_faces*0.7));
			m_mType = QSLIM;
			char pszInfo[_MAX_FNAME + 1];
			sprintf(pszInfo, "decimate %d faces using:%lf seconds", dF,slim_time);
			MessageBox(NULL, LPCSTR(pszInfo), "Info", MB_OK);
			fprintf(outputFile, "简化%d面片用时：%lf秒\n", dF, slim_time);
			// put buffer data to disk [6/8/2012 Han]
			fflush(outputFile);

		}
		break;
	}
	pView->Invalidate();
}
// 选择不同的简化方法,使用基于分块重要度的边折叠简化策略
void CMeshSimpDoc::OnSimplificationalgorithmQuadricweightedbysegment()
{
	if (m_edgemethod != PMesh::QUADRIC_SEGMENT)
	{
		m_edgemethod = PMesh::QUADRIC_SEGMENT;
		if (m_pProgMesh)
			m_pProgMesh->_cost = PMesh::QUADRIC_SEGMENT;
		changeSimplificationAlgorithm();
	}
}
// 选择Melax简化方法
void CMeshSimpDoc::OnSimpMelax()
{
	if (m_edgemethod != PMesh::MELAX)
	{
		m_edgemethod = PMesh::MELAX;
		if (m_pProgMesh)
			m_pProgMesh->_cost = PMesh::MELAX;
		changeSimplificationAlgorithm();
	}
}
// 使用加权qradric简化方法
void CMeshSimpDoc::OnSimpQuadricWeighted()
{
	if (m_edgemethod != PMesh::QUADRICTRI)
	{
		m_edgemethod = PMesh::QUADRICTRI;
		if (m_pProgMesh)
			m_pProgMesh->_cost = PMesh::QUADRICTRI;
		changeSimplificationAlgorithm();
	}
}
// 使用Quadric简化方法
void CMeshSimpDoc::OnSimpQuadric()
{
	if (m_edgemethod != PMesh::QUADRIC)
	{
		m_edgemethod = PMesh::QUADRIC;
		if (m_pProgMesh)
			m_pProgMesh->_cost = PMesh::QUADRIC;
		changeSimplificationAlgorithm();
	}

	// 使用Garland的算法进行模型简化 [11/24/2011 Han]
	//slim_init();
}
// 使用最短边简化方法
void CMeshSimpDoc::OnSimpShorttest()
{
	if (m_edgemethod != PMesh::SHORTEST)
	{
		m_edgemethod = PMesh::SHORTEST;
		if (m_pProgMesh)
			m_pProgMesh->_cost = PMesh::SHORTEST;
		changeSimplificationAlgorithm();
	}
}
// 使用基于motiong blur加权的边折叠简化方法
void CMeshSimpDoc::OnSimplificationalgorithmQuadricweightedbymotionblur()
{
	if (m_edgemethod != PMesh::BLUR)
	{
		m_edgemethod = PMesh::BLUR;
		if (m_pProgMesh)
			m_pProgMesh->_cost = PMesh::BLUR;
		changeSimplificationAlgorithm();
	}
}
// User has selected a new mesh simplification algorithm
void CMeshSimpDoc::changeSimplificationAlgorithm()
{
	if (PMesh::QUADRIC_SEGMENT == m_edgemethod && false == m_pProgMesh->HasSegment())
		return;

	// Need not to reload the mesh, cause the collapse list will be calced realtime [10/13/2011 Han Honglei]
	//if (m_pMesh == NULL)
	//{
	//	if (0 == strlen(m_sFilePath))
	//		return;

	//	m_pMesh = new jmsMesh(m_sFilePath);
	//	// 测试，不进行标准化，保持原来的大小
	//	//if (m_pMesh) m_pMesh->Normalize();
	//}
	//delete m_pProgMesh;
	//m_pProgMesh = new PMesh(m_pMesh, m_edgemethod);

	if (m_pProgMesh)
		m_pProgMesh->createEdgeCollapseList();
	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	pMain-> GetActiveView()->Invalidate(); 	
}
// 是否考虑视点观察距离
void CMeshSimpDoc::OnSimplificationalgorithmSimpbydistance()
{
	m_SimpByMethod = (m_SimpByMethod == SIMP_BY_DIST)? NONE:SIMP_BY_DIST;
}
void CMeshSimpDoc::OnFilterradius1pixel()
{
	m_filterRadius = PIXEL1;
}
void CMeshSimpDoc::OnFilterradius2pixel()
{
	m_filterRadius = PIXEL2;
}
void CMeshSimpDoc::OnFilterradiusCsf()
{
	m_filterRadius = CSF;
}
void CMeshSimpDoc::OnFilterradius2()
{
	m_filterRadius = CSF2;
}
// 使用基于msdm2误差控制的模型简化
void CMeshSimpDoc::OnSimplificationalgorithmSimpbymsdm2()
{
	// test delete [7/23/2012 Han]
	double es = MSDM2::MSDM2_computation(m_pOriginalQMesh, m_pQslimMesh, true);
	return;

	//////////////////////////////////////////////////////////////////////////
	m_SimpByMethod = (m_SimpByMethod == SIMP_BY_MSDM)? NONE:SIMP_BY_MSDM;
	if (m_SimpByMethod == NONE)
		return;

	slim->m_bUpdateMSDM = !slim->m_bUpdateMSDM;

	double slim_time = 0;
	int decimateF = 1000;
	while (slim->valid_faces > 700)
	{
		fprintf(outputFile, "---------------------------\n面片数:%d\n", slim->valid_faces);
		//int dF = slim->valid_faces - slim->valid_faces*0.7;
		TIMING(slim_time, slim->decimate(slim->valid_faces - decimateF));
		fprintf(outputFile, "简化%d面片用时：%lf秒\n", decimateF, slim_time);
		double msdm2Error, minMsdm, maxMsdm;
		MSDM2::GetMSDMInfo(m_pQslimMesh, msdm2Error, minMsdm, maxMsdm);
		fprintf(outputFile, "全局MSDM2：%lf；局部最小MSDM2：%lf；局部最大MSDM2：%lf\n", msdm2Error, minMsdm, maxMsdm);
		// put buffer data to disk [6/8/2012 Han]
		fflush(outputFile);

	}
	m_pQslimMesh->VertexColor(1);
}
// 使用距离作为依据,按照csf得到误差阈值,对原始模型进行相应简化
void CMeshSimpDoc::OnResampleReducebycsf()
{
	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 
	SimpByMethod t = m_SimpByMethod;
	m_SimpByMethod = SIMP_BY_DIST;
	//////////////////////////////////////////////////////////////////////////
	pView->resetOrientation();		// Rest eye pos and rotate arcball
	m_mType = QSLIM;

	float dv = pView->GetViewDistance();

	AdjustLOD(dv);
	m_SimpByMethod = t;
	pView->Invalidate();
}
// 按照用户自定义的目标面片数进行简化
void CMeshSimpDoc::OnResampleReduceuserdef()
{
	// 暂时借用这个对话框完成简化到指定面片数目的目的
	SegmentParam dlg;
	float nTargetFaceNum = 100;
	dlg.SetValues(sigma, k, nTargetFaceNum,segFace);
	if (dlg.DoModal() == IDOK)
	{
		dlg.GetValues(sigma, k, nTargetFaceNum,segFace);
	}
	else
		return;

	// test delete [7/23/2012 Han]
	// 为了计算msdm2的值，在简化时保留对应关系 [7/23/2012 Han]
	//slim->m_bUpdateMSDM = true;
	//////////////////////////////////////////////////////////////////////////

	slim->decimate(nTargetFaceNum);
}
// ↑ 为关于模型简化的函数 [8/13/2013 Han]
//////////////////////////////////////////////////////////////////////////

