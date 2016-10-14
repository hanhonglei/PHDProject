// ******************************
// pmesh.h
//
// Progressive mesh class.
// This mesh can be simplified by
// removing edges & triangles, while
// retaining the same shape.
//
// Jeff Somers
// Copyright (c) 2002
//
// jsomers@alumni.williams.edu
// March 27, 2002
// ******************************

#ifndef __PMesh_h
#define __PMesh_h



#if defined (_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#pragma warning(disable:4710) // function not inlined
#pragma warning(disable:4702) // unreachable code
#pragma warning(disable:4514) // unreferenced inline function has been removed
#pragma warning(disable:4786) /* disable "identifier was truncated to '255' characters in the browser information" warning in Visual C++ 6*/
#endif

#define WIN32_LEAN_AND_MEAN
#include "..\stdafx.h"
#include "..\Option.h"
#include <vector>
#include <list>
#include "vertex.h"
#include "triangle.h"
#include "jmsmesh.h"
#include "..\SparseMatrix.h"
#include "..\CSMatrix.h"
#include "..\MultigridContractionSolver.h"
#include "Matrix.h"
#include <afx.h>
using namespace std;


// The edge collapse structure.  The "from vertex" will
// be collapsed to the "to vertex."  This may flatten some
// triangles, which will be removed, and will affect those
// triangles which contain the "from vertex".  Those triangles
// will be updated with the new vertex.
struct EdgeCollapse
{
	int _vfrom;
	int _vto;
	set<int> _trisRemoved;
	set<int> _trisAffected;

	bool _center;


	// Used for debugging
	void dumpEdgeCollapse()
	{
		std::cout << "**** Edge Collapse Dump ****" << std::endl;

		std::cout << "\tFrom Vert# " << _vfrom << " to Vert# " << _vto << std::endl;
		cout << "\tTris removed:";
		set<int>::iterator pos;
		for (pos = _trisRemoved.begin(); pos != _trisRemoved.end(); ++pos) 
		{
			std::cout << " " << *pos;
		}
		cout << std::endl << "\tTris affected:";
		for (pos = _trisAffected.begin(); pos != _trisAffected.end(); ++pos) 
		{
			std::cout << " " << *pos;
		}
		std::cout  << std::endl << "**** End of Edge Collapse Dump ****" << std::endl;
	}
};

// This is a "pointer" to a vertex in a given mesh
struct vertexPtr
{
	jmsMesh* _meshptr;
	int _index; // ptr to vertex position in mesh
	double _center;

	bool operator<(const vertexPtr& vp) const 
	{
		return (_meshptr->getVertex(_index) < vp._meshptr->getVertex(vp._index));
	}
};


typedef multiset<vertexPtr, less<vertexPtr> > vertexPtrSet;

// 骨架点数据结构 [7/16/2011 Han Honglei]
struct VertexRecord
{
	bool center;
	set<int> collapseFrom;			// 收缩到这个骨架点的原始模型顶点编号
	int colorIndex;
	double nodeSize;
	vec3 pos;						// 骨架点位置
	int pqIndex;
	double radius;
	int vIndex;						// 本骨架点在原始模型中的顶点编号
};

// LOD控制的结构体，用来对LOD进行选择 [9/28/2011 Han Honglei]
struct LODController
{
	float startDist, vanishDist;		// 分别表示开始简化和只剩一个面片的距离，一个模型有可能是过细节的，则startDist为0
	int startTri, startVert;			// 分别表示最高细节的三角形和顶点个数
	float collapseVertPerDist;			// 每增加一个单位距离需要执行顶点合并的次数
};
// Progressive jmsMesh class.  This class will calculate and keep track
// of which vertices and triangles should be removed from/added to the
// mesh as it's simplified (or restored).
class PMesh
{
public:
	// Type of progress mesh algorithm
	enum EdgeCost {SHORTEST, MELAX, QUADRIC, QUADRICTRI, QUADRIC_SEGMENT, HAN, BLUR, MAX_EDGECOST};	// QUADRIC_SEGMENT表示添加分割的权重
	enum DrawMeshType {NONE, SKEL_MAP, COLL_DIST, SEGMENT};

	EdgeCost _cost; // Type of progressive mesh algorithm
	
	PMesh(jmsMesh* mesh, EdgeCost ec);
	PMesh(jmsMesh* mesh, EdgeCost ec, list<int> selList);

	// New collapse function, calc collapse realtime [9/20/2011 Han Honglei]
	bool collapseEdgeRealtime();
	// Collapse one vertex to another.
	bool collapseEdge();
	//bool collapseEdgeOnSel();
	
	// One vertex will be split into two vertices -- this
	// is the opposite of a collapse
	bool splitVertex();

	// number of edge collapses
	int numCollapses() {return _edgeCollList.size();}
	int numEdgeCollapses() {return _edgeCollList.size();}

	// number of triangles, and visible triangles in mesh
	int numTris() {return _newmesh.getNumTriangles();}
	int numVisTris() {return _nVisTriangles;}
	int numVertex() {return _newmesh.getNumVerts();}


	// Create the list of the edge collapses used
	// to simplify the mesh.
	void createEdgeCollapseList();


	bool getTri(int i, triangle& t) {
		t = _newmesh.getTri(i);
		return true;
	}
	triangle&  getTri(int i) {
		return _newmesh.getTri(i);
	}

	// Return a short text description of the current Edge Cost method
	char* getEdgeCostDesc();

	//  [8/18/2010 admin]
	bool SaveFile(CArchive& ar);
	bool LoadSegFiles(char *fileName);
	bool LoadBlurFiles(char *fileName);	// 载入模型顶点运动模糊的距离信息  [7/27/2011 Han Honglei]

	~PMesh();

	FILE* outputFile;							// 打开文件用于输出调试信息 [7/14/2011 Han Honglei]

	LONGLONG StartTime();						// 用于计时的函数
	double GetElapseTime(LONGLONG startTime);
	void Normalize() {_newmesh.Normalize();}
	bool AdjustLOD(float fDist);
//////////////////////////////////////////////////////////////////////////
	// 骨架抽取相关成员
	int collapseIterNum;						// 模型收缩的次数
	bool isCollapsing;							// 当前程序是否处于模型收缩阶段
	bool bCollapsed;							// 模型是否完成了收缩
	bool isSkeling;								// 当前程序是否处于骨架点抽取阶段（即模型收缩完毕，对收缩后的模型进行简化获得骨架点）
	bool bSimplified;							// 模型是否完成了简化，即是否已经提取了骨架点

	vector<VertexRecord> simplifiedVertexRec;	// 骨架点列表
	vector<set<int>> simplifiedVertexNeighbour;	// 骨架点之间的邻接关系

	bool ChangeColl(bool bNext);				// 将模型切换为每次迭代收缩后的结果
	void Simplification(CSkeOption *opt);		// 对收缩后的模型进行简化，得到骨架点
	void DrawSimplifiedVertices();				// 绘制经过简化后得到的骨架信息
	int DrawOriginalMesh(DrawMeshType drawType = NONE, bool bSmooth = true);	// 绘制原始模型
	void CalcCollapseDist();					// 计算各个模型每次收缩每个顶点移动的距离，以及模型中顶点的最大和最小移动距离
	bool HasSegment() {return _newmesh._bHasSegment;}// 模型是否有分隔信息
	bool RestoreMesh();							// 将模型恢复为原来的样子
	bool GeometryCollapse(CSkeOption *opt);		// 对模型进行递归式地折叠，直到体积收缩到一定程度，返回false
//////////////////////////////////////////////////////////////////////////
	// calc motion blur weights of model's vertexes
	bool CalcBlurWeights(/*vec3 eye,vec3 gaze,vec3 up*/);
	void meshTransform( const Matrix & mx );
private:
	// If motion blur weights of vertexes are loaded from file(.blr) [11/2/2011 Han]
	bool bIsBlurWeightsLoaded;
//////////////////////////////////////////////////////////////////////////
	// 用于抽取骨架的成员变量
	 int* adjSegmentVertex;						// 未使用
	 double* lap[3];
	 double* lapWeight;							// 进行拉普拉斯收缩的因子
	 double* posWeight;							// 在收缩过程中保持形状属性的因子
	 MultigridContractionSolver *multigridSolver;// 求解器
	 double* oldAreaRatio;						// 多边形面积比例
	 double* originalFaceArea;					// 原来的模型面积
	 double originalArea;						// 保存模型表面的原始面积
	 void* solver;								// 求解器
	 void* symbolicSolver;						// 	 
	 vector<vector<vec3>> allCollPos;			// 用来保存每次收缩后的结果 [7/6/2011 Han Honglei]
	 int currCollPos;							// 当前的收缩程度
	 vector<vector<double>> allCollapseDist;	// 保存每次收缩时，模型顶点的移动距离
	 vector<double> minCollapseDist;			// 保存每次收缩过程中，模型中顶点移动的最小距离
	 vector<double> maxCollapseDist;			// 最大距离
	 vector<set<int>> allCollapseFrom;			// 用来保存收缩后进行边折叠后的结果 [7/10/2011 Han Honglei]

	 void InitCollapse(const CSkeOption *opt);	// 在第一次收缩之前进行一些初始化工作
	 void InitSkelValues();						// 这个类构造的时候初始化一些成员变量			
	 SparseMatrix BuildMatrixA(const CSkeOption *otp);// 创建用于收缩的矩阵
	 void Dispose();							// 对动态创建的变量进行删除
	 // 用于骨架点提取的函数
	 void calcShapeMatrices(jmsMesh &mesh);		// 按照文章的方法计算每个顶点的形状矩阵，类似于Q
	 double quadricCollapseCostSkel(jmsMesh& m, vertex& v, CSkeOption* opt);// 计算折叠代价
	 void calcEdgeCollapseCostsSkel(vertexPtrSet &vertSet, vector<vertexPtrSet::iterator> &vertSetVec, 
		 int nVerts, jmsMesh &mesh, EdgeCost &cost, CSkeOption* opt);// 计算每个顶点的折叠代价，并按照代价排序，结果保存于vertSet中
	 void buildEdgeCollapseListSkel(jmsMesh &mesh, const EdgeCost &cost, 
		 list<EdgeCollapse> &edgeCollList,
		 vertexPtrSet &vertSet, 
		 vector<vertexPtrSet::iterator> &vertSetVec, CSkeOption* opt);// 按照折叠代价逐步进行边折叠
	 void insureEdgeCollapseValidSkel(EdgeCollapse &ec, vertex &vc, jmsMesh &mesh,  // 确保边折叠合法
		 const EdgeCost &cost, bool &bBadVertex, CSkeOption *opt);	
	 void updateTrianglesSkel(EdgeCollapse &ec, vertex &vc, set<int> &affectedVerts, jmsMesh &mesh);
	 void MatrixXMatrix(double m1[4][4], double m2[4][4], double result[4][4]);// 计算两个矩阵相乘
	 void MatrixXVector( double m[4][4],double v[4], double result[4]);		   // 计算矩阵和向量相乘
	 // 下面的变量是原作者中的变量，暂时不使用或者已经被替换
	 //double* originalVertexPos;
	 //jmsMesh mesh;
	 //CSkeOption opt;
	 //object myDisplayLock = new object();
	 //int remainingVertexCount;
	 //int FaceCount;
	 //int* faceIndex;
	 //double* collapsedLength;
	 //double* collapsedVertexPos;				// 这两个变量使用其他变量代替
	 //CCSMatrix ccsA;
	 //CCSMatrix ccsATA;
	 //bool displayIntermediateMesh;
	 //bool displayNodeSphere;
	 //bool displayOriginalMesh;
	 //bool displaySimplifiedMesh;
	 //int displaySimplifiedMeshIndex;
	 //bool displaySkeleton;
	 //VertexRecord rootNode;
	 //List<VertexRecord> simplifiedVertexRec;
	 //float skeletonNodeSize;
	 //VertexRecord[] vRec;
//////////////////////////////////////////////////////////////////////////
	jmsMesh* _mesh; // original mesh - not changed
	jmsMesh _newmesh; // we change this one

	bool _bSel;		// 当前简化是否处于“选定面片简化”模式

	LODController _lod;	// 对模型的lod级别进行控制

	list<int> _selList;

	// This is a set of vertex pointers, ordered by edge collapse cost.
	vertexPtrSet _vertSet;
	vector<vertexPtrSet::iterator> _vertSetVec;

	list<EdgeCollapse> _edgeCollList; // list of edge collapses
	list<EdgeCollapse>::iterator _edgeCollapseIter;

	// functions used to calculate edge collapse costs.  Different
	// methods can be used, depending on user preference.
	double shortEdgeCollapseCost(jmsMesh& m, vertex& v);
	double melaxCollapseCost(jmsMesh& m, vertex& v);
	double quadricCollapseCost(jmsMesh& m, vertex& v);

	int _nVisTriangles; // # of triangles, after we collapse edges

	// Used in the QEM edge collapse methods.
	void calcAllQMatrices(jmsMesh& mesh, bool bUseTriArea); // used for quadric method
	double calcQuadricError(double Qsum[4][4], vertex& v, double triArea); // used for quadric method

	enum {BOUNDARY_WEIGHT = 1000}; // used to weight border edges so they don't collapse
	void applyBorderPenalties(set<border> &borderSet, jmsMesh &mesh);

	PMesh(const PMesh&); // don't allow copy ctor -- too expensive
	PMesh& operator=(const PMesh&); // don't allow assignment op.
	bool operator==(const PMesh&); // don't allow op==

#ifndef NDEBUG
	// used in debugging
	void assertEveryVertActive(int nVerts, int nTri, jmsMesh &mesh);
#endif
	// helper function for edge collapse costs
	void calcEdgeCollapseCosts(vertexPtrSet &vertSet, vector<vertexPtrSet::iterator> &vertSetVec, 
								  int nVerts, jmsMesh &mesh, EdgeCost &cost);

	// Calculate the QEM matrices used to computer edge
	// collapse costs.
	void calcQuadricMatrices(EdgeCost &cost, jmsMesh &mesh);

	// We can't collapse Vertex1 to Vertex2 if Vertex2 is invalid.
	// This can happen if Vertex2 was previously collapsed to a
	// separate vertex.
	void insureEdgeCollapseValid(EdgeCollapse &ec, vertex &vc, jmsMesh &mesh, 
									const EdgeCost &cost, bool &bBadVertex);

	// Calculate the QEM for the "to vertex" in the edge collapse.
	void setToVertexQuadric(vertex &to, vertex &from, const EdgeCost &cost);

	// At this point, we have an edge collapse.  We're collapsing the "from vertex"
	// to the "to vertex."  For all the surrounding triangles which use this edge, 
	// update "from vertex" to the "to vertex".  Also keep track of the vertices
	// in the surrounding triangles. 
	void updateTriangles(EdgeCollapse &ec, vertex &vc, set<int> &affectedVerts, jmsMesh &mesh);


	// These affected vertices are not in the current collapse, 
	// but are in the triangles which share the collapsed edge.
	void updateAffectedVertNeighbors(vertex &vert, const EdgeCollapse &ec, 
		set<int> &affectedVerts);

	// Reset the edge collapse costs of vertices which were
	// affected by a previous edge collapse.
	void resetAffectedVertCosts(const EdgeCost &cost, jmsMesh &newmesh, vertex &vert);

	// If this vertex has no active triangles (i.e. triangles which have
	// not been removed from the mesh) then set it to inactive.
	void removeVertIfNecessary(vertex &vert, vertexPtrSet &vertSet, 
								  vector<vertexPtrSet::iterator> &vertSetVec, 
								  jmsMesh &mesh, const EdgeCost &cost, 
									set<int> &affectedQuadricVerts);

	// Update the vertices affected by the most recent edge collapse
	void updateAffectedVerts(jmsMesh &_newmesh, vector<vertexPtrSet::iterator> &vertSetVec, 
							vertexPtrSet &vertSet, const EdgeCollapse &ec, 
							set<int> &affectedVerts, const EdgeCost &cost, 
							set<int> &affectedQuadricVerts);

	// Recalculate the QEM matrices (yeah, that's redundant) if we're
	// using the Quadrics to calculate edge collapse costs.
	void recalcQuadricCollapseCosts(set<int> &affectedQuadricVerts, 
								   jmsMesh &mesh, const EdgeCost &cost);

	// Calculate the list of edge collapses.  Each edge collapse
	// consists of two vertices:  a "from vertex" and a "to vertex".
	// The "from vertex" is collapsed to the "to vertex".  The
	// "from vertex" is removed from the mesh.
	void buildEdgeCollapseList(jmsMesh &mesh, const EdgeCost &cost, 
							  list<EdgeCollapse> &_edgeCollList,
								vertexPtrSet &vertSet, 
								vector<vertexPtrSet::iterator> &vertSetVec);

	// Helper function for melaxCollapseCost().  This function
	// will loop through all the triangles to which this vertex
	// belongs.
	void calcMelaxMaxValue(jmsMesh &mesh, set<int> &adjfaces, 
							  vertex &v, set<int> &tneighbors,
								float &retmaxValue, 
								bool &bMaxValueFound);
public:
	void InitOpt(CSkeOption *opt);
};

#endif // __PMesh_h
