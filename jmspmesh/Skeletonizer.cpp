#include <assert.h>
#include <float.h>

#if defined (_MSC_VER) && (_MSC_VER >= 1020)
#pragma warning(disable:4710) // function not inlined
#pragma warning(disable:4702) // unreachable code
#pragma warning(disable:4514) // unreferenced inline function has been removed
#pragma warning(disable:4786) /* disable "identifier was truncated to '255' characters in the browser information" warning in Visual C++ 6*/
#endif

#include <set>
#include <map>
#include <ostream>

using namespace std;

#include "pmesh.h"
//#include "..\MainFrm.h"
#include <float.h>
#include <math.h>
#include <gl/glut.h>

// 用来初始化一些参数信息 [6/21/2011 Han Honglei]
void PMesh::InitOpt(CSkeOption *opt)
{
	opt->laplacianConstraintWeight = 1.0 / (10.0 * sqrt(_newmesh.AverageFaceArea()));

}
void PMesh::InitCollapse(const CSkeOption *opt)
{
	//displaySkeleton = true;
	//skeletonNodeSize = 6.0f;
	solver = NULL;
	symbolicSolver = NULL;
	collapseIterNum = 0;
	isCollapsing = true;

	int n = _newmesh.getNumVerts();
	int fn = _newmesh.getNumTriangles();
	//opt = opt;
	//if (opt->AddNoise)
	//{
	//	AddNoise();
	//}
	if (lapWeight == NULL)
		lapWeight = new double[n];
	if (posWeight == NULL)
		posWeight = new double[n];
	if (originalFaceArea == NULL)
		originalFaceArea = new double[fn];
	//originalVertexPos = (double[]) mesh.VertexPos.Clone();		// 暂时不需要 [6/21/2011 Han Honglei]
	for (int i = 0; i < fn; i++)
	{
		originalFaceArea[i] = abs(_newmesh.ComputeFaceArea(i));
	}
	// 将收缩结果保存 [7/6/2011 Han Honglei]
	vector<vec3> vecPos;
	vecPos.resize(n);
	allCollPos.clear();
	allCollPos.push_back(vecPos);
	currCollPos = allCollPos.size()-1;
	for (int i = 0; i < n; i++)
	{
		lapWeight[i] = opt->laplacianConstraintWeight;
		posWeight[i] = opt->positionalConstraintWeight;
		allCollPos[currCollPos][i] = _newmesh.getVertex(i).getXYZ(); // 将最初的状态保存
	}
	if (opt->useIterativeSolver)
	{
		if (multigridSolver != NULL)
			delete multigridSolver;
		multigridSolver = new MultigridContractionSolver(&_newmesh);
	}
}

// 来源自sig08文章 [6/14/2011 Han Honglei]

void PMesh::InitSkelValues()
{
	//////////////////////////////////////////////////////////////////////////
	// 做一些骨架抽取的初始化工作
	adjSegmentVertex = NULL;
	for (int i = 0; i < 3; i++)
	{
		lap[i] = NULL;
	}
	//lap = NULL;
	lapWeight = NULL;
	oldAreaRatio = NULL;
	originalFaceArea = NULL;
	posWeight = NULL;
	multigridSolver =NULL;
	solver =NULL;
	symbolicSolver =NULL;
	collapseIterNum = 0;
	isCollapsing = false;
	currCollPos = 0;

	bSimplified = false;
	bCollapsed = false;
	isSkeling = false;
	outputFile = fopen("骨架信息.txt", "wt");
}

void PMesh::Dispose()
{
	// 利用动态链接库使用函数 [6/26/2011 Han Honglei]
	if (multigridSolver != NULL)
		if (multigridSolver->hDllInst == NULL)
			multigridSolver->hDllInst = LoadLibrary("taucs.DLL");
	if(multigridSolver->hDllInst)
	{
		typedef int (* MYFUNC)(void*);
		MYFUNC FreeSolver = NULL; 
		FreeSolver = (MYFUNC)GetProcAddress
			(multigridSolver->hDllInst,"FreeSolver");
		if(FreeSolver)
		{
			if (solver != NULL)
			{
				FreeSolver(solver);	// 这个是taucs库函数
				solver = NULL;
			}
		}

		MYFUNC FreeSymbolicSolver = NULL; 
		FreeSymbolicSolver = (MYFUNC)GetProcAddress
			(multigridSolver->hDllInst,"FreeSymbolicSolver");
		if(FreeSymbolicSolver)
		{
			if (symbolicSolver != NULL)
			{
				FreeSymbolicSolver(symbolicSolver);	// 这个是taucs库函数
				symbolicSolver = NULL;
			}

		}
	}
}

// 这个函数完成收缩权重和位置权重的更新
SparseMatrix PMesh::BuildMatrixA(const CSkeOption *opt)
{
	int n = _newmesh.getNumVerts()/*mesh.VertexCount*/;
	int fn = _newmesh.getNumTriangles()/*mesh.FaceCount*/;
	SparseMatrix A(n, n);
	double* areaRatio = new double[fn];
	bool* collapsed = new bool[n];
	if (oldAreaRatio == NULL)
	{
		oldAreaRatio = new double[fn];
		for (int i = 0; i < fn; i++)
		{
			oldAreaRatio[i] = 0.4;
		}
	}
	int i = 0;
	for (int j = 0; i < fn; j += 3)
	{
		int c1 = getTri(i).getVert1Index();
		int c2 = getTri(i).getVert2Index();
		int c3 = getTri(i).getVert3Index();
		vec3 v1 = vec3(getTri(i).getVert1());
		vec3 v2 = vec3(getTri(i).getVert2());
		vec3 v3 = vec3(getTri(i).getVert3());

		areaRatio[i] = /*Math.Abs*/abs(/*mesh*/_newmesh.ComputeFaceArea(i)) / originalFaceArea[i];
		double areaRatioThreshold = opt->areaRatioThreshold;
		double num2 = areaRatio[i];

		vec3 CS0 = v2 - v1;
		vec3 CS1 = v2 - v1;
		double cot1 = CS0.dot(v3 - v1) / CS1.cross(v3 - v1).length();
		vec3 CS3 = v3 - v2;
		vec3 CS4 = v3 - v2;
		double cot2 = CS3.dot(v1 - v2) / CS4.cross(v1 - v2).length();
		vec3 CS6 = v1 - v3;
		vec3 CS7 = v1 - v3;
		double cot3 = CS6.dot(v2 - v3) / CS7.cross(v2 - v3).length();
		//try
		//{
		//	if (double.IsNaN(cot1))
		//	{
		//		throw new Exception();
		//	}
		//	if (double.IsNaN(cot2))
		//	{
		//		throw new Exception();
		//	}
		//	if (double.IsNaN(cot3))
		//	{
		//		throw new Exception();
		//	}
		//}
		//catch (Exception e)
		//{
		//	Program.PrintText(e.Message);
		//	Program.PrintText(string.Concat(new object[] { "!!! ", cot1, " ", cot2, " ", cot3 }));
		//	cot1 = cot2 = cot3 = 0.0;
		//}
		A.AddValueTo(c2, c2, -cot1);
		A.AddValueTo(c2, c3, cot1);
		A.AddValueTo(c3, c3, -cot1);
		A.AddValueTo(c3, c2, cot1);
		A.AddValueTo(c3, c3, -cot2);
		A.AddValueTo(c3, c1, cot2);
		A.AddValueTo(c1, c1, -cot2);
		A.AddValueTo(c1, c3, cot2);
		A.AddValueTo(c1, c1, -cot3);
		A.AddValueTo(c1, c2, cot3);
		A.AddValueTo(c2, c2, -cot3);
		A.AddValueTo(c2, c1, cot3);
		i++;
	}
	double count = 0.0;
	for (int i = 0; i < n; i++)
	{
		double totRatio = 0.0;
		double oldTotRatio = 0.0;
		// adjFF、adjVF、adjVV应该分别代表面面相邻、顶点面相邻和顶点顶点相邻矩阵
		set<int>& triset =_newmesh.getVertex(i).getTriNeighbors();
		set<int>::iterator iter;
		for (iter = triset.begin(); iter != triset.end(); ++iter)
		{
			totRatio += areaRatio[*iter];
			oldTotRatio += oldAreaRatio[*iter];
		}

		totRatio /= (double) triset.size()/*mesh.AdjVF[i].Length*/;
		oldTotRatio /= (double) triset.size()/*mesh.AdjVF[i].Length*/;
		double tot = 0.0;
		tot = 0.0;

		vector<vector<Element>>::iterator itARow = A.Rows();
		for (vector<Element>::iterator it = itARow[i].begin(); it < itARow[i].end(); it++)
		{
			if (it->i != it->j)
			{
				tot += it->value;
			}
		}

		if (tot > 10000.0)
		{
			collapsed[i] = true;
			//_newmesh.Flag[i] = 1;	// 暂时删除
			for (vector<Element>::iterator it = itARow[i].begin(); it < itARow[i].end(); it++)
			{
				it->value /= tot / 10000.0;
			}
			//mesh.Flag[i] = 1;
			//foreach (SparseMatrix.Element e in A.Rows[i])
			//{
			//	e.value /= tot / 10000.0;
			//}
		}
		for (vector<Element>::iterator it = itARow[i].begin(); it < itARow[i].end(); it++)
		{
			it->value *= lapWeight[i];
		}

		lapWeight[i] *= /*opt->*/opt->laplacianConstraintScale;
		if (lapWeight[i] > 2048.0)
		{
			lapWeight[i] = 2048.0;
		}
		double d = (1.0 / /*Math.Sqrt*/sqrt(totRatio)) * /*opt->*/opt->positionalConstraintWeight;
		if (!_isnan(d)/*double.IsNaN(d)*/)
		{
			posWeight[i] = d;
		}
		if (posWeight[i] > 10000.0)
		{
			posWeight[i] = 10000.0;
		}
		count++;
		bool ok = true;
		for (vector<Element>::iterator it = itARow[i].begin(); it < itARow[i].end(); it++)
		{
			if (_isnan(it->value))
			{
				ok = false;
			}
		}
		if (!ok)
		{
			for (vector<Element>::iterator it = itARow[i].begin(); it < itARow[i].end(); it++)
			{
				if (it->i == it->j)
				{
					it->value = -1.0;
				}
				else
				{
					it->value = 1.0 / ((double)_newmesh.getVertex(i).getVertNeighbors().size() /*mesh.AdjVV[i].Length*/);
				}
			}
		}
	}
	double* x = new double[n];
	for (int i = 0; i < 3; i++)
	{
		if (lap[i] == NULL)
			lap[i] = new double[n];
		int j = 0;
		for (int k = i; j < n; k += 3)
		{
			switch (i)
			{
			case 0:
				x[j] = _newmesh.getVertex(j).getXYZ().x/*mesh.VertexPos[k]*/;
				break;
			case 1:
				x[j] = _newmesh.getVertex(j).getXYZ().y/*mesh.VertexPos[k]*/;
				break;
			default:
				x[j] = _newmesh.getVertex(j).getXYZ().z/*mesh.VertexPos[k]*/;
				break;
			}
			j++;
		}
		A.Multiply(x, lap[i]);
	}
	delete[] x;
	for (int i = 0; i < fn; i++)
	{
		oldAreaRatio[i] = areaRatio[i];
	}
	for (int i = 0; i < n; i++)
	{
		A.AddRow();
		A.AddElement(i + n, i, posWeight[i]);
	}
	for (int i = 0; i < n; i++)
	{
		A.AddRow();
		A.AddElement((i + n) + n, i, opt->originalPositionalConstraintWeight);
	}
	A.SortElement();

// 做清理工作 [7/3/2011 Han Honglei]
	delete[] areaRatio;
	delete[] collapsed;

	return A;
}

// 开始计时，确定初始时刻
LONGLONG PMesh::StartTime()
{
	LONGLONG   t1 = 0;
	QueryPerformanceCounter((LARGE_INTEGER   *)&t1);  
	return t1;
}

double PMesh::GetElapseTime(LONGLONG startTime)
{
	LONGLONG   t2 = 0;   
	LONGLONG   persecond = 0;   
	QueryPerformanceFrequency((LARGE_INTEGER   *)&persecond);//询问系统一秒钟的频率   
	double perms = persecond;
 
	QueryPerformanceCounter((LARGE_INTEGER   *)&t2);  
	double   time=double(t2-startTime)/perms;			// 求得栅格化所用时间
	return time;	
}

// 通过拉普拉斯平滑，对原始模型进行塌陷 [6/13/2011 Han Honglei]
bool PMesh::GeometryCollapse(CSkeOption *opt)
{
	LONGLONG sTime = StartTime();

	bool bOk = false;
	// 添加初始化操作 [6/26/2011 Han Honglei]
	if (opt->firstTime)
	{
		InitCollapse((CSkeOption *)opt);
		opt->originalVolume = _newmesh.Volume();
		for (int i = 0; i < /*mesh.FaceCount*/_newmesh.getNumTriangles(); i++)
		{
			originalArea += /*mesh.*//*_newmesh.ComputeFaceArea(i)*/originalFaceArea[i];
		}
		opt->firstTime = false;
	}

	double	currentVolume = _newmesh.Volume();
	// 暂时删除，未用到 [7/5/2011 Han Honglei]
	//double	currentArea = 0.0;
	//for (int i = 0; i < _newmesh.getNumTriangles()/*mesh.FaceCount*/; i++)
	//{
	//	//currentArea += mesh.ComputeFaceArea(i);
	//	currentArea += _newmesh.ComputeFaceArea(i);
	//}
	if (((currentVolume / opt->originalVolume) > opt->volumneRatioThreashold) && (opt->maxIterations > 0))
	{
			if (opt->useIterativeSolver)
			{
				double** pos = multigridSolver->SolveSystem(lapWeight, posWeight);
				int i = 0;
				vector<vec3> vecPos;
				vecPos.resize(_newmesh.getNumVerts());
				allCollPos.push_back(vecPos);
				currCollPos = allCollPos.size()-1;
				// 将收缩结果保存 [7/6/2011 Han Honglei]
				for (int j = 0; i < _newmesh.getNumVerts(); j += 3)
				{
					_newmesh.setVerPos(i,pos[0][i],pos[1][i],pos[2][i]);
					allCollPos[currCollPos][i].x = pos[0][i];
					allCollPos[currCollPos][i].y = pos[1][i];
					allCollPos[currCollPos][i].z = pos[2][i];
					i++;
				}
				BuildMatrixA(opt);
				_newmesh.calcTriNormals();
				_newmesh.calcVertNormals();
				// 输出用时 [7/14/2011 Han Honglei]
				if (outputFile != NULL)
				{
					double elapseT = GetElapseTime(sTime);
					fprintf(outputFile, "第%d次收缩用时为：%lf秒\n", collapseIterNum, elapseT);
				}
				opt->maxIterations--;
				collapseIterNum++;
				isCollapsing = true;
				for (int i = 0; i < 3; i++)
				{
					delete[] pos[i];
				}
				delete[] pos;
			}
			// 下面语句暂时删除，用来测试前面基本代码的正确性 [6/21/2011 Han Honglei]
			//else
			//{
			//	SparseMatrix A = BuildMatrixA();
			//	ccsA = new CCSMatrix(A);
			//	ccsATA = MultiplyATA(ccsA);
			//	A = NULL;
			//	//GC.Collect();
			//	if (opt->UseSymbolicSolver)
			//	{
			//		symbolicSolver = SymbolicFactorization(ccsATA);
			//		//Program.PrintText("Symbolic solver: " + ((symbolicSolver != NULL)).ToString());
			//		int ret = NumericFactorization(symbolicSolver, ccsATA);
			//		//Program.PrintText("Numeric solver: " + ((ret == 0)).ToString());
			//	}
			//	else
			//	{
			//		if (solver != NULL)
			//		{
			//			FreeSolver(solver);										// 这个是taucs库函数
			//		}
			//		solver = Factorization(ccsATA);
			//		if (solver == NULL)
			//		{
			//			//throw new Exception();
			//			return;		// 发生错误
			//		}
			//		//GC.Collect();
			//	}
			//	ImplicitSmooth();
			//}

			//timer.Stop();

			//Program.PrintText(string.Concat(new object[] { iter, " Area: ", currentArea / originalArea, " Vol: ", currentVolume / originalVolume, " CPU Time: ", timer.Duration }));
			//Thread.Sleep(0);
			//CMainFrame *p=(CMainFrame*)AfxGetMainWnd();
			//CView *pv=p->GetActiveView();
			//pv->Invalidate();
		Dispose();
		bOk = true;
	}
	return bOk;
}
bool PMesh::ChangeColl(bool bNext)
{
	if (allCollPos.size() == 0 || currCollPos < 0 || currCollPos >=allCollPos.size())
		return false;
	if (bNext)
	{
		if ((currCollPos + 1) >= allCollPos.size())
			return false;
		currCollPos++;

	}
	else
	{
		if ((currCollPos - 1) < 0)
			return false;
		currCollPos--;
	}
	for (int i = 0; i < _newmesh.getNumVerts(); i++)
	{
		_newmesh.setVerPos(i, allCollPos[currCollPos][i]);
	}
	_newmesh.calcTriNormals();
	_newmesh.calcVertNormals();

}

void PMesh::calcShapeMatrices(jmsMesh &mesh)
{

	double matrix[4][4];

	for (int i = 0; i < mesh.getNumVerts(); i++)
	{
		vertex &v = mesh.getVertex(i);
		memset(matrix, 0.0, 16*sizeof(double));

		set<int>& neighbors = v.getVertNeighbors();
		set<int>::iterator pos;
		for (pos = neighbors.begin(); pos != neighbors.end(); ++pos) 
		{
			vertex& n = mesh.getVertex(*pos);
			if (!n.isActive()) continue;
			if (n == v) continue;

			// 计算矩阵
			vec3 u = n.getXYZ() - v.getXYZ();
			u.normalize();
			vec3 w = u.cross(v.getXYZ());
			double m[4][4];
			memset(m, 0.0, 16*sizeof(double));
			m[0][1] = -u.z;
			m[0][2] = u.y;
			m[0][3] = -w.x;
			m[1][0] = u.z;
			m[1][2] = -u.x;
			m[1][3] = -w.y;
			m[2][0] = -u.y;
			m[2][1] = u.x;
			m[2][3] = -w.z;

			double mT[4][4];
			memset(mT, 0.0, 16*sizeof(double));
			mT[1][0] = -u.z;
			mT[2][0] = u.y;
			mT[3][0] = -w.x;
			mT[0][1] = u.z;
			mT[2][1] = -u.x;
			mT[3][1] = -w.y;
			mT[0][2] = -u.y;
			mT[1][2] = u.x;
			mT[3][2] = -w.z;

			double resultM[4][4];

			MatrixXMatrix(mT, m, resultM);

			for (int i = 0; i < 4; i ++)
				for (int j = 0; j < 4; j++)
					matrix[i][j] += resultM[i][j];

		}
		v.setQuadric(matrix);
	}

}

// 重写这个函数来得到sig08的折叠误差：e = Wa*Fa+Wb*Fb，Fa类似于Q这样的形状误差；而Fb是避免出现长边的误差
double PMesh::quadricCollapseCostSkel(jmsMesh& m, vertex& v, CSkeOption* opt)
{
	// get list of all active neighbors
	// calculate quadric cost
	double mincost = FLT_MAX; // from float.h
	bool bNeighborFound = false;
	double cost = 0;

	double Q1[4][4];
	v.getQuadric(Q1);
	set<int>& neighbors = v.getVertNeighbors();
	set<int>::iterator pos;

	v.setEdgeRemoveCost(FLT_MAX);
	v.setMinCostEdgeVert(-1);
	if (neighbors.size() <= 1 || v.getTriNeighbors().size() <= 0)
		return FLT_MAX;

	double totLength = 0.0;
	for (pos = neighbors.begin(); pos != neighbors.end(); ++pos) 
	{
		vertex& n = m.getVertex(*pos);
		if (!n.isActive()) continue;
		if (n == v) continue;
		vec3 CS0001 = v.getXYZ() - n.getXYZ();
		totLength += CS0001.length();
	}
	totLength /= (double) v.getVertNeighbors().size();

	for (pos = neighbors.begin(); pos != neighbors.end(); ++pos) 
	{
		bool found = false;
		for (set<int>::iterator index = v.getTriNeighbors().begin(); index != v.getTriNeighbors().end(); index++)
		{
			triangle &t = m.getTri(*index);
			if (t.hasVertex(*pos))
			{
				found = true;
			}
		}
		if (found)
		{
			vertex& n = m.getVertex(*pos);
			if (!n.isActive()) continue;
			if (n == v) continue;

			if ((n.getTriNeighbors()).size() != 0)
			{
				cost = 0;
				if (opt->useSamplingEnergy)
				{
					vec3 CS0004 = v.getXYZ() - n.getXYZ();
					cost = CS0004.length() * totLength;
				}
				if (opt->useShapeEnergy)
				{
					double v1[4] = {n.getXYZ().x, n.getXYZ().y,n.getXYZ().z,1.0};
					vec3 tmp = ((v.getXYZ() + n.getXYZ()) / 2.0);
					double v2[4] = {tmp.x, tmp.y, tmp.z, 1.0};
					double mat1[4][4];
					double mat2[4][4];
					double mat[4][4];
					v.getQuadric(mat1);
					n.getQuadric(mat2);
					for (int i = 0; i < 4; i ++)
						for (int j = 0; j < 4; j++)
							mat[i][j] = mat1[i][j]+mat2[i][j];
					double R1[4], R2[4];
					MatrixXVector(mat, v1, R1);
					MatrixXVector(mat, v2, R2);
					double tD1 = 0, tD2 = 0;
					for (int i = 0; i < 4; i++)
					{
						tD1 += R1[i] * v1[i];
						tD2 += R2[i] * v2[i];
					}
					double e1 = tD1 * opt->shapeEnergyWeight;
					double e2 = tD2 * opt->shapeEnergyWeight;

					if (e1 < e2)
					{
						cost += e1;
					}
					else
					{
						cost += e2;
						opt->reserved = true;		// 用于传递给外面的顶点记录，表明这个点是center
						v._center = true;
					}
				}
			}

			if (cost < mincost)
			{
				bNeighborFound = true;
				mincost = cost;
				v.setEdgeRemoveCost(cost);
				v.setMinCostEdgeVert(*pos);
				assert(v.minCostEdgeVert() >= 0 && v.minCostEdgeVert() < m.getNumVerts());
			}
		}
	}
	if (bNeighborFound) {
		return mincost;
	} else {
		return FLT_MAX; // vertex not connected to an edge
	}
}

// Calculate edge collapse costs.  Edges with low costs
// are collapsed first.
void PMesh::calcEdgeCollapseCostsSkel(vertexPtrSet &vertSet, vector<vertexPtrSet::iterator> &vertSetVec, 
									  int nVerts, jmsMesh &mesh, EdgeCost &cost, CSkeOption* opt)
{
	int i;
	for (i = 0; i < nVerts; ++i)
	{
		vertex& currVert = mesh.getVertex(i);

		opt->reserved = currVert._center;			// 用来判断当前折叠点是否是center

		quadricCollapseCostSkel(mesh, currVert, opt);

		vertexPtr v;
		v._index = i;
		v._meshptr = &mesh;
		v._center = opt->reserved;

		vertSetVec[i] = vertSet.insert(v); // inserts a copy
	}

#ifdef PRINT_DEBUG_INFO
	int count=0; // for debug
	std::cout << "---- Initial State ----" << std::endl;
	mesh.dump();
	dumpset(vertSet);
	std::cout << "---- End Initial State ----" << std::endl;
#endif 

}
// At this point, we have an edge collapse.  We're collapsing the "from vertex"
// to the "to vertex."  For all the surrounding triangles which use this edge, 
// update "from vertex" to the "to vertex".  Also keep track of the vertices
// in the surrounding triangles. 
void PMesh::updateTrianglesSkel(EdgeCollapse &ec, vertex &vc, set<int> &affectedVerts, jmsMesh &mesh)
{
	//set<int>& triNeighbors = vc.getTriNeighbors();
	//set<int>& verNeighbors = vc.getVertNeighbors();
	//set<int>::iterator verIter = verNeighbors.begin();
	//set<int>::iterator pos;

	//for (pos = triNeighbors.begin(); pos != triNeighbors.end(); ++pos) 
	//{
	//	// get triangle
	//	int triIndex = *pos;
	//	triangle& t = mesh.getTri(triIndex);
	//	if (!t.isActive()) continue;
	//	bool bRemoveTri = false;
	//	int vert1, vert2, vert3;
	//	t.getVerts(vert1, vert2, vert3);
	//	if (t.hasVertex(ec._vto) || vert1 == vert2 || vert2 == vert3 || vert3 == vert1)
	//	{
	//		while (verIter != verNeighbors.end())
	//		{
	//			mesh.getVertex(*verIter)->removeTriNeighbor(triIndex);
	//			this.vRec[index2].adjF.Remove(index);
	//			verIter++;
	//		}
	//		facesLeft--;
	//	}
	//	else
	//	{
	//		t.changeVertex(ec._vfrom, ec._vto); // update the vertex of this triangle
	//		t.calcNormal(); // reset the normal for the triangle

	//		// make sure the "to" vertex knows about this triangle
	//		mesh.getVertex(ec._vto).addTriNeighbor(triIndex);
	//	}
	//}
}

// Calculate the list of edge collapses.  Each edge collapse
// consists of two vertices:  a "from vertex" and a "to vertex".
// The "from vertex" is collapsed to the "to vertex".  The
// "from vertex" is removed from the mesh.
void PMesh::buildEdgeCollapseListSkel(jmsMesh &mesh, const EdgeCost &cost, 
								  list<EdgeCollapse> &edgeCollList,
								  vertexPtrSet &vertSet, 
								  vector<vertexPtrSet::iterator> &vertSetVec, CSkeOption* opt)
{
	int facesLeft = _newmesh.getNumTriangles();/*this.mesh.FaceCount;*/
	int vertexLeft = _newmesh.getNumVerts();/*this.mesh.VertexCount;*/
	int edgeLeft = 0;
	int remainingVertexCount = vertexLeft;
	for (vertexPtrSet::iterator iter = vertSet.begin(); iter != vertSet.end(); iter++)
	{
		vertex &v = _newmesh.getVertex(iter->_index);
		edgeLeft += v.getVertNeighbors().size();
	}
	edgeLeft /= 2;

	for (;;)
	{
		if (0 == vertSet.size() || facesLeft <= 0||vertexLeft <= opt->targetVertexCount)
		{
			// we're done
			break;
		}
		const vertexPtr vp = *(vertSet.begin()); // This is a copy of the first element
		vertex &vc = mesh.getVertex(vp._index);
		assert(vp._index == vc.getIndex());

		EdgeCollapse ec; // create EdgeCollapse structure
		ec._center = vp._center;

		ec._vfrom = vc.getIndex(); // we'll collapse this vertex...
		ec._vto = vc.minCostEdgeVert(); // to this one

		if (-1 == vc.minCostEdgeVert())
			continue;

		mesh.getVertex(ec._vfrom).setActive(false);
		vertSet.erase(vertSet.begin());

		vertex& to = mesh.getVertex(ec._vto);
		vertex& from = mesh.getVertex(ec._vfrom);

		// 添加对center的支持 [7/10/2011 Han Honglei]
		if (vp._center)
		{
			// 这个会影响最终生成的骨架点位置
			to.setPos((to.getXYZ() + from.getXYZ())/2.0);
		}

		setToVertexQuadric(to, from, QUADRIC);

		// 保留折叠记录 [7/10/2011 Han Honglei]
		if (bCollapsed && isSkeling)
		{
			allCollapseFrom[ec._vto].insert(ec._vfrom);
			set<int>::iterator fromIter;
			for (fromIter = allCollapseFrom[ec._vfrom].begin(); fromIter != allCollapseFrom[ec._vfrom].end(); ++fromIter) 
				allCollapseFrom[ec._vto].insert(*fromIter);
			allCollapseFrom[ec._vfrom].clear();
			allCollapseFrom[ec._vfrom].insert(-1);		// 以此证明这个点被删除了
		}

		set<int>& triNeighbors = from.getTriNeighbors();
		set<int>& verNeighbors = from.getVertNeighbors();
		set<int>::iterator verIter = verNeighbors.begin();
		set<int>::iterator pos;
		for (pos = triNeighbors.begin(); pos != triNeighbors.end(); ++pos) 
		{
			// get triangle
			int triIndex = *pos;
			triangle& t = mesh.getTri(triIndex);
			bool bRemoveTri = false;
			int vert1, vert2, vert3;
			t.getVerts(vert1, vert2, vert3);
			if (t.hasVertex(ec._vto) || vert1 == vert2 || vert2 == vert3 || vert3 == vert1)
			{
				for (verIter = verNeighbors.begin(); verIter != verNeighbors.end(); verIter++)
				{
					mesh.getVertex(*verIter).removeTriNeighbor(triIndex);
				}
				facesLeft--;
			}
			else
			{
				t.changeVertex(ec._vfrom, ec._vto); // update the vertex of this triangle
				//t.calcNormal(); // reset the normal for the triangle

				// make sure the "to" vertex knows about this triangle
				to.addTriNeighbor(triIndex);
			}
		}

		for (verIter = from.getVertNeighbors().begin(); verIter != from.getVertNeighbors().end(); verIter++/*int index in rec1.adjV*/)
		{
			if (mesh.getVertex(*verIter).hasVertNeighbor(from.getIndex()))
			{
				mesh.getVertex(*verIter).removeVertNeighbor(from.getIndex());
				edgeLeft--;
			}
			if (*verIter != to.getIndex()/*index != r2*/)
			{
				mesh.getVertex(*verIter).addVertNeighbor(to.getIndex());
				to.addVertNeighbor(*verIter);
			}
		}
		for (verIter = to.getVertNeighbors().begin(); verIter != to.getVertNeighbors().end(); verIter++)
		{
			vertex& currV = mesh.getVertex(*verIter);
			opt->reserved = currV._center;			// 用来判断当前折叠点是否是center

			vertSet.erase(vertSetVec[*verIter]);

			quadricCollapseCostSkel(mesh, currV, opt);

			vertexPtr v;
			v._index = *verIter;
			v._meshptr = &mesh;
			v._center = opt->reserved;

			vertSetVec[*verIter] = vertSet.insert(v); // inserts a copy
		}
		int nI = to.getIndex();
		vertSet.erase(vertSetVec[nI]);
		opt->reserved = to._center;			// 用来判断当前折叠点是否是center
		quadricCollapseCostSkel(mesh, to, opt);
		vertexPtr v;
		v._index = nI;
		v._meshptr = &mesh;
		v._center = opt->reserved;
		vertSetVec[nI] = vertSet.insert(v); // inserts a copy

		for (pos = to.getTriNeighbors().begin(); pos != to.getTriNeighbors().end(); pos++)
		{
			int vert1, vert2, vert3;
			mesh.getTri(*pos).getVerts(vert1, vert2, vert3);
			if (vert3 == vert1 || vert1 == vert2 || vert2 == vert3)
			{
				to.removeTriNeighbor(*pos);
				for (verIter = to.getVertNeighbors().begin(); verIter != to.getVertNeighbors().end(); verIter++)
				{
					mesh.getVertex(*verIter).removeTriNeighbor(*pos);
				}
				facesLeft--;
			}
		}

		vertexLeft--;

		remainingVertexCount = vertexLeft;
	}
}

bool PMesh::RestoreMesh()
{
	_newmesh = *_mesh;
	if (_bSel)
	{
		list<int>::iterator selIter;

		for (selIter = _selList.begin(); selIter != _selList.end(); ++selIter) 
		{
			triangle &t = _newmesh.getTri(*selIter);
			t.setSel(true);
		}
	}
	return true;

}

void PMesh::insureEdgeCollapseValidSkel(EdgeCollapse &ec, vertex &vc, jmsMesh &mesh, 
	const EdgeCost &cost, bool &bBadVertex, CSkeOption *opt)
{
	int nLoopCount = 0;
	for (;;) // at most this will loop twice -- the "to vertex" could have been removed, so it may need to be recalculated.
	{
		++nLoopCount;
		ec._vfrom = vc.getIndex(); // we'll collapse this vertex...
		ec._vto = vc.minCostEdgeVert(); // to this one

		if (-1 == vc.minCostEdgeVert())
		{
			// this point is isolated -- it's not connected to another point by an edge
			// Erase this point & get the next one
			bBadVertex = true;
			break;
		}

		if (nLoopCount > 2 || !mesh.getVertex(ec._vfrom).isActive())
		{
			bBadVertex = true;
			break;
		}

		// If not vertex active, recalc
		if (!mesh.getVertex(ec._vto).isActive())
				quadricCollapseCostSkel(mesh, vc, opt);
		else
			break;
	}
}

void PMesh::CalcCollapseDist()
{
	int nColl = allCollPos.size();
	allCollapseDist.clear();
	allCollapseDist.resize(nColl-1);
	minCollapseDist.clear();
	maxCollapseDist.clear();
	for (int level = 0; level < allCollPos.size()-1; level++)
	{
		double fMin = FLT_MAX, fMax = 0.0;
		double fDist = 0.0;
		for (int i = 0; i < allCollPos[0].size(); i++)
		{
			fDist = (allCollPos[level+1][i] - allCollPos[level][i]).length();
			allCollapseDist[level].push_back(fDist);
			if (fDist > fMax)
				fMax = fDist;
			if (fDist < fMin)
				fMin = fDist;
		}
		minCollapseDist.push_back(fMin);
		maxCollapseDist.push_back(fMax);
	}
}

void PMesh::Simplification(CSkeOption *opt)
{
	LONGLONG sTime = StartTime(); // 开始计时

	int nVerts = _newmesh.getNumVerts();
	int nTri = _newmesh.getNumTriangles();

	_nVisTriangles = nTri; // number of visible triangles

	// 用来记录骨架点对应的表面顶点编号
	isSkeling = true;
	allCollapseFrom.clear();
	allCollapseFrom.resize(nVerts);
	// 保存骨架点信息
	simplifiedVertexRec.clear();
	simplifiedVertexRec.reserve(1024);
	// assert each vert is active -- sanity check

#ifndef NDEBUG
	assertEveryVertActive(nVerts, nTri, _newmesh);
#endif

	// calculate all 4x4 Q matrices for each vertex 
	// if using the Quadric method
	calcShapeMatrices(_newmesh);

	// This is a set of vertex pointers, ordered by edge collapse cost.
	vertexPtrSet vertSet;
	vector<vertexPtrSet::iterator> vertSetVec(nVerts);

	// Go through, calc cost here for all vertices
	calcEdgeCollapseCostsSkel(vertSet, vertSetVec, nVerts, _newmesh, _cost, opt);

	// For all vertices:
	//	find lowest cost
	//	store the edge collapse structure
	//	update all verts, triangles affected by the edge collapse
	buildEdgeCollapseListSkel(_newmesh, QUADRIC, _edgeCollList,
		vertSet, vertSetVec, opt);

	// set iterator to point to beginning
	//_edgeCollapseIter = _edgeCollList.begin();

	// 输出用时 [7/14/2011 Han Honglei]
	if (outputFile != NULL)
	{
		double elapseT = GetElapseTime(sTime);
		fprintf(outputFile, "获取骨架点用时为：%lf秒\n", elapseT);
	}

	for (int i = 0; i < nVerts; i++)
	{
		set<int>::iterator fromIter = allCollapseFrom[i].begin();
		if (allCollapseFrom[i].size() == 1 && *fromIter == -1)	// 被删除的点
			continue;

		VertexRecord vR;
		vR.center = false;
		vR.vIndex = i;
		vR.pos = _newmesh.getVertex(vR.vIndex).getXYZ();
		vR.collapseFrom = allCollapseFrom[i];

		while(fromIter != allCollapseFrom[i].end())
		{
			if (*fromIter < 0 || *fromIter > nVerts-1)
			{
				fromIter++;
				continue;
			}
			vertex &v = _newmesh.getVertex(*fromIter);
			vertex &v2 = _mesh->getVertex(*fromIter);
			v.skeletonNode = i;
			v2.skeletonNode = i;
			fromIter++;
		}
		vertex &v = _newmesh.getVertex(vR.vIndex);
		vertex &v2 = _mesh->getVertex(vR.vIndex);
		v.skeletonNode = i;
		v2.skeletonNode = i;

		simplifiedVertexRec.push_back(vR);
	}	
	simplifiedVertexNeighbour.clear();
	simplifiedVertexNeighbour.resize(simplifiedVertexRec.size());
	// 保存邻接信息，以便绘制骨架线 [7/14/2011 Han Honglei]
	for (int i = 0; i < simplifiedVertexRec.size(); i++)
	{
		int nIndex = simplifiedVertexRec[i].vIndex;
		vertex &v = _newmesh.getVertex(nIndex);
		const set<int> &neigV = v.getVertNeighbors();
		set<int>::iterator neigIter = neigV.begin();
		while (neigIter != neigV.end())
		{
			for (int j = 0; j < simplifiedVertexRec.size(); j++)
			{
				if (*neigIter == simplifiedVertexRec[j].vIndex)
					simplifiedVertexNeighbour[i].insert(j);
			}
			neigIter++;
		}
	}
	bSimplified = true;
	isSkeling = false;

	currCollPos = 0;// 回归到以前的网格// 骨架数据已经获得，将网格回复到以前的状态
	RestoreMesh();
}

void PMesh::MatrixXMatrix(double m1[4][4], double m2[4][4], double result[4][4])
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			result[i][j] = 0.0;
			for (int k = 0; k < 4; k++)
			{
				result[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
}

void PMesh::MatrixXVector( double m[4][4],double v[4], double result[4])
{
	result[0] = (m[0][0] * v[0]) + (m[0][1] * v[1]) + (m[0][2] * v[2]) + (m[0][3] * v[3]);
	result[1] = (m[1][0] * v[0]) + (m[1][1] * v[1]) + (m[1][2] * v[2]) + (m[1][3] * v[3]);
	result[2] = (m[2][0] * v[0]) + (m[2][1] * v[1]) + (m[2][2] * v[2]) + (m[2][3] * v[3]);
	result[3] = (m[3][0] * v[0]) + (m[3][1] * v[1]) + (m[3][2] * v[2]) + (m[3][3] * v[3]);

}
int nSeparateColor = 12;
float mat_separate[12][4] =  {
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
	{0.0f, 1.0f, 0.0f, 0.5},
	{0.0f, 0.0f, 1.0f, 0.5}
};	
void PMesh::DrawSimplifiedVertices()
{
	// 绘制骨架点
	glPointSize(7);	
	glBegin(GL_POINTS);
	for (int i = 0; i < simplifiedVertexRec.size(); i++)
	{
		int colorIndex = simplifiedVertexRec[i].vIndex % nSeparateColor;
		glColor4fv(mat_separate[colorIndex]);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_separate[colorIndex]);
		glVertex3d(simplifiedVertexRec[i].pos.x,simplifiedVertexRec[i].pos.y, simplifiedVertexRec[i].pos.z);
	}
	glEnd();
	
	// 绘制骨架点之间的连线
	glColor4fv(mat_separate[11]);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_separate[11]);
	glLineWidth(2);
	glBegin(GL_LINES);
	for (int i = 0; i < simplifiedVertexRec.size(); i++)
	{
		set<int>::iterator ite = simplifiedVertexNeighbour[i].begin();
		while (ite != simplifiedVertexNeighbour[i].end()) 
		{
			glVertex3d(simplifiedVertexRec[i].pos.x,simplifiedVertexRec[i].pos.y, simplifiedVertexRec[i].pos.z);
			glVertex3d(simplifiedVertexRec[*ite].pos.x,simplifiedVertexRec[*ite].pos.y, simplifiedVertexRec[*ite].pos.z);
			ite++;
		}
	}
	glEnd();
}

int PMesh::DrawOriginalMesh(DrawMeshType drawType, bool bSmooth)
{
	int colorIndex;
	int activeTriNum = 0;

	for (int i=0; i<_newmesh.getNumTriangles(); i++)
	{
		triangle t = _newmesh.getTri(i);
		if (t.isActive())
		{
			activeTriNum++;
			glLoadName(i);										// Assign Object A Name (ID)

			if (t.isSel())
			{
				glPushAttrib(GL_LIGHTING_BIT);
				glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_separate[9]);
			}

			glBegin(GL_TRIANGLES);

			const float* pFltArrayN = t.getNormal();
			glNormal3fv(pFltArrayN);

			float col[4] = {0.5, 0.5, 0.5, 0.5};
			vertex& v1 = t.getVert1vertex();
			// 测试 [7/27/2011 Han Honglei]
			glColor3f(v1._blur, 0, 0);
			switch(drawType)
			{
			case SKEL_MAP:
				if (bSimplified)	// 如果已经得到骨架
				{
					colorIndex = v1.skeletonNode%nSeparateColor;
					glColor4fv(mat_separate[colorIndex]);
					glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_separate[colorIndex]);
				}
				break;
			case COLL_DIST:
				if (bSimplified && currCollPos < allCollPos.size()-1)
				{
					float c = allCollapseDist[currCollPos][v1.getIndex()] / maxCollapseDist[currCollPos];
					col[0] = c;
					col[1] = abs(sin(3.141592654f * c));
					col[2] = 1 - c;
					glColor4fv(col);
					glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, col);
				}
				break;
			case SEGMENT:
				if (HasSegment())
				{
					colorIndex = v1.getLayer()%nSeparateColor;
					glColor4fv(mat_separate[colorIndex]);
					glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_separate[colorIndex]);
				}
				break;
			default:
				break;
			}
			if (bSmooth)
			{
				pFltArrayN = v1.getArrayVertNorms();
				glNormal3fv(pFltArrayN);
			}
			const float* pFltArray1 = v1.getArrayVerts();
			glVertex3fv(pFltArray1);

			vertex& v2 = t.getVert2vertex();
			// 测试 [7/27/2011 Han Honglei]
			glColor3f(v1._blur, 0, 0);

			switch(drawType)
			{
			case SKEL_MAP:
				if (bSimplified)	// 如果已经得到骨架
				{
					colorIndex = v2.skeletonNode%nSeparateColor;
					glColor4fv(mat_separate[colorIndex]);
					glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_separate[colorIndex]);
				}
				break;
			case COLL_DIST:
				if (bSimplified&& currCollPos < allCollPos.size()-1)
				{
					float c = allCollapseDist[currCollPos][v2.getIndex()] / maxCollapseDist[currCollPos];
					col[0] = c;
					col[1] = abs(sin(3.141592654f * c));
					col[2] = 1 - c;
					glColor4fv(col);
					glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, col);
				}		
				break;
			case SEGMENT:
				if (HasSegment())
				{
					colorIndex = v2.getLayer()%nSeparateColor;
					glColor4fv(mat_separate[colorIndex]);
					glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_separate[colorIndex]);
				}			
				break;
			default:
				break;
			}
			if (bSmooth)
			{
				pFltArrayN = v2.getArrayVertNorms();
				glNormal3fv(pFltArrayN);
			}
			const float* pFltArray2 = v2.getArrayVerts();
			glVertex3fv(pFltArray2);

			vertex& v3 = t.getVert3vertex();
			// 测试 [7/27/2011 Han Honglei]
			glColor3f(v3._blur, 0, 0);

			switch(drawType)
			{
			case SKEL_MAP:
				if (bSimplified)	// 如果已经得到骨架
				{
					colorIndex = v3.skeletonNode%nSeparateColor;
					glColor4fv(mat_separate[colorIndex]);
					glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_separate[colorIndex]);
				}
				break;
			case COLL_DIST:
				if (bSimplified&& currCollPos < allCollPos.size()-1)
				{
					float c = allCollapseDist[currCollPos][v3.getIndex()] / maxCollapseDist[currCollPos];
					col[0] = c;
					col[1] = abs(sin(3.141592654f * c));
					col[2] = 1 - c;
					glColor4fv(col);
					glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, col);
				}
				break;
			case SEGMENT:
				if (HasSegment())
				{
					colorIndex = v3.getLayer()%nSeparateColor;
					glColor4fv(mat_separate[colorIndex]);
					glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_separate[colorIndex]);
				}
				break;
			default:
				break;
			}
			if (bSmooth)
			{
				pFltArrayN = v3.getArrayVertNorms();
				glNormal3fv(pFltArrayN);
			}
			const float* pFltArray3 = v3.getArrayVerts();
			glVertex3fv(pFltArray3);

			glEnd();
			if (t.isSel())
				glPopAttrib();
		}
	}
	return activeTriNum;
}

/*
CCSMatrix PMesh::MultiplyATA(CCSMatrix A)
{
	int[] last = new int[A.RowSize];
	int[] next = new int[A.NumNonZero];
	int[] colIndex = new int[A.NumNonZero];
	for (int i = 0; i < last.Length; i++)
	{
		last[i] = -1;
	}
	for (int i = 0; i < next.Length; i++)
	{
		next[i] = -1;
	}
	for (int i = 0; i < A.ColumnSize; i++)
	{
		for (int j = A.ColIndex[i]; j < A.ColIndex[i + 1]; j++)
		{
			int k = A.RowIndex[j];
			if (last[k] != -1)
			{
				next[last[k]] = j;
			}
			last[k] = j;
			colIndex[j] = i;
		}
	}
	last = NULL;
	CCSMatrix ATA = new CCSMatrix(A.ColumnSize, A.ColumnSize);
	Set<int> set = new Set<int>();
	double[] tmp = new double[A.ColumnSize];
	List<int> ATA_RowIndex = new List<int>();
	List<double> ATA_Value = new List<double>();
	for (int i = 0; i < A.ColumnSize; i++)
	{
		tmp[i] = 0.0;
	}
	for (int j = 0; j < A.ColumnSize; j++)
	{
		for (int col = A.ColIndex[j]; col < A.ColIndex[j + 1]; col++)
		{
			int num1 = A.RowIndex[col];
			double val = A.Values[col];
			int curr = col;
			while (true)
			{
				int i = colIndex[curr];
				set.Add(i);
				tmp[i] += val * A.Values[curr];
				if (next[curr] == -1)
				{
					break;
				}
				curr = next[curr];
			}
		}
		int[] s = set.ToArray();
		Array.Sort<int>(s);
		int count = 0;
		foreach (int k in s)
		{
			if (tmp[k] != 0.0)
			{
				ATA_RowIndex.Add(k);
				ATA_Value.Add(tmp[k]);
				tmp[k] = 0.0;
				count++;
			}
		}
		ATA.ColIndex[j + 1] = ATA.ColIndex[j] + count;
		set.Clear();
	}
	ATA.RowIndex = ATA_RowIndex.ToArray();
	ATA_RowIndex = NULL;
	ATA.Values = ATA_Value.ToArray();
	ATA_Value = NULL;
	return ATA;
}


void* PMesh::SymbolicFactorization(CCSMatrix C)
{
	fixed (int* ri = C.RowIndex)
	{
		fixed (int* ci = C.ColIndex)
		{
			fixed (double* val = C.Values)
			{
				return CreaterSymbolicSolver(C.ColumnSize, C.NumNonZero, ri, ci, val);   // 这个是taucs库函数
			}
		}
	}
}

int PMesh::NumericFactorization(void* symoblicSolver, CCSMatrix C)
{
	fixed (int* ri = C.RowIndex)
	{
		fixed (int* ci = C.ColIndex)
		{
			fixed (double* val = C.Values)
			{
				return NumericFactor(symoblicSolver, C.ColumnSize, C.NumNonZero, ri, ci, val);// 这个是taucs库函数
			}
		}
	}
}

void* PMesh::Factorization(CCSMatrix C)
{
	fixed (int* ri = C.RowIndex)
	{
		fixed (int* ci = C.ColIndex)
		{
			fixed (double* val = C.Values)
			{
				return CreaterCholeskySolver(C.ColumnSize, C.NumNonZero, ri, ci, val);// 这个是taucs库函数
			}
		}
	}
}

void PMesh::ImplicitSmooth()
{
	int n = mesh.VertexCount;
	double[] x = new double[n];
	double[] b = new double[n * 3];
	double[] ATb = new double[n];
	double[] oldPos = (double[]) mesh.VertexPos.Clone();
	for (int i = 0; i < 3; i++)
	{
		int j = 0;
		for (int k = 0; j < n; k += 3)
		{
			b[j] = 0.0;
			b[j + n] = mesh.VertexPos[k + i] * posWeight[j];
			b[(j + n) + n] = 0.0;
			j++;
		}
		ccsA.PreMultiply(b, ATb);
		if (opt->UseSymbolicSolver)
		{
			fixed (double* _x = x)
			{
				fixed (double* _ATb = ATb)
				{
					NumericSolve(symbolicSolver, _x, _ATb);
				}
			}
		}
		else
		{
			double[] CS$0$0002;
			if (((CS$0$0002 = x) == NULL) || (CS$0$0002.Length == 0))
			{
				_x = NULL;
				goto Label_0124;
			}
			fixed (double* _x = CS$0$0002)
			{
				double[] CS$0$0003;
Label_0124:
				if (((CS$0$0003 = ATb) == NULL) || (CS$0$0003.Length == 0))
				{
					_ATb = NULL;
					goto Label_0140;
				}
				fixed (double* _ATb = CS$0$0003)
				{
Label_0140:
					Solve(solver, _x, _ATb);
				}
			}
		}
		lock (mesh.VertexPos)
		{
			int j = 0;
			for (int k = 0; j < n; k += 3)
			{
				mesh.VertexPos[k + i] = x[j];
				j++;
			}
		}
	}
	iter++;
	mesh.ComputeFaceNormal();
	mesh.ComputeVertexNormal();
	int i = 0;
	for (int j = 0; i < n; j += 3)
	{
		double d1 = mesh.VertexPos[j] - oldPos[j];
		double d2 = mesh.VertexPos[j + 1] - oldPos[j + 1];
		double d3 = mesh.VertexPos[j + 2] - oldPos[j + 2];
		collapsedLength[i] += Math.Sqrt(((d1 * d1) + (d2 * d2)) + (d3 * d3));
		i++;
	}
}
*/
