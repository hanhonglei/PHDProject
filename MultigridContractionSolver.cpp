#include "MultigridContractionSolver.h"
#include <math.h>
#include <float.h>
//extern "C" {
//#include <taucs.h>
//}

MultigridContractionSolver::MultigridContractionSolver(void)
{
	convergeRatio = 0.001;
	mesh = NULL;
	//resolutions = new List<Resolution>();
	targetVertexCount = 0x7d0;
	// 利用动态链接库使用函数 [6/26/2011 Han Honglei]
	hDllInst = NULL;

}

MultigridContractionSolver::~MultigridContractionSolver(void)
{
	//if (hDllInst != NULL)
	//{
	//	FreeLibrary(hDllInst);
	//	hDllInst = NULL;
	//}
}



MultigridContractionSolver::MultigridContractionSolver(	jmsMesh* _newmesh)
{
	// 利用动态链接库使用函数 [6/26/2011 Han Honglei]
	hDllInst = NULL;

	convergeRatio = 0.001;
	mesh = NULL;
	//resolutions = new List<Resolution>();
	targetVertexCount = 0x7d0;

	//Program.PrintText("[Create Multigrid Solver]");
	mesh = _newmesh;
	Resolution res(mesh);
	resolutions.reserve(1024);
	resolutions.push_back(res);
	//resolutions.Add(new Resolution(mesh));
	while (resolutions[resolutions.size() - 1].vertexList.size() > targetVertexCount)
	{
		//Program.PrintText(r.vertexList.Length.ToString());
		//resolutions.Add(GetNextLevel(r));
		Resolution res;
		resolutions.push_back(res);
		GetNextLevel(resolutions[resolutions.size() - 2], resolutions[resolutions.size() - 1]);
	}
	//Program.PrintText(resolutions[resolutions.Count - 1].vertexList.Length.ToString());
	//resolutions.Reverse();
	// 由于反向的速度很慢，所以，不进行反向处理 [7/6/2011 Han Honglei]
	//reverse(resolutions.begin(),resolutions.end());
}

//[DllImport("taucs.dll")]
//private static extern unsafe void* CreaterCholeskySolver(int n, int nnz, int* rowIndex, int* colIndex, double* value);


void* MultigridContractionSolver::Factorization(CCSMatrix &C)
{
	void* p = NULL;
	
	// 利用动态链接库使用函数 [6/26/2011 Han Honglei]
	if (hDllInst == NULL)
	{
		hDllInst = LoadLibrary(_T("taucs.DLL"));
	}
	if(hDllInst)
	{
		typedef void* (* MYFUNC)(int, int , int* , int* , double*);
		MYFUNC CreaterCholeskySolver = NULL; 
		CreaterCholeskySolver = (MYFUNC)GetProcAddress
			(hDllInst,"CreaterCholeskySolver");
		if(CreaterCholeskySolver)
		{
			p = CreaterCholeskySolver(C.ColumnSize(), C.NumNonZero(), C.rowIndex, C.colIndex, C.values/*C.ColumnSize(), C.NumNonZero(), C.RowIndex(), C.ColIndex(), C.Values()*/);
		}
		//FreeLibrary(hDllInst);
	}

	return p;
}


//[DllImport("taucs.dll")]
//private static extern unsafe int Solve(void* solver, double* x, double* b);
//[DllImport("taucs.dll")]
//private static extern unsafe double SolveEx(void* solver, double* x, int xIndex, double* b, int bIndex);
double** MultigridContractionSolver::SolveSystem(double* lapWeight, double* posWeight)
{
	int vn = mesh->getNumVerts();
	int faceCount = mesh->getNumTriangles();
	double** pos = new double* [3];
	for (int i = 0; i < 3; i++)
	{
		pos[i] = new double[vn];
	}
	int i = 0;
	for (int j = 0; i < vn; j += 3)
	{
		pos[0][i] = mesh->getVertex(i).getXYZ().x;
		pos[1][i] = mesh->getVertex(i).getXYZ().y;
		pos[2][i] = mesh->getVertex(i).getXYZ().z;
		i++;
	}
	for (int levelR = 0; levelR < resolutions.size()/*.Count*/; levelR++)
	{
		int level = resolutions.size() - levelR - 1;
		ColMatrix colA = BuildMatrixA(resolutions[level], lapWeight, posWeight);
		CCSMatrix* ccsATA = MultiplyATA(colA);
		//string s = "#: " + r.vertexList.Length + " iter: ";
		int n = resolutions[level].vertexList.size();
		double* x = new double[n];
		double* b = new double[n * 2];
		double* ATb = new double[n];
		void* solver = NULL;
		if (level == 0)
		{
			solver = Factorization(*ccsATA);
		}
		delete ccsATA;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < n; j++)
			{
				b[j] = 0.0;
				int k = resolutions[level].vertexList[j];
				double tmp[3] = {mesh->getVertex(k).getXYZ().x, mesh->getVertex(k).getXYZ().y, mesh->getVertex(k).getXYZ().z};
				b[j + n] = /*mesh->VertexPos[(k * 3) + i]*/tmp[i] * posWeight[k];
				x[j] = pos[i][k];
			}
			colA.PreMultiply(b, ATb);
			if (level == 0)
			{
				// 利用动态链接库使用函数 [6/26/2011 Han Honglei]
				if (hDllInst == NULL)
				{
					hDllInst = LoadLibrary("taucs.DLL");
				}
				if(hDllInst)
				{
					typedef int (* MYFUNC)(void* , double* , double* );
					MYFUNC Solve = NULL; 
					Solve = (MYFUNC)GetProcAddress
						(hDllInst,"Solve");

					if(Solve)
					{
						int r = Solve(solver, x, ATb);
					}
					//FreeLibrary(hDllInst);
				}
			}
			else
			{
				int iter = colA.ATACG(x, ATb, convergeRatio, n);
				//s = s + " " + iter;
			}
			for (int j = 0; j < n; j++)
			{
				pos[i][resolutions[level].vertexList[j]] = x[j];
			}
		}
		//Program.PrintText(s);
		if (level < (resolutions.size()/*.Count*/ - 1))
		{
			for (int i = 0; i < resolutions[level].collapsedIndex.size(); i++)
			{
				int index = resolutions[level].collapsedIndex[i];
				//Vector3d p = new Vector3d();
				vec3 p;
				double totWeight = 0.0;
				map<int,double>::iterator it;
				for ( it=resolutions[level].weight[i].begin() ; it != resolutions[level].weight[i].end(); it++ )
				{
					double w = (*it).second;
					p.x += pos[0][(*it).first] * w;
					p.y += pos[1][(*it).first] * w;
					p.z += pos[2][(*it).first] * w;
					totWeight += w;

				}
				//foreach (int adj in r.weight[i].Keys)
				//{
				//	double w = r.weight[i][adj];
				//	p.x += pos[0][adj] * w;
				//	p.y += pos[1][adj] * w;
				//	p.z += pos[2][adj] * w;
				//	totWeight += w;
				//}
				pos[0][index] = p.x /= totWeight;
				pos[1][index] = p.y /= totWeight;
				pos[2][index] = p.z /= totWeight;
			}
		}
		// 删除动态分配的空间 [6/28/2011 Han Honglei]
		delete[] x;
		delete[] b;
		delete[] ATb;
		if (hDllInst)
		{
			typedef int (* MYFUNC)(void*);
			MYFUNC FreeSolver = NULL; 
			FreeSolver = (MYFUNC)GetProcAddress
				(hDllInst,"FreeSolver");
			if(FreeSolver)
			{
				if (solver != NULL)
				{
					FreeSolver(solver);	// 这个是taucs库函数
					solver = NULL;
				}
			}
		}
	}
	return pos;
}

ColMatrix MultigridContractionSolver::BuildMatrixA(Resolution &r, double* lapWeight, double* posWeight)
{
	int* maps = new int[mesh->getNumVerts()];
	for (int i = 0; i < r.vertexList.size(); i++)
	{
		maps[r.vertexList[i]] = i;
	}
	int n = r.vertexList.size();
	int fn = r.faceList.size() / 3;
	SparseMatrix A(2 * n, n);
	int i = 0;
	for (int j = 0; i < fn; j += 3)
	{
		int c1 = r.faceList[j];
		int c2 = r.faceList[j + 1];
		int c3 = r.faceList[j + 2];
		//Vector3d v1 = new Vector3d(mesh->VertexPos, c1 * 3);
		//Vector3d v2 = new Vector3d(mesh->VertexPos, c2 * 3);
		//Vector3d v3 = new Vector3d(mesh->VertexPos, c3 * 3);
		vec3 v1 = vec3(mesh->getVertex(c1).getXYZ());
		vec3 v2 = vec3(mesh->getVertex(c2).getXYZ());
		vec3 v3 = vec3(mesh->getVertex(c3).getXYZ());

		vec3 CS0 = v2 - v1;
		vec3 CS1 = v2 - v1;
		double cot1 = CS0.dot(v3 - v1) / CS1.cross(v3 - v1).length();
		vec3 CS3 = v3 - v2;
		vec3 CS4 = v3 - v2;
		double cot2 = CS3.dot(v1 - v2) / CS4.cross(v1 - v2).length();
		vec3 CS6 = v1 - v3;
		vec3 CS7 = v1 - v3;
		double cot3 = CS6.dot(v2 - v3) / CS7.cross(v2 - v3).length();
		//Vector3d CS0000 = v2 - v1;
		//Vector3d CS0001 = v2 - v1;
		//double cot1 = CS0000.Dot(v3 - v1) / CS0001.Cross(v3 - v1).Length();
		//Vector3d CS0003 = v3 - v2;
		//Vector3d CS0004 = v3 - v2;
		//double cot2 = CS0003.Dot(v1 - v2) / CS0004.Cross(v1 - v2).Length();
		//Vector3d CS0006 = v1 - v3;
		//Vector3d CS0007 = v1 - v3;
		//double cot3 = CS0006.Dot(v2 - v3) / CS0007.Cross(v2 - v3).Length();

		//if (double.IsNaN(cot1))
		//{
		//	throw new Exception();
		//}
		//if (double.IsNaN(cot2))
		//{
		//	throw new Exception();
		//}
		//if (double.IsNaN(cot3))
		//{
		//	throw new Exception();
		//}
		c1 = maps[c1];
		c2 = maps[c2];
		c3 = maps[c3];
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
	vector<vector<Element>>::iterator itARow = A.Rows();
	for (int i = 0; i < n; i++)
	{
		double tot = 0.0;

		for (vector<Element>::iterator it = itARow[i].begin(); it < itARow[i].end(); it++)
		{
			if(it->i != it->j)
				tot += it->value;				
		}
		//foreach (SparseMatrix.Element e in A.Rows[i])
		//{
		//	if (e.i != e.j)
		//	{
		//		tot += e.value;
		//	}
		//}
		if (tot > 10000.0)
		{
			//foreach (SparseMatrix.Element e in A.Rows[i])
			//{
			//	e.value /= tot / 10000.0;
			//}
			for (vector<Element>::iterator it = itARow[i].begin(); it < itARow[i].end(); it++)
			{
				it->value /= tot / 10000.0;
			}

		}
		//foreach (SparseMatrix.Element e in A.Rows[i])
		//{
		//	e.value *= lapWeight[r.vertexList[i]];
		//}
		for (vector<Element>::iterator it = itARow[i].begin(); it < itARow[i].end(); it++)
			it->value *= lapWeight[r.vertexList[i]];
	}
	for (int i = 0; i < n; i++)
	{
		A.AddElement(i + n, i, posWeight[r.vertexList[i]]);
	}
	A.SortElement();
	ColMatrixCreator cmA(2 * n, n);
	for (int i = 0; i < A.RowSize(); i++)
	{
		for (vector<Element>::iterator it = itARow[i].begin(); it < itARow[i].end(); it++)
			cmA.AddValueTo(it->i, it->j, it->value);
	}
	//foreach (List<SparseMatrix.Element> row in A.Rows)
	//{
	//	foreach (SparseMatrix.Element e in row)
	//	{
	//		cmA.AddValueTo(e.i, e.j, e.value);
	//	}
	//}
	// delete [7/3/2011 Han Honglei]
	delete []maps;

	return cmA.ToColMatrix();
}


//[DllImport("taucs.dll")]
//private static extern unsafe int FreeSolver(void* sp);
void MultigridContractionSolver::GetNextLevel(Resolution &r, Resolution &result)
{
	//Vector3d CS0002;
	vec3 CS0002;
	int vn = mesh->getNumVerts();
	int* oldMap = new int[vn];
	for (int i = 0; i < vn; i++)
	{
		oldMap[i] = -1;
	}
	for (int i = 0; i < r.vertexList.size(); i++)
	{
		oldMap[r.vertexList[i]] = i;
	}
	bool* marked = new bool[vn];
	for (int i = 0; i < vn; i++)
	{
		marked[i] = false;
	}
// 修改为不进行拷贝的方式 [7/6/2011 Han Honglei]
	//Resolution res;
	//res.vertexList = remainingIndex;
	//res.faceList = newFaceIndex;
	//res.adjVertex = newAdjVertexIndex;
	//res.adjFace = newAdjFaceIndex;
	//res.collapsedIndex = collapsedIndex;
	//res.weight = interpolateWeight ;

	//vector<int> collapsedIndex;
	//vector<int> remainingIndex;
	for (int i = 0; i < r.vertexList.size(); i++)
	{
		int index = r.vertexList[i];
		bool found = false;
		for (set<int>::iterator it = r.adjVertex[i].begin(); it != r.adjVertex[i].end(); it++)
		{
			int count = 0;
			for (set<int>::iterator it2 = r.adjVertex[oldMap[*it]].begin(); it2 != r.adjVertex[oldMap[*it]].end(); it2++)
			{
				for (set<int>::iterator it3 = r.adjVertex[oldMap[*it2]].begin(); it3 != r.adjVertex[oldMap[*it2]].end(); it3++)
				{
					if (*it3 == index)
					{
						count++;
					}
				}
			}
			if (count == 2)
			{
				found = true;
				break;
			}
		}
		//foreach (int adj in r.adjVertex[i])
		//{
		//	int count = 0;
		//	foreach (int adj2 in r.adjVertex[oldMap[adj]])
		//	{
		//		if (r.adjVertex[oldMap[adj2]].Contains(index))
		//		{W
		//			count++;
		//		}
		//	}
		//	if (count == 2)
		//	{
		//		found = true;
		//		break;
		//	}
		//}
		if (!found)
		{
			marked[index] = true;
		}
		if (marked[index])
		{
			//remainingIndex.Add(index);
			result.vertexList.push_back(index);
		}
		else
		{
			result.collapsedIndex.push_back(index);
			for (set<int>::iterator it = r.adjVertex[i].begin(); it != r.adjVertex[i].end(); it++)
			{
				marked[*it] = true;
			}
			//collapsedIndex.Add(index);
			//foreach (int adj in r.adjVertex[i])
			//{
			//	marked[adj] = true;
			//}
		}
	}
	int* collapseTo = new int[vn];
	for (int i = 0; i < vn; i++)
	{
		collapseTo[i] = i;
	}
	for (int i = 0; i < result.collapsedIndex.size()/*Count*/; i++)
	{
		int index = result.collapsedIndex[i];
		vec3 p(mesh->getVertex(index).getXYZ());
		double minLen = DBL_MAX;
		//Vector3d p = new Vector3d(mesh->VertexPos, index * 3);
		//double minLen = double.MaxValue;
		int minAdj = -1;
		for (set<int>::iterator it = r.adjVertex[oldMap[index]].begin(); it != r.adjVertex[oldMap[index]].end(); it++)
		{
			int count = 0;
			for (set<int>::iterator it2 = r.adjVertex[oldMap[*it]].begin(); it2 != r.adjVertex[oldMap[*it]].end(); it2++)
			{
				for (set<int>::iterator it3 = r.adjVertex[oldMap[*it2]].begin(); it3 != r.adjVertex[oldMap[*it2]].end(); it3++)
				{
					if (*it3 == index)
					{
						count ++;
					}
				}
			}
				//foreach (int adj in r.adjVertex[oldMap[index]])
				//{
				//	int count = 0;
				//	foreach (int adj2 in r.adjVertex[oldMap[adj]])
				//	{
				//		if (r.adjVertex[oldMap[adj2]].Contains(index))
				//		{
				//			count++;
				//		}
				//	}
			if (count == 2)
			{
				if (/*adj*/*it == index)
				{
					//Program.PrintText("?");
					int test = 1;
				}
				else
				{
					//Vector3d q = new Vector3d(mesh->VertexPos, adj * 3);
					vec3 q(mesh->getVertex(*it).getXYZ());
					CS0002 = p - q;
					double len = CS0002.length();
					if (len < minLen)
					{
						minLen = len;
						minAdj = *it;
					}
					if (*it == index)
					{
						throw *it;
						//throw new Exception();
					}
				}
			}
		}
		if (minAdj == -1)
		{
			//throw new Exception();
			throw minAdj;
		}
		collapseTo[index] = minAdj;
	}


	//Set<int>[] tmpAdjList = new Set<int>[r.vertexList.Length];
	vector<set<int>> tmpAdjList;
	tmpAdjList.reserve(r.vertexList.size());

	for (int i = 0; i < r.vertexList.size(); i++)
	{
		//Set<int> s = new Set<int>(6);
		set<int> s;
		for (set<int>::iterator it = r.adjVertex[i].begin(); it != r.adjVertex[i].end(); it++)
		{
			s.insert(*it);
		}
		//foreach (int adj in r.adjVertex[i])
		//{
		//	s.Add(adj);
		//}
		//tmpAdjList[i] = s;
		tmpAdjList.push_back(s);
	}
	for (int tmp = 0; tmp < result.collapsedIndex.size(); tmp++)
	{
		int index = result.collapsedIndex[tmp];
	//foreach (int index in collapsedIndex)
	//{
		int i = oldMap[index];
		int to = collapseTo[index];
		int k = oldMap[to];
		for (set<int>::iterator it = r.adjVertex[i].begin(); it != r.adjVertex[i].end(); it++)
		//foreach (int adj in r.adjVertex[i])
		{
			int j = oldMap[*it/*adj*/];
			//tmpAdjList[j].Remove(index);
			tmpAdjList[j].erase(index);
			if ((*it != to) && (j != -1))
			{
				//tmpAdjList[j].Add(to);
				//tmpAdjList[k].Add(adj);
				tmpAdjList[j].insert(to);
				tmpAdjList[k].insert(*it);
			}
		}
	}
	//vector<set<int>> newAdjVertexIndex;
	//List<Set<int>> newAdjVertexIndex = new List<Set<int>>();
	for (int i = 0; i < r.vertexList.size(); i++)
	{
		int index = r.vertexList[i];
		if (marked[index])
		{
			result.adjVertex.push_back(tmpAdjList[i]);
			//newAdjVertexIndex.Add(tmpAdjList[i]);
		}
	}
	//List<int> newFaceIndex = new List<int>();
	for (int i = 0; i < r.faceList.size(); i += 3)
	{
		int c1 = r.faceList[i];
		int c2 = r.faceList[i + 1];
		int c3 = r.faceList[i + 2];
		c1 = collapseTo[c1];
		c2 = collapseTo[c2];
		c3 = collapseTo[c3];
		if (((c1 != c2) && (c2 != c3)) && (c3 != c1))
		{
			//newFaceIndex.Add(c1);
			//newFaceIndex.Add(c2);
			//newFaceIndex.Add(c3);
			result.faceList.push_back(c1);
			result.faceList.push_back(c2);
			result.faceList.push_back(c3);
		}
	}
	int* maps = new int[vn];
	for (int i = 0; i < vn; i++)
	{
		maps[i] = -1;
	}
	//vector<set<int>> newAdjFaceIndex;
	//newAdjFaceIndex.reserve(remainingIndex.size());
	result.adjFace.resize(result.vertexList.size());

	for (int i = 0; i < result.vertexList.size(); i++)
	{
		maps[result.vertexList[i]] = i;
		//newAdjFaceIndex[i] = new Set<int>();
	}
	int i = 0;
	for (int j = 0; i < result.faceList.size(); j++)
	{
		int c1 = maps[result.faceList[i]];
		int c2 = maps[result.faceList[i + 1]];
		int c3 = maps[result.faceList[i + 2]];
		if (((c1 == -1) || (c2 == -1)) || (c3 == -1))
		{
			throw -1;
		}
		//newAdjFaceIndex[c1].Add(j);
		//newAdjFaceIndex[c2].Add(j);
		//newAdjFaceIndex[c3].Add(j);
		result.adjFace[c1].insert(j);
		result.adjFace[c2].insert(j);
		result.adjFace[c3].insert(j);
		i += 3;
	}
	//vector<map<int, double>> interpolateWeight;
	result.weight.reserve(result.collapsedIndex.size());
	//Dictionary<int, double>[] interpolateWeight = new Dictionary<int, double>[collapsedIndex.Count];
	for (int i = 0; i < result.collapsedIndex.size(); i++)
	{
		//Dictionary<int, double> dict = new Dictionary<int, double>(6);
		map<int, double> dict;
		int index = result.collapsedIndex[i];
		int j = oldMap[index];
		for (set<int>::iterator it = r.adjFace[j].begin(); it != r.adjFace[j].end(); it++)
		//foreach (int adj in r.adjFace[j])
		{
			map<int, double> CS0004;
			int CS0005;
			int k = *it/*adj*/ * 3;
			int c1 = r.faceList[k];
			int c2 = r.faceList[k + 1];
			int c3 = r.faceList[k + 2];
			if (c2 == index)
			{
				c2 = c3;
				c3 = c1;
				c1 = index;
			}
			if (c3 == index)
			{
				c3 = c2;
				c2 = c1;
				c1 = index;
			}
			//Vector3d v1 = new Vector3d(mesh->VertexPos, c1 * 3);
			//Vector3d v2 = new Vector3d(mesh->VertexPos, c2 * 3);
			//Vector3d v3 = new Vector3d(mesh->VertexPos, c3 * 3);
			vec3 v1(mesh->getVertex(c1).getXYZ());
			vec3 v2(mesh->getVertex(c2).getXYZ());
			vec3 v3(mesh->getVertex(c3).getXYZ());
			CS0002 = v2 - v1;
			CS0002 = v2 - v1;
			double cot1 = CS0002.dot(v3 - v1) / CS0002.cross(v3 - v1).length();
			CS0002 = v3 - v2;
			CS0002 = v3 - v2;
			double cot2 = CS0002.dot(v1 - v2) / CS0002.cross(v1 - v2).length();
			CS0002 = v1 - v3;
			CS0002 = v1 - v3;
			double cot3 = CS0002.dot(v2 - v3) / CS0002.cross(v2 - v3).length();
			
			if (_isnan(cot1))
			{
				throw FALSE;
			}
			if (_isnan(cot2))
			{
				throw FALSE;
			}
			if (_isnan(cot3))
			{
				throw FALSE;
			}
			if (/*dict.ContainsKey(c3)*/dict.find(c3) != dict.end())
			{
				(CS0004 = dict)[CS0005 = c3] = CS0004[CS0005] + cot2;
			}
			else
			{
				//dict.Add(c3, cot2);
				dict[c3] = cot2;

			}
			if (/*dict.ContainsKey(c2)*/dict.find(c2) != dict.end())
			{
				(CS0004 = dict)[CS0005 = c2] = CS0004[CS0005] + cot3;
			}
			else
			{
				//dict.Add(c2, cot3);
				dict[c2] = cot3;
			}
		}
		result.weight.push_back(dict);
	}
	// 清理工作 [7/3/2011 Han Honglei]
	delete[] maps;
	delete[] oldMap;
	delete[] marked;
	delete[] collapseTo;

	// 由于vector是拷贝的方式，所以尽量少一些相等操作 [7/6/2011 Han Honglei]
	//Resolution res;
	//res.vertexList = remainingIndex;
	//res.faceList = newFaceIndex;
	//res.adjVertex = newAdjVertexIndex;
	//res.adjFace = newAdjFaceIndex;
	//res.collapsedIndex = collapsedIndex;
	//res.weight = interpolateWeight ;
	//return res;
}

CCSMatrix* MultigridContractionSolver::MultiplyATA(ColMatrix A)
{
	// 这个函数中引起的数值有问题 [6/27/2011 Han Honglei]
	vector<int> count;
	count.resize(A.PRowSize(), 0);
	//for (int i = 0; i < count.size(); i++)
	//{
	//	//count[i] = 0;
	//}
	for (vector<vector<int>>::iterator it = A.rowIndex.begin(); it != A.rowIndex.end(); it++)
	{
		for (vector<int>::iterator iti = it->begin(); iti != it->end(); iti++)
		{
			count[*iti]++;
		}
	}
	//foreach (int[] r in A.rowIndex)
	//{
	//	foreach (int ri in r)
	//	{
	//		count[ri]++;
	//	}
	//}
	//int[][] colIndex = new int[A.RowSize][];
	//int[][] listIndex = new int[A.RowSize][];
	vector<vector<int>> colIndex;
	colIndex.resize(A.PRowSize());
	vector<vector<int>> listIndex;
	listIndex.resize(A.PRowSize());
	for (int i = 0; i < A.PRowSize(); i++)
	{
		colIndex[i].resize(count[i]) ;
		listIndex[i].resize(count[i]);
	}
	for (int i = 0; i < count.size(); i++)
	{
		count[i] = 0;
	}
	for (int i = 0; i < A.values.size(); i++)
	{
		///*int[]*/vector<int> row = A.rowIndex[i];
		for (int j = 0; j < A.rowIndex[i].size(); j++)
		{
			int r = A.rowIndex[i][j];
			int c = count[r];
			colIndex[r][c] = i;
			listIndex[r][c] = j;
			count[r]++;
		}
	}
	//count = null;
	count.clear();

	CCSMatrix *ATA = new CCSMatrix(A.PColumnSize(), A.PColumnSize());
	//Set<int> sets = new Set<int>();
	set<int> sets;
	double* tmp = new double[A.PColumnSize()];
	
	//List<int> ATA_RowIndex = new List<int>();
	//List<double> ATA_Value = new List<double>();
	vector<int> ATA_RowIndex;
	vector<double> ATA_Value;
	for (int i = 0; i < A.PColumnSize(); i++)
	{
		tmp[i] = 0.0;
	}
	for (int j = 0; j < A.PColumnSize(); j++)
	{
		for (int ri = 0; ri < A.rowIndex[j].size(); ri++)
		{
			int k = A.rowIndex[j][ri];
			double val = A.values[j][ri];
			for (int k2 = 0; k2 < colIndex[k].size(); k2++)
			{
				int i = colIndex[k][k2];
				if (i >= j)
				{
					//sets.Add(i);
					sets.insert(i);
					tmp[i] += val * A.values[i][listIndex[k][k2]];
				}
			}
		}
		//sort(sets.begin(), sets.end());	// sorts里面的值本身就是排序的，所以删除排序
		//int[] s = sets.ToArray();
		//Array.Sort<int>(s);
		int cc = 0;
		for (set<int>::iterator it = sets.begin(); it != sets.end(); it++)
		{
			if (tmp[*it] != 0.0)
			{
				ATA_RowIndex.push_back(*it);
				ATA_Value.push_back(tmp[*it]);
				tmp[*it] = 0.0;
				cc++;
			}
		}
		//foreach (int k in s)
		//{
		//	if (tmp[k] != 0.0)
		//	{
		//		ATA_RowIndex.Add(k);
		//		ATA_Value.Add(tmp[k]);
		//		tmp[k] = 0.0;
		//		cc++;
		//	}
		//}
		ATA->ColIndex()[j + 1] = ATA->ColIndex()[j] + cc;
		sets.clear();
	}
	//ATA->RowIndex = ATA_RowIndex.ToArray();
	//ATA->Values = ATA_Value.ToArray();
	if (ATA->rowIndex == NULL)
		ATA->rowIndex = new int [ATA_RowIndex.size()];
	if (ATA->values == NULL)
		ATA->values = new double [ATA_Value.size()];

	for (int i = 0; i < ATA_RowIndex.size(); i++)
	{
		ATA->RowIndex()[i] = ATA_RowIndex[i];
		ATA->Values()[i] = ATA_Value[i];
	}
	ATA->nValuesNum = ATA_Value.size();
	//int * p;
	//p=ATA_RowIndex.get_allocator().allocate(ATA_RowIndex.size());
	//memcpy(ATA->RowIndex(), p, ATA_RowIndex.size()*sizeof(int));
	//ATA_RowIndex.get_allocator().deallocate(p, ATA_RowIndex.size());
	//ATA_RowIndex.clear()/* = null*/;

	//double *q;
	//q=ATA_Value.get_allocator().allocate(ATA_Value.size());
	//memcpy(ATA->Values(), q, ATA_Value.size()*sizeof(double));
	//ATA_Value.get_allocator().deallocate(q,ATA_Value.size());
	//ATA_Value.clear()/* = null*/;
	// clear [7/3/2011 Han Honglei]
	delete []tmp;

	return ATA;
}

//////////////////////////////////////////////////////////////////////////
// ColMatrix中的成员函数实现

ColMatrix::ColMatrix(int _m, int _n)
{
	m = _m;
	n = _n;
	rowIndex.resize(n);
	values.resize(n);
	//for (int i = 0; i < n; i++)
	//{
	//	rowIndex[i] = null;
	//	values[i] = null;
	//}
}

int ColMatrix::ATACG(double* x, double* b, double eps, int maxIter)
{
	double* inv = new double[n];
	double* r = new double[n];
	double* d = new double[n];
	double* q = new double[n];
	double* s = new double[n];
	double errNew = 0.0;
	double err = 0.0;
	double errOld = 0.0;
	double* tmp = new double[m];
	for (int i = 0; i < n; i++)
	{
		inv[i] = 0.0;
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < rowIndex[i].size(); j++)
		{
			double val = values[i][j];
			inv[i] += val * val;
		}
	}
	for (int i = 0; i < n; i++)
	{
		inv[i] = 1.0 / inv[i];
	}
	int iter = 0;
	Multiply(x, tmp);
	PreMultiply(tmp, r);
	for (int i = 0; i < n; i++)
	{
		r[i] = b[i] - r[i];
	}
	for (int i = 0; i < n; i++)
	{
		d[i] = inv[i] * r[i];
	}
	for (int i = 0; i < n; i++)
	{
		errNew += r[i] * d[i];
	}
	err = errNew;
	//Program.PrintText("err: " + err.ToString());
	while ((iter < maxIter) && (errNew > (eps * err)))
	{
		Multiply(d, tmp);
		PreMultiply(tmp, q);
		double alpha = 0.0;
		for (int i = 0; i < n; i++)
		{
			alpha += d[i] * q[i];
		}
		alpha = errNew / alpha;
		for (int i = 0; i < n; i++)
		{
			x[i] += alpha * d[i];
		}
		if ((iter % 50) == 0)
		{
			Multiply(x, tmp);
			PreMultiply(tmp, r);
			for (int i = 0; i < n; i++)
			{
				r[i] = b[i] - r[i];
			}
		}
		else
		{
			for (int i = 0; i < n; i++)
			{
				r[i] -= alpha * q[i];
			}
		}
		for (int i = 0; i < n; i++)
		{
			s[i] = inv[i] * r[i];
		}
		errOld = errNew;
		errNew = 0.0;
		for (int i = 0; i < n; i++)
		{
			errNew += r[i] * s[i];
		}
		double beta = errNew / errOld;
		for (int i = 0; i < n; i++)
		{
			d[i] = s[i] + (beta * d[i]);
		}
		iter++;
	}
// clear [7/3/2011 Han Honglei]
	delete[] inv;
	delete[] r;
	delete[] d;
	delete[] q;
	delete[] s;
	delete[] tmp;

	return iter;
}

void ColMatrix::Multiply(double* xIn, double* xOut)
{
	//if ((xIn.Length < n) || (xOut.Length < m))
	//{
	//	throw new ArgumentException();
	//}
	for (int i = 0; i < m; i++)
	{
		xOut[i] = 0.0;
	}
	for (int i = 0; i < n; i++)
	{
		//int[] r = rowIndex[i];
		//double[] v = values[i];
		for (int j = 0; j < rowIndex[i].size()/*r.Length*/; j++)
		{
			//int ri = r[j];
			//xOut[ri] += v[j] * xIn[i];
			int ri = rowIndex[i][j];
			xOut[ri] += values[i][j] * xIn[i];
		}
	}
}

int ColMatrix::NumOfElements()
{
	int c = 0;
	for (vector<vector<int>>::iterator it = rowIndex.begin(); it != rowIndex.end(); it++)
	{
		c += it->size();
	}
	//foreach (int[] r in rowIndex)
	//{
	//	c += r.Length;
	//}
	return c;
}

//public static implicit operator MultigridContractionSolver.ColMatrix(MultigridContractionSolver.ColMatrixCreator CC)
//{
//	return CC.ToColMatrix();
//}

void ColMatrix::PreMultiply(double* xIn, double* xOut)
{
	//if ((xIn.Length < m) || (xOut.Length < n))
	//{
	//	throw new ArgumentException();
	//}
	for (int i = 0; i < n; i++)
	{
		xOut[i] = 0.0;
	}
	for (int i = 0; i < n; i++)
	{
		double sum = 0.0;
		//int[] r = rowIndex[i];
		//double[] v = values[i];
		for (int j = 0; j < rowIndex[i].size()/*r.Length*/; j++)
		{
			int ri = rowIndex[i][j];
			sum += values[i][j] * xIn[ri];
		}
		xOut[i] = sum;
	}
}

int ColMatrix::PColumnSize()
{
		return n;
}

int ColMatrix::PRowSize()
{
		return m;
}

//////////////////////////////////////////////////////////////////////////
// ColMatrixCreator的成员函数实现

ColMatrixCreator::ColMatrixCreator(int _m, int _n)
{
	m = _m;
	n = _n;
	rowIndex.resize(n)/* = new List<int>[n]*/;
	values.resize(n)/* = new List<double>[n]*/;
	//for (int i = 0; i < n; i++)
	//{
	//	rowIndex[i] = new List<int>();
	//	values[i] = new List<double>();
	//}
}

void ColMatrixCreator::AddValueTo(int i, int j, double value)
{
	//List<int> r = rowIndex[j];
	if (i >= m)
	{
		m = i + 1;
	}
	int ri = -1;
	for (int k = 0; k < rowIndex[j].size()/*r.Count*/; k++)
	{
		if (rowIndex[j][k] == i)
		{
			ri = k;
			break;
		}
	}
	if (ri == -1)
	{
		rowIndex[j].push_back(i);
		values[j].push_back(value);
	}
	else
	{
		//List<double> CS0000;
		//int CS0001;
		//(CS0000 = values[j])[CS0001 = ri] = CS0000[CS0001] + value;
		values[j][ri] += value;
	}
}

ColMatrix ColMatrixCreator::ToColMatrix()
{
	ColMatrix C(m, n);
	for (int i = 0; i < n; i++)
	{
		//try
		//{
			C.rowIndex[i] = rowIndex[i]/*.ToArray()*/;
		//}
		//catch (OutOfMemoryException)
		//{
		//	GC.Collect();
		//	C.rowIndex[i] = rowIndex[i].ToArray();
		//}
		//rowIndex[i] = null;
		//try
		//{
			C.values[i] = values[i]/*.ToArray()*/;
		//}
		//catch (OutOfMemoryException)
		//{
		//	GC.Collect();
		//	C.values[i] = values[i].ToArray();
		//}
		//values[i] = null;
	}
	return C;
}

int ColMatrixCreator::ColumnSize()
{
		return n;
}

int ColMatrixCreator::RowSize()
{
		return m;
}
//
//public void PreMultiply(double[] xIn, double[] xOut, int[] index)
//{
//	if ((xIn.Length < m) || (xOut.Length < n))
//	{
//		throw new ArgumentException();
//	}
//	for (int i = 0; i < n; i++)
//	{
//		xOut[i] = 0.0;
//	}
//	foreach (int i in index)
//	{
//		double sum = 0.0;
//		int[] r = rowIndex[i];
//		double[] v = values[i];
//		for (int j = 0; j < r.Length; j++)
//		{
//			int ri = r[j];
//			sum += v[j] * xIn[ri];
//		}
//		xOut[i] = sum;
//	}
//}
//
//public void PreMultiplyOffset(double[] xIn, double[] xOut, int startOut, int offsetOut)
//{
//	for (int i = startOut; i < (n + offsetOut); i += offsetOut)
//	{
//		xOut[i] = 0.0;
//	}
//	int i = 0;
//	for (int k = startOut; i < n; k += offsetOut)
//	{
//		double sum = 0.0;
//		int[] r = rowIndex[i];
//		double[] v = values[i];
//		for (int j = 0; j < r.Length; j++)
//		{
//			int ri = r[j];
//			sum += v[j] * xIn[ri];
//		}
//		xOut[k] = sum;
//		i++;
//	}
//}

//public MultigridContractionSolver.ColMatrix Transpose()
//{
//	MultigridContractionSolver.ColMatrixCreator C = new MultigridContractionSolver.ColMatrixCreator(n, m);
//	for (int i = 0; i < n; i++)
//	{
//		int[] r = rowIndex[i];
//		double[] v = values[i];
//		for (int j = 0; j < r.Length; j++)
//		{
//			C.rowIndex[r[j]].Add(i);
//			C.values[r[j]].Add(v[j]);
//		}
//	}
//	return C;
//}

//public int Levels
//{
//	get
//	{
//		return resolutions.Count;
//	}
//}

//public List<Resolution> Resolutions
//{
//	get
//	{
//		return resolutions;
//	}
//}


