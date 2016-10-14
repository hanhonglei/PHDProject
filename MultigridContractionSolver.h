#pragma once
#include "jmspmesh/jmsmesh.h"
#include "CSMatrix.h"
#include <map>
#include <vector>
#include <set>
using namespace std;

class Resolution
{
	//Set<int>[] adjFace;
	//Set<int>[] adjVertex;
public:
	vector<set<int>> adjFace;
	vector<set<int>> adjVertex;
	vector<int> collapsedIndex;
	vector<int> faceList;
	vector<int> vertexList;
	//Dictionary<int, double>[] weight;
	vector<map<int, double>> weight;
public:
	Resolution(){ }
	~Resolution()
	{
		//delete []adjFace;
		//delete []adjVertex;
		//delete []collapsedIndex;
		//delete []faceList;
		//delete []vertexList;
	}
	Resolution(jmsMesh *mesh)
	{
		int vn = mesh->getNumVerts();
		vertexList.reserve(vn);
		//faceList = (int[]) mesh.FaceIndex.Clone();
		faceList.reserve(mesh->getNumTriangles()*3);
		for (int i = 0; i < mesh->getNumTriangles(); i++)
		{
			faceList.push_back(mesh->getTri(i).getVert1Index());
			faceList.push_back(mesh->getTri(i).getVert2Index());
			faceList.push_back(mesh->getTri(i).getVert3Index());

			//faceList[i*3] = mesh->getTri(i).getVert1Index();
			//faceList[i*3+1] = mesh->getTri(i).getVert2Index();
			//faceList[i*3+2] = mesh->getTri(i).getVert3Index();
		}
		//adjVertex = new Set<int>[vn];
		//adjFace = new Set<int>[vn];
		adjVertex.reserve(vn);
		adjFace.reserve(vn);		
		for (int i = 0; i < vn; i++)
		{
			vertexList.push_back(i);
			//vertexList[i] = i;
			//adjVertex[i] = new Set<int>(8);
			//adjFace[i] = new Set<int>(8);
			adjVertex.push_back(mesh->getVertex(i).getVertNeighbors());
			adjFace.push_back(mesh->getVertex(i).getTriNeighbors());
		}
	}
};


class ColMatrix
{
public:
	int m;
	int n;
	//public int[][] rowIndex;
	//public double[][] values;
	vector<vector<int>> rowIndex;
	vector<vector<double>> values;


	ColMatrix(int _m, int _n);

	int ATACG(double* x, double* b, double eps, int maxIter);

	void Multiply(double* xIn, double* xOut);

	int NumOfElements();

	void PreMultiply(double* xIn, double* xOut);

	int PColumnSize();
	int PRowSize();
};


class ColMatrixCreator
{
public:
	int m;
	int n;
	//List<int>[] rowIndex;
	//List<double>[] values;
	vector<vector<int>> rowIndex;
	vector<vector<double>> values;

	ColMatrixCreator(int _m, int _n);
	void AddValueTo(int i, int j, double value);
	ColMatrix ToColMatrix();

	int ColumnSize();

	int RowSize();
};

class MultigridContractionSolver
{
public:
	MultigridContractionSolver(void);
	~MultigridContractionSolver(void);

	double** SolveSystem(double* lapWeight, double* posWeight);

	MultigridContractionSolver(	jmsMesh* _newmesh);
	void* Factorization(CCSMatrix &C);

	ColMatrix BuildMatrixA(Resolution &r, double* lapWeight, double* posWeight);

	void GetNextLevel(Resolution &r, Resolution &result);

	CCSMatrix* MultiplyATA(ColMatrix A);
	// 利用动态链接库使用函数 [6/26/2011 Han Honglei]
	HINSTANCE hDllInst;
private:
	double convergeRatio;
	jmsMesh *mesh;
	//List<Resolution> resolutions;
	vector<Resolution> resolutions;
	int targetVertexCount;

};
