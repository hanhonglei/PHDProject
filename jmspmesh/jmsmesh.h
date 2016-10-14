// ******************************
// mesh.h
//
// jmsMesh class, which stores a list
// of vertices & a list of triangles.
//
// Jeff Somers
// Copyright (c) 2002
//
// jsomers@alumni.williams.edu
// March 27, 2002
// ******************************

#ifndef __mesh_h
#define __mesh_h

#if defined (_MSC_VER) && (_MSC_VER >= 1020)

#pragma once
#pragma warning(disable:4710) // function not inlined
#pragma warning(disable:4702) // unreachable code
#pragma warning(disable:4514) // unreferenced inline function has been removed
#define WIN32_LEAN_AND_MEAN
//#include <windows.h>
#include <afx.h>

#endif



#include <vector>
#include "vertex.h"
#include "triangle.h"
#include "Matrix.h"
using namespace std;


// jmsMesh class.  This stores a list of vertices &
// another list of triangles (which references the vertex list)
class jmsMesh
{
public:
	// Constructors and Destructors
	jmsMesh() {_numVerts = _numTriangles = 0; _numActiveTris = _numTriangles; _numActiveVerts = _numVerts; _vNewIndexList = NULL;};
	jmsMesh(char* filename); // passed name of mesh file
	~jmsMesh();

	jmsMesh(const jmsMesh&); // copy ctor
	jmsMesh& operator=(const jmsMesh&); // assignment op

	// Get list of vertices, triangles
	vertex& getVertex(int index) {return _vlist[index];};
	const vertex& getVertex(int index) const {return _vlist[index];};
	triangle& getTri(int index) {return _plist[index];};
	const triangle& getTri(int index) const {return _plist[index];};

	int getNumVerts() {return _numVerts;};
	void setNumVerts(int n) {_numVerts = n;};
	int getNumTriangles() {return _numTriangles;};
	void setNumTriangles(int n) {_numTriangles = n;};

	void Normalize();// center mesh around the origin & shrink to fit in [-1, 1]
	void meshTransform( const Matrix & mx );

	void calcOneVertNormal(unsigned vert); // recalc normal for one vertex

	void dump(); // print mesh state to cout
	bool SaveFile(char* filename); // load from PLY file
	bool SaveFile(CArchive& ar); // load from PLY file

	int _numActiveVerts;		// 经过编辑以后仍然存活的顶点和面片数目
	int _numActiveTris;

	vector<vertex> _vSimplist; // list of vertices in mesh
	int* _vNewIndexList;

	bool RefreshData();

	// 用于骨架抽取 [6/21/2011 Han Honglei]
	double Volume();
	double ComputeFaceArea(int fIndex);
	void calcTriNormals();
	double AverageFaceArea();
	void calcVertNormals(); // Calculate the vertex normals after loading the mesh

	bool setVerPos(int index, vec3 pos);
	bool setVerPos(int index,double x, double y, double z);

	bool _bHasSegment;

	

private:
	vector<vertex> _vlist; // list of vertices in mesh
	vector<triangle> _plist; // list of triangles in mesh

	int _numVerts;
	int _numTriangles;


	bool operator==(const jmsMesh&); // don't allow op== -- too expensive
	
	bool loadFromFile(char* filename); // load from PLY file

	void ChangeStrToLower(char* pszUpper)
	{
		for(char* pc = pszUpper; pc < pszUpper + strlen(pszUpper); pc++) {
			*pc = (char)tolower(*pc);
		}
	}

	// get bounding box for mesh
	void setMinMax(float min[3], float max[3]);

	void calcAllQMatrices(jmsMesh& mesh); // used for Quadrics method


	// Helper function for reading PLY mesh file
	bool readNumPlyVerts(FILE *&inFile, int& nVerts);
	bool readNumPlyTris(FILE *&inFile, int& nTris);
	bool readPlyHeader(FILE *&inFile);
	bool readPlyVerts(FILE *&inFile);
	bool readPlyTris(FILE *&inFile);

	bool savePlyHeader(FILE *&inFile);
	bool savePlyVerts(FILE *&inFile);
	bool savePlyTris(FILE *&inFile);

	bool savePlyHeader(CArchive& ar);
	bool savePlyVerts(CArchive& are);
	bool savePlyTris(CArchive& ar);

};

#endif // __mesh_h
