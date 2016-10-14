// ******************************
// mesh.cpp
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

#include <assert.h>
#include <float.h>

#if defined (_MSC_VER) && (_MSC_VER >= 1020)
#pragma warning(disable:4710) // function not inlined
#pragma warning(disable:4702) // unreachable code
#pragma warning(disable:4514) // unreferenced inline function has been removed
#endif

#include <iostream>

#include "jmsmesh.h"


jmsMesh::jmsMesh(char* filename)
{
	_numVerts = _numTriangles = 0;
	if (!loadFromFile(filename))
	{
		// we failed to load mesh from the file
		_numVerts = _numTriangles = 0;
		_vlist.clear();
		_plist.clear();
	}
	_numActiveTris = _numTriangles;
	_numActiveVerts = _numVerts;

	_vNewIndexList = NULL;
	_bHasSegment = false;
}

jmsMesh::jmsMesh(const jmsMesh& m)
{
	_numVerts = m._numVerts;
	_numTriangles = m._numTriangles;

	_numActiveTris = m._numActiveTris;
	_numActiveVerts = m._numActiveVerts;

	_vlist = m._vlist; // NOTE: triangles are still pointing to original mesh
	_plist = m._plist;

	_bHasSegment = m._bHasSegment;

	// NOTE: should reset tris in _vlist, _plist
	_vNewIndexList = NULL;

}

jmsMesh& jmsMesh::operator=(const jmsMesh& m)
{
	if (this == &m) return *this; // don't assign to self
	_numVerts = m._numVerts;
	_numTriangles = m._numTriangles;
	
	_numActiveTris = m._numActiveTris;
	_numActiveVerts = m._numActiveVerts;

	_vlist = m._vlist; // NOTE: triangles are still pointing to original mesh
	_plist = m._plist;
	// 将每个三角形的网格指针也进行更新，这样，彻底实现了重新赋值 [7/6/2011 Han Honglei]
	for (unsigned i = 0; i < _plist.size(); ++i)
		_plist[i].updateMeshPointer(this);

	_bHasSegment = m._bHasSegment;
	// NOTE: should reset tris in _vlist, _plist
	return *this;
}

jmsMesh::~jmsMesh()
{
	_numVerts = _numTriangles = 0;
	_vlist.erase(_vlist.begin(), _vlist.end());
	_plist.erase(_plist.begin(), _plist.end());

	_vSimplist.erase(_vSimplist.begin(), _vSimplist.end());
	if (_vNewIndexList != NULL)
	{
		delete []_vNewIndexList;
		_vNewIndexList = NULL;
	}
}

// Helper function for reading PLY mesh file
bool jmsMesh::readNumPlyVerts(FILE *&inFile, int& nVerts)
{
	// Read # of verts
	bool bElementFound = false;
	/* Get number of vertices in mesh*/
	for(;;)
	{
		char tempStr[1024];
		fscanf(inFile, "%s", tempStr);
		if (feof(inFile))
		{
			MessageBox(NULL, LPCSTR("Reached End of File and string \"element vertex\" not found!\n"),
				NULL, MB_ICONEXCLAMATION);
			return false;
		}

		/* change tempStr to lower case */
		ChangeStrToLower(tempStr);

		if (bElementFound && !strncmp(tempStr, "vertex", 6))
		{
			break;
		}

		if (!strncmp(tempStr, "element", 7))
		{
			bElementFound = true;
			continue;
		}
	}

	fscanf(inFile, "%d", &nVerts); 
	if (feof(inFile))
	{
		MessageBox(NULL, LPCSTR("Reached End of File before \"element face\" found!\n"),
			NULL, MB_ICONEXCLAMATION);
		return false;
	}
	return true;
}

// Helper function for reading PLY mesh file
bool jmsMesh::readNumPlyTris(FILE *&inFile, int& nTris)
{
	bool bElementFound = false;
	/* Get number of faces in mesh*/
	for(;;)
	{
		char tempStr[1024];
		fscanf(inFile, "%s", tempStr);
		if (feof(inFile))
		{
			MessageBox(NULL, LPCSTR("Reached End of File and string \"element face\" not found!\n"),
				NULL, MB_ICONEXCLAMATION);
			return false;
		}

		/* change tempStr to lower case */
		ChangeStrToLower(tempStr);

		if (bElementFound && !strncmp(tempStr, "face", 4))
		{
			break;
		}

		if (!strncmp(tempStr, "element", 7))
		{
			bElementFound = true;
			continue;
		}
	}

	fscanf(inFile, "%d", &nTris);
	if (feof(inFile))
	{
		MessageBox(NULL, LPCSTR("Reached End of File before list of vertices found!\n"),
			NULL, MB_ICONEXCLAMATION);
		return false;
	}
	return true;
}

// Helper function for reading PLY mesh file
bool jmsMesh::readPlyHeader(FILE *&inFile)
{
	char tempStr[1024];

	// Read "ply" string
	do
	{
		fscanf(inFile, "%s", tempStr);
		if (feof(inFile))
		{
			MessageBox(NULL, LPCSTR("Reached End of File and the string \"ply\" NOT FOUND!!\n"),
				NULL, MB_ICONEXCLAMATION);
			return false;
		}
		ChangeStrToLower(tempStr); // change tempStr to lower case 
	} while (strncmp(tempStr, "ply", 3));

	// Read # of verts
	if (!readNumPlyVerts(inFile, _numVerts))
	{
		return false;
	}

	// Read # of triangles
	if (!readNumPlyTris(inFile, _numTriangles))
	{
		return false;
	}

	// get end_header
	do
	{
		fscanf(inFile, "%s", tempStr);
		if (feof(inFile))
		{
			MessageBox(NULL, LPCSTR("Reached End of File and string \"end_header\" not found!\n"),
				NULL, MB_ICONEXCLAMATION);
			return false;
		}

		/* change tempStr to lower case */
		ChangeStrToLower(tempStr);
	} while (strncmp(tempStr, "end_header", 10));

	////////// end of header
	return true;
}

// Helper function for reading PLY mesh file
bool jmsMesh::readPlyVerts(FILE *&inFile)
{
	int i;
	// read vertices
	for (i = 0; i < _numVerts; i++)
	{
		char tempStr[1024];

#pragma warning(disable:4244)		/* disable double -> float warning */
		fscanf(inFile, "%s", tempStr);
		float x = atof(tempStr); 
		fscanf(inFile, "%s", tempStr);
		float y = atof(tempStr); 
		fscanf(inFile, "%s", tempStr);
		float z = atof(tempStr); 
#pragma warning(default:4244)		/* double -> float */

		vertex v(x, y, z);
		v.setIndex(i);

		_vlist.push_back(v); // push_back puts a *copy* of the element at the end of the list
		if (feof(inFile))
		{
			MessageBox(NULL, LPCSTR("Reached End of File before all vertices found!\n"),
				NULL, MB_ICONEXCLAMATION);
			return false;
		}

		// read until end of line
		while (fgetc(inFile) != '\n');
	}
	return true;
}

// Helper function for reading PLY mesh file
bool jmsMesh::readPlyTris(FILE *&inFile)
{
	int i;
	// read triangles
	for (i = 0; i < _numTriangles; i++)
	{
		int v1, v2, v3;
		int nVerts;
		fscanf(inFile, "%d", &nVerts);
		if (3 != nVerts)
		{
			MessageBox(NULL, LPCSTR("Error:  Ply file contains polygons which are not triangles!\n"),
				NULL, MB_ICONEXCLAMATION);
			return false;
		}
		fscanf(inFile, "%d", &v1);   // get value for vertex A
		fscanf(inFile, "%d", &v2);   // get value for vertex B
		fscanf(inFile, "%d", &v3);   // get value for vertex C

		// make sure verts in correct range
		assert(v1 < _numVerts && v2 < _numVerts && v3 < _numVerts);

		triangle t(this, v1, v2, v3);
		t.setIndex(i);

		_plist.push_back(t); // push_back puts a *copy* of the element at the end of the list

		// update each vertex w/ its neighbors (vertrices & triangles)
		_vlist[v1].addTriNeighbor(i);
		_vlist[v1].addVertNeighbor(v2);
		_vlist[v1].addVertNeighbor(v3);

		_vlist[v2].addTriNeighbor(i);
		_vlist[v2].addVertNeighbor(v1);
		_vlist[v2].addVertNeighbor(v3);

		_vlist[v3].addTriNeighbor(i);
		_vlist[v3].addVertNeighbor(v1);
		_vlist[v3].addVertNeighbor(v2);

		if (feof(inFile))
		{
			MessageBox(NULL, LPCSTR("Reached End of File before all faces found!\n"),
				NULL, MB_ICONEXCLAMATION);
			return false;
		}
		// read until end of line
		while (fgetc(inFile) != '\n');
	}
	return true;
}


// Load mesh from PLY file
bool jmsMesh::loadFromFile(char* filename)
{
	FILE* inFile = fopen(filename, "rt");
	if (inFile == NULL)
	{
		char pszError[_MAX_FNAME + 1];
		sprintf(pszError, "%s does not exist!\n", filename);
		MessageBox(NULL, LPCSTR(pszError), NULL, MB_ICONEXCLAMATION);
		return FALSE;
	}

	// read header to PLY file
	if (!readPlyHeader(inFile))
	{
		return false;
	}

	// read vertex data from PLY file
	if (!readPlyVerts(inFile))
	{
		return false;
	}

	// read triangle data from PLY file
	if (!readPlyTris(inFile))
	{
		return false;
	}

	fclose(inFile); // close the file

	calcVertNormals();

	return true;
}


// Recalculate the normal for one vertex
void jmsMesh::calcOneVertNormal(unsigned vert)
{
	vertex& v = getVertex(vert);
	/*const */set<int>& triset = v.getTriNeighbors();

	set<int>::iterator iter;

	vec3 vec;

	for (iter = triset.begin(); iter != triset.end(); ++iter)
	{
		// get the triangles for each vertex & add up the normals.
		vec += getTri(*iter).getNormalVec3();
	}

	vec.normalize(); // normalize the vertex	
	v.setVertNomal(vec);
}


// Calculate the vertex normals after loading the mesh.
void jmsMesh::calcVertNormals()
{
	// Iterate through the vertices
	for (unsigned i = 0; i < _vlist.size(); ++i)
	{
		calcOneVertNormal(i);
	}
}


// Used for debugging
void jmsMesh::dump()
{
	std::cout << "*** jmsMesh Dump ***" << std::endl;
	std::cout << "# of vertices: " << _numVerts << std::endl;
	std::cout << "# of triangles: " << _numTriangles << std::endl;
	for (unsigned i = 0; i < _vlist.size(); ++i)
	{
		std::cout << "\tVertex " << i << ": " << _vlist[i] << std::endl;
	}
	std::cout << std::endl;
	for (unsigned i = 0; i < _plist.size(); ++i)
	{
		std::cout << "\tTriangle " << i << ": " << _plist[i] << std::endl;
	}
	std::cout << "*** End of jmsMesh Dump ***" << std::endl;
	std::cout << std::endl;
}

// Get min, max values of all verts
void jmsMesh::setMinMax(float min[3], float max[3])
{
	max[0] = max[1] = max[2] = -FLT_MAX;
	min[0] = min[1] = min[2] = FLT_MAX;

	for (unsigned int i = 0; i < _vlist.size(); ++i)
	{
		const float* pVert = _vlist[i].getArrayVerts();
		if (pVert[0] < min[0]) min[0] = pVert[0];
		if (pVert[1] < min[1]) min[1] = pVert[1];
		if (pVert[2] < min[2]) min[2] = pVert[2];
		if (pVert[0] > max[0]) max[0] = pVert[0];
		if (pVert[1] > max[1]) max[1] = pVert[1];
		if (pVert[2] > max[2]) max[2] = pVert[2];
	}
}


void jmsMesh::meshTransform( const Matrix & mx )
{
	for (unsigned int i = 0; i < _vlist.size(); ++i)
	{
		_vlist[i].getXYZ() = transformLoc(mx, _vlist[i].getXYZ());
	}
	//if (vertextype == 0)
	//{	
	//	for(unsigned int i=0;i<nverts;i++)
	//	{
	//		verts[i] = transformLoc(mx,verts[i]);
	//	}
	//}

	//if (vertextype == 3)
	//{	
	//	for(unsigned int i=0;i<nverts;i++)
	//	{
	//		vertUVNs[i].vertex = transformLoc(mx,vertUVNs[i].vertex);
	//		vertUVNs[i].normal = transformVec(mx,vertUVNs[i].normal).getNormalized();
	//	}
	//}
}
// Center mesh around origin.
// Fit mesh in box from (-1, -1, -1) to (1, 1, 1)
void jmsMesh::Normalize()  
{
	float min[3], max[3], Scale;

	setMinMax(min, max);

	vec3 minv(min);
	vec3 maxv(max);

	vec3 dimv = maxv - minv;
	
	//if (dimv.x >= dimv.y && dimv.x >= dimv.z) Scale = 2.0f/dimv.x;
	//else if (dimv.y >= dimv.x && dimv.y >= dimv.z) Scale = 2.0f/dimv.y;
	//else Scale = 2.0f/dimv.z;

	//vec3 transv = minv + maxv;

	//transv *= 0.5f;

	//for (unsigned int i = 0; i < _vlist.size(); ++i)
	//{
	//	_vlist[i].getXYZ() -= transv;
	//	_vlist[i].getXYZ() *= Scale;
	//}

	float f_scale = 1.0 / dimv.length() * 1.5;
	vec3 trans = -f_scale * (minv + 0.5 * (maxv - minv));

	meshTransform(scale(f_scale ,f_scale , f_scale));
	meshTransform(translate(trans.x , trans.y , trans.z));
}

// 保存ply文件 [8/16/2010 Leo Han]
//ply
//format ascii 1.0           { ascii/binary, format version number }
//comment made by anonymous  { comments keyword specified, like all lines }
//comment this file is a cube
//element vertex 8           { define "vertex" element, 8 of them in file }
//property float32 x         { vertex contains float "x" coordinate }
//property float32 y         { y coordinate is also a vertex property }
//property float32 z         { z coordinate, too }
//element face 6             { there are 6 "face" elements in the file }
//property list uint8 int32 vertex_index { "vertex_indices" is a list of ints }
//end_header                 { delimits the end of the header }
bool jmsMesh::savePlyHeader(FILE *&inFile)
{
	char tempStr[1024];

	// header
	fprintf(inFile, "ply\n");
	fprintf(inFile, "format ascii 1.0\n");
	fprintf(inFile, "element vertex %d\n", _numVerts);
	fprintf(inFile, "property float32 x\nproperty float32 y\nproperty float32 z\n");
	fprintf(inFile, "element face %d\n", _numTriangles);
	fprintf(inFile, "property list uint8 int32 vertex_index\n");
	fprintf(inFile, "end_header\n");

	////////// end of header
	return true;
}

// Helper function for reading PLY mesh file
bool jmsMesh::savePlyVerts(FILE *&inFile)
{
	int i;
	vec3 v;

	// read vertices
	for (i = 0; i < _numVerts; i++)
	{
		v = _vlist[i].getXYZ();
		fprintf(inFile, "%f %f %f\n", v.x, v.y, v.z);
	}
	return true;
}

// Helper function for reading PLY mesh file
bool jmsMesh::savePlyTris(FILE *&inFile)
{
	int i;
	// read triangles
	for (i = 0; i < _numTriangles; i++)
		fprintf(inFile, "3 %d %d %d\n", _plist[i].getVert1Index(), _plist[i].getVert2Index(),_plist[i].getVert3Index());
	return true;
}

// save mesh to PLY file
bool jmsMesh::SaveFile(char* filename)
{
	FILE* inFile;
	errno_t err;

	// Open for read (will fail if file "crt_fopen_s.c" does not exist)
	err  = fopen_s( &inFile, filename, "wt");

	if (err != 0)
	{
		char pszError[_MAX_FNAME + 1];
		sprintf(pszError, "%s does not exist!\n", filename);
		MessageBox(NULL, LPCSTR(pszError), NULL, MB_ICONEXCLAMATION);
		return FALSE;
	}

	// read header to PLY file
	if (!savePlyHeader(inFile))
	{
		return false;
	}

	// read vertex data from PLY file
	if (!savePlyVerts(inFile))
	{
		return false;
	}

	// read triangle data from PLY file
	if (!savePlyTris(inFile))
	{
		return false;
	}

	fprintf(inFile, "\n");

	fclose(inFile); // close the file

	return true;
}


// 保存ply文件 [8/16/2010 Leo Han]
//ply
//format ascii 1.0           { ascii/binary, format version number }
//comment made by anonymous  { comments keyword specified, like all lines }
//comment this file is a cube
//element vertex 8           { define "vertex" element, 8 of them in file }
//property float32 x         { vertex contains float "x" coordinate }
//property float32 y         { y coordinate is also a vertex property }
//property float32 z         { z coordinate, too }
//element face 6             { there are 6 "face" elements in the file }
//property list uint8 int32 vertex_index { "vertex_indices" is a list of ints }
//end_header                 { delimits the end of the header }
bool jmsMesh::savePlyHeader(CArchive& ar)
{
	char tempStr[1024];

	CString header;
	header.Format(
		"ply\r\nformat ascii 1.0\r\nelement vertex %d\r\nproperty float32 x\r\nproperty float32 y\r\nproperty float32 z\r\nelement face %d\r\nproperty list uint8 int32 vertex_index\r\nend_header\r\n"
		,_numActiveVerts/*_numVerts*/,_numActiveTris);

	ar.WriteString(header);
	return true;
}

// Helper function for reading PLY mesh file
// Problem: about 1 to 3 vertexes are missing in the saved file [11/18/2011 Han]
bool jmsMesh::savePlyVerts(CArchive& ar)
{
	int i;
	vec3 v;

	CString verStr;
	// read vertices
	for (i = 0; i < _vSimplist.size(); i++)
	{
		v = _vSimplist[i].getXYZ();
		verStr.Format("%f %f %f\r\n", v.x,v.y , v.z);
		ar.WriteString(verStr);
	}
	return true;
}

// Helper function for reading PLY mesh file
bool jmsMesh::savePlyTris(CArchive& ar)
{
	int i;
	CString triStr;
	// read triangles
	for (i = 0; i < _numTriangles; i++)
	{
		if (!_plist[i].isActive())
			continue;

		triStr.Format("3 %d %d %d\r\n", _vNewIndexList[_plist[i].getVert1Index()]
		,_vNewIndexList[_plist[i].getVert2Index()],_vNewIndexList[_plist[i].getVert3Index()]);
		ar.WriteString(triStr);
	}
	return true;
}

// save mesh to PLY file
bool jmsMesh::SaveFile(CArchive& ar)
{
	RefreshData();
	// read header to PLY file
	if (!savePlyHeader(ar))
	{
		return false;
	}

	// read vertex data from PLY file
	if (!savePlyVerts(ar))
	{
		return false;
	}

	// read triangle data from PLY file
	if (!savePlyTris(ar))
	{
		return false;
	}

	return true;
}

bool jmsMesh::RefreshData()
{
	int nNewIndex = 0;
	vertex vSimp ;

	_vSimplist.clear();
	if (_vNewIndexList != NULL)
		delete[] _vNewIndexList;
	_vNewIndexList = new int[_numVerts];
	memset(_vNewIndexList, -1, _numVerts*sizeof(int));
	int nCurVer = -1;

	int nActiveVerTest = 0;

	for (int i=0;i <  _plist.size(); i++)
	{
		if (!_plist[i].isActive())
			continue;
		
		while (1)
		{
			nCurVer = -1;
			if(_vNewIndexList[_plist[i].getVert1Index()] == -1)
			{
				nCurVer = _plist[i].getVert1Index();
			}
			else if (_vNewIndexList[_plist[i].getVert2Index()] == -1)
			{
				nCurVer = _plist[i].getVert2Index();
			}
			else if (_vNewIndexList[_plist[i].getVert3Index()] == -1)
			{
				nCurVer = _plist[i].getVert3Index();
			}
			if (nCurVer != -1)
			{
				vSimp = _vlist[nCurVer];
				nNewIndex = _vSimplist.size();
				vSimp.setIndex(nNewIndex);
				_vNewIndexList[nCurVer] = nNewIndex;

				_vSimplist.push_back(vSimp);
				nActiveVerTest++;
			}
			else
				break;
		}
	}

	return true;
}

// 以下是用于求解骨架所使用的函数 [6/13/2011 Han Honglei]
void jmsMesh::calcTriNormals()
{
	// Iterate through the vertices
	for (unsigned i = 0; i < getNumTriangles(); ++i)
	{
		getTri(i).calcNormal();
	}
}

// 得到模型的体积 [6/13/2011 Han Honglei]
double jmsMesh::Volume()
{
	double totVolume = 0.0;
	int i = 0;
	for (int j = 0; i < /*faceCount*/getNumTriangles(); j += 3)
	{
		//int c1 = faceIndex[j] * 3;
		//int c2 = faceIndex[j + 1] * 3;
		//int c3 = faceIndex[j + 2] * 3;
		//Vector3d a = new Vector3d(vertexPos, c1);
		//Vector3d b = new Vector3d(vertexPos, c2);
		//Vector3d c = new Vector3d(vertexPos, c3);
		vec3 a = vec3(getTri(i).getVert1());
		vec3 b = vec3(getTri(i).getVert2());
		vec3 c = vec3(getTri(i).getVert3());

		totVolume += ((((((a.x * b.y) * c.z) - ((a.x * b.z) * c.y)) - ((a.y * b.x) * c.z)) + ((a.y * b.z) * c.x)) + ((a.z * b.x) * c.y)) - ((a.z * b.y) * c.x);
		i++;
	}
	return totVolume;
}


double jmsMesh::ComputeFaceArea(int fIndex)
{
	vec3 v1 = vec3(getTri(fIndex).getVert1());
	vec3 v2 = vec3(getTri(fIndex).getVert2());
	vec3 v3 = vec3(getTri(fIndex).getVert3());

	vec3 v = v2 - v1;
	return(v.cross(v3 - v1).length() / 2.0);
	//int b = fIndex * 3;
	//Vector3d v1 = new Vector3d(VertexPos, faceIndex[b] * 3);
	//Vector3d v2 = new Vector3d(VertexPos, faceIndex[b + 1] * 3);
	//Vector3d v3 = new Vector3d(VertexPos, faceIndex[b + 2] * 3);
	//Vector3d CS$0$0000 = v2 - v1;
	//return (CS$0$0000.Cross(v3 - v1).Length() / 2.0);
}
double jmsMesh::AverageFaceArea()
{
	double tot = 0.0;
	int i = 0;
	//for (int j = 0; i < faceCount; j += 3)
	//{
	//	Vector3d v1 = new Vector3d(VertexPos, faceIndex[j] * 3);
	//	Vector3d v2 = new Vector3d(VertexPos, faceIndex[j + 1] * 3);
	//	Vector3d v3 = new Vector3d(VertexPos, faceIndex[j + 2] * 3);
	//	Vector3d CS0000 = v2 - v1;
	//	tot += CS0000.Cross(v3 - v1).Length() / 2.0;
	//	i++;
	//}
	for (int j = 0; i < /*faceCount*/getNumTriangles(); j += 3)
	{
		tot += ComputeFaceArea(i);
		i++;
	}

	return (tot / ((double)getNumTriangles()));

}

bool jmsMesh::setVerPos(int index, vec3 pos)
{
	 if (index >= getNumVerts())
	 {
		 return false;
	 }
	 _vlist[index].setPos(pos);
	 return true;
}

bool jmsMesh::setVerPos(int index,double x, double y, double z)
{
	if (index >= getNumVerts())
	{
		return false;
	}
	_vlist[index].setPos(vec3(x, y, z));
	return true;	
}
