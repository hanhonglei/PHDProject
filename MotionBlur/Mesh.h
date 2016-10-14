
/*****************************************************************************\
 *                                                                           *
 * filename : Mesh.h                                                         *
 * last mod : 08/06/02                                                       *
 * author   : R. Keith Morley / Peter Shirley                                *
 *                                                                           *
\*****************************************************************************/

#ifndef _MESH_H_
#define _MESH_H_ 1

#include "Vertex.h"
#include "Surface.h"
#include "Material.h"
//#include "Matrix.h"

class Matrix;

class Mesh 
{
public:  
	Mesh(){}
	Mesh( const Mesh & ms);
   ~Mesh();
   Material* getMaterial(int index)const;
   
   bool ReadObjFile(char* filename);

   void meshTransform( const Matrix & mx);
	void calTexCoord( const char & type );
	
	//calculate every vertex's shift in the image plane.
	void calVertexShift( const Matrix curmx ,const Matrix & mx); 

	//calculate all the area of the triangles
	void calTriangleArea();
	
   // data members
   Material**  mptr;
   
   Vector3*   verts;//points
   VertexN*   vertNs;
   VertexUV*  vertUVs;//vertex index
   VertexUVN* vertUVNs;//vertex


public:
	unsigned int ntris, nverts;
	int *vertexIndex;
	float smptime;
	float * shifts;
	float * triAreas;
	float * triDepths;
	Vector3 m_min,m_max;

	short vertextype;//0 for vector3,1 for vertexN, 2 for vertexUV, 3for vertexUVN
};

inline
Material* Mesh::getMaterial(int index)const { return mptr[index]; }

#endif // _MESH_H_
