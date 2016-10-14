// ******************************
// triangle.cpp
//
// Triangle class.
// Contains vetices & normal to
// the place defined by the 3
// vertices.
//
// Jeff Somers
// Copyright (c) 2002
//
// jsomers@alumni.williams.edu
// March 27, 2002
// ******************************

#if defined (_MSC_VER) && (_MSC_VER >= 1020)
#pragma warning(disable:4710) // function not inlined
#pragma warning(disable:4702) // unreachable code
#pragma warning(disable:4514) // unreferenced inline function has been removed
#endif

#include <assert.h>
#include "triangle.h"
#include "jmsmesh.h"

// retrieve vertices as an array of floats
const float* triangle::getVert1() {return (_mesh->getVertex(_vert1)).getArrayVerts();}
const float* triangle::getVert2() {return (_mesh->getVertex(_vert2)).getArrayVerts();}
const float* triangle::getVert3() {return (_mesh->getVertex(_vert3)).getArrayVerts();}

// retrieve vertices as a vertex object
vertex& triangle::getVert1vertex() const {return _mesh->getVertex(_vert1);};
vertex& triangle::getVert2vertex() const {return _mesh->getVertex(_vert2);};
vertex& triangle::getVert3vertex() const {return _mesh->getVertex(_vert3);};

void triangle::setVertexesLayer(short nLayer)
{
	_mesh->getVertex(_vert1).setLayer(nLayer); 
	_mesh->getVertex(_vert2).setLayer(nLayer); 
	_mesh->getVertex(_vert3).setLayer(nLayer);
}

// Calculate normal of triangle
void
triangle::calcNormal()
{
	assert(_mesh);

	vec3 vec1 = (_mesh->getVertex(_vert1)).getXYZ();
	vec3 vec2 = (_mesh->getVertex(_vert2)).getXYZ();
	vec3 vec_3 = (_mesh->getVertex(_vert3)).getXYZ();

	vec3 veca = vec2 - vec1;
	vec3 vecb = vec_3 - vec2;

	_normal = veca.unitcross(vecb);
	// Note that if the triangle is degenerate (all vertices lie in a line),
	// the normal will be <0,0,0>

	// This is the "d" from the plane equation ax + by + cz + d = 0;

	_d = -_normal.dot(vec1);
}

// Calculate area of triangle
float
triangle::calcArea()
{
	assert(_mesh);

	// If a triangle is defined by 3 points, say p, q and r, then
	// its area is 0.5 * length of ((p - r) cross (q - r))
	// See Real-Time Rendering book, Appendix A
	vec3 vec1 = (_mesh->getVertex(_vert1)).getXYZ();
	vec3 vec2 = (_mesh->getVertex(_vert2)).getXYZ();
	vec3 vec_3 = (_mesh->getVertex(_vert3)).getXYZ();
	vec3 vecA = vec1 - vec2;
	vec3 vecB = vec_3 - vec2;

	vec3 cross = vecA.cross(vecB);
	float area = float(0.5 * cross.length());
	return area;
}


// Used for output
std::ostream&
operator<<(std::ostream& os, const triangle& to)
{
	os << "vert1: " << to._vert1 << " vert2: " << to._vert2 << " vert3: " << to._vert3; // for some reason this isn't working as a friend function, not sure why
	os << " Normal: " << to._normal << " Active? " << to.isActive();
	os << " Index: " << to._index;

	// it is pulling ostream from the STL typedef, not the regular ostream, though.
	return os;
}

