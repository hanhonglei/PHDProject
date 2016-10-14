// ******************************
// vec3.h
//
// Vector class
// 
// Jeff Somers
// Copyright (c) 2002
//
// jsomers@alumni.williams.edu
// March 27, 2002
// ******************************

#ifndef __Vec3_h
#define __Vec3_h

#if defined (_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#pragma warning(disable:4710) // function not inlined
#pragma warning(disable:4702) // unreachable code
#pragma warning(disable:4514) // unreferenced inline function has been removed
#endif

#include <math.h>
#include <iostream>

// vec3 is a vector class
// I couldn't call it vec3, because that already
// seems to be defined.
class vec3 
{
public:
	float x, y, z;

	enum {MAX_INPUT_STRING = 40};

	// Constructors and Destructors
	vec3() {x = y = z = 0.0f;};
	vec3(float x1, float y1, float z1) {x = x1; y = y1; z = z1;};
	vec3(float av[3]) {x = av[0]; y = av[1]; z = av[2];};
	vec3(const vec3& v) {x = v.x; y = v.y; z = v.z;};
	//  [6/21/2011 Han Honglei]
	vec3(const float* av) {x = av[0]; y = av[1]; z = av[2];};

	~vec3() {}; // Destructor intentially does nothing

	// Assignment operator
	vec3& operator=(const vec3& v) {x = v.x; y = v.y; z = v.z; return *this;};

	// Comparision operators
	bool operator==(const vec3& v) {return (x == v.x && y == v.y && z == v.z);};
	bool operator!=(const vec3& v) {return (x != v.x || y != v.y || z != v.z);};

	// Scalar operations
	vec3 operator+(float f) const {return vec3(x + f, y + f, z + f);};
	vec3 operator-(float f) const {return vec3(x - f, y - f, z - f);};
	vec3 operator*(float f) const {return vec3(x * f, y * f, z * f);};
	vec3 operator/(float f) const {vec3 v1(x,y,z); if (f != 0.0f) {v1.x /= f; v1.y /= f; v1.z /= f;}; return v1;};

	vec3& operator+=(float f) {x += f; y += f; z += f; return *this;};
	vec3& operator-=(float f) {x -= f; y -= f; z -= f; return *this;};
	vec3& operator*=(float f) {x *= f; y *= f; z *= f; return *this;};
	vec3& operator/=(float f) {if(f!=0.0f){ x /= f; y /= f; z /= f;}; return *this;};


	// Vector operations
	vec3 operator+(const vec3& v) const {return vec3(x + v.x, y + v.y, z + v.z);};
	vec3& operator+=(const vec3& v) {x += v.x; y += v.y; z += v.z; return *this;};
	vec3 operator-(const vec3& v) const {return vec3(x - v.x, y - v.y, z - v.z);};
	vec3& operator-=(const vec3& v) {x -= v.x; y -= v.y; z -= v.z; return *this;};

	// Unary operators
	vec3 operator-() const {return vec3 (-x, -y, -z); };

	// Dot and Cross Products
	float dot(const vec3& v) const {return (x * v.x + y * v.y + z * v.z);};
	vec3 cross(const vec3& v) const {return vec3(y * v.z - z * v.y,
											 z * v.x - x * v.z,
											 x * v.y - y * v.x);};
	vec3 unitcross(const vec3& v) const {vec3 vr(y * v.z - z * v.y,
											 z * v.x - x * v.z,
											 x * v.y - y * v.x); vr.normalize(); return vr;};

	// Miscellaneous
	void normalize() {float a = float(sqrt(x*x + y*y + z*z)); if (a!=0.0f) {x/=a; y/=a; z/=a;};};
	void setZero() {x = y = z = 0.0f;};
	float length() {return float(sqrt(x*x + y*y + z*z));};

	// Friend functions
	friend vec3 operator*(float a, const vec3& v) {return vec3 (a * v.x, a * v.y, a * v.z);};

	// dot and cross products
	float dot(const vec3& v1, const vec3& v2) {return (v1.x * v2.x + v1.y * v2.y +v1. z * v2.z);};
	vec3 cross(const vec3& v1, const vec3& v2)  {return vec3 (v1.y * v2.z - v1.z * v2.y,
														v1.z * v2.x - v1.x * v2.z,
															v1.x * v2.y - v1.y * v2.x);};
	vec3 unitcross(const vec3& v1, const vec3& v2)  {vec3 vr(v1.y * v2.z - v1.z * v2.y,
													v1.z * v2.x - v1.x * v2.z,
													v1.x * v2.y - v1.y * v2.x); 
													vr.normalize(); return vr;};


	// Input and Output
	friend std::ostream& operator<<(std::ostream& os, const vec3& vo);

	// han [4/28/2012 Han]
	float DistToTriangle( const vec3 *triangle);
	float Dist2ToTriangle( const vec3 *triangle);
	float length2() {return (x*x + y*y + z*z);};


//	friend istream& operator>>(istream& is, vec3& vi);

private:
};

#endif // #ifndef __Vec3_h