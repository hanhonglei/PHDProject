// ******************************
// vec3.cpp
//
// Vector class
// 
// Jeff Somers
// Copyright (c) 2002
//
// jsomers@alumni.williams.edu
// March 27, 2002
// ******************************

// Friend functions for input/output

#if defined (_MSC_VER) && (_MSC_VER >= 1020)
#pragma warning(disable:4710) // function not inlined
#pragma warning(disable:4702) // unreachable code
#pragma warning(disable:4514) // unreferenced inline function has been removed
#endif

#include "vec3.h"

std::ostream&
operator<<(std::ostream& os, const vec3& vo)
{
	return os << "<" << vo.x << ", " << vo.y << ", " << vo.z << ">";
}

float clamp( float v, float minV, float maxV)
{
	v = v < minV ? minV : v;
	v = v > maxV ? maxV : v;
	return v;
}
// 添加支持点到三角形最近距离的方法 [4/28/2012 Han]
// source: http://www.arongranberg.com/astar/docs/index.php
vec3 closesPointOnTriangle( const vec3 *triangle, const vec3 &sourcePosition )
{
	vec3 edge0 = triangle[1] - triangle[0];
	vec3 edge1 = triangle[2] - triangle[0];
	vec3 v0 = triangle[0] - sourcePosition;

	float a = edge0.dot( edge0 );
	float b = edge0.dot( edge1 );
	float c = edge1.dot( edge1 );
	float d = edge0.dot( v0 );
	float e = edge1.dot( v0 );

	float det = a*c - b*b;
	float s = b*e - c*d;
	float t = b*d - a*e;

	if ( s + t < det )
	{
		if ( s < 0.f )
		{
			if ( t < 0.f )
			{
				if ( d < 0.f )
				{
					s = clamp( -d/a, 0.f, 1.f );
					t = 0.f;
				}
				else
				{
					s = 0.f;
					t = clamp( -e/c, 0.f, 1.f );
				}
			}
			else
			{
				s = 0.f;
				t = clamp( -e/c, 0.f, 1.f );
			}
		}
		else if ( t < 0.f )
		{
			s = clamp( -d/a, 0.f, 1.f );
			t = 0.f;
		}
		else
		{
			float invDet = 1.f / det;
			s *= invDet;
			t *= invDet;
		}
	}
	else
	{
		if ( s < 0.f )
		{
			float tmp0 = b+d;
			float tmp1 = c+e;
			if ( tmp1 > tmp0 )
			{
				float numer = tmp1 - tmp0;
				float denom = a-2*b+c;
				s = clamp( numer/denom, 0.f, 1.f );
				t = 1-s;
			}
			else
			{
				t = clamp( -e/c, 0.f, 1.f );
				s = 0.f;
			}
		}
		else if ( t < 0.f )
		{
			if ( a+d > b+e )
			{
				float numer = c+e-b-d;
				float denom = a-2*b+c;
				s = clamp( numer/denom, 0.f, 1.f );
				t = 1-s;
			}
			else
			{
				s = clamp( -e/c, 0.f, 1.f );
				t = 0.f;
			}
		}
		else
		{
			float numer = c+e-b-d;
			float denom = a-2*b+c;
			s = clamp( numer/denom, 0.f, 1.f );
			t = 1.f - s;
		}
	}

	return triangle[0] + s * edge0 + t * edge1;
}

float vec3::DistToTriangle( const vec3 *triangle)
{
	return (*this - closesPointOnTriangle(triangle, *this)).length();
}

float vec3::Dist2ToTriangle( const vec3 *triangle)
{
	return (*this - closesPointOnTriangle(triangle, *this)).length2();
}
/* NOT IMPLEMENTED
istream& operator>>(istream &io, vec3 &vi)
{
	char inBuf[vec3::MAX_INPUT_STRING];
	io >> inBuf; // operator>>(ostream&, char*); -- or is it istream?
	//!FIX need to convert string to vector here
//	vi = inBuf;	// String::operator=(const char*)
	return io;
}

*/

