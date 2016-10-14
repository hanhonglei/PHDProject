#include "pmesh.h"

const float PI = 3.141592654f;

float getpixelsize(float fovy, float planedist , float imageHeight)
{
	return tan(fovy/2.0f/180.0f*PI)*planedist/imageHeight*2.0f;
}

bool PMesh::CalcBlurWeights(/*vec3 eye,vec3 gaze,vec3 up*/)
{
	//jmsMesh m;
	Matrix curmx = translate(0,0,0);
	Matrix mx1 = translate(1,1,1);
	Matrix mx2 = rotate( vec3(1,1,1) , 60);
	double areasum = 0 , perareasum = 0;
	Matrix projectmx , transfermx , permx;
	float imgplanedist=0.5;
	vec3 eye = vec3(0,0,2);
	vec3 gaze = vec3(0,0,-1);
	vec3 up = vec3(0,1,0);

	//m.ReadObjFile("Wei_Bunny01(0.00081).obj");
	//m.calTriangleArea();
	//vec3 maxscale = m.m_max - m.m_min;
	//float f_scale = max(max( maxscale.x() , maxscale.y()) , maxscale.z());
	//m.meshTransform(scale(1.f/f_scale , 1.f/f_scale , 1.f/f_scale));
	//vec3 diag = m.m_max - m.m_min;
	//float f_scale = 1.0 / diag.length() * 1.5;
	//vec3 trans = -f_scale * (m.m_min + 0.5 * (m.m_max - m.m_min));
	//m.meshTransform(translate(trans.x , trans.y , trans.z));
	//m.meshTransform(scale(f_scale ,f_scale , f_scale));
	//m.calTriangleArea();
	//std::ofstream fin("d:\\Wei_Bunny01(0.00081).blr"); 

	//
	projectmx = viewMatrix( eye , gaze , up );
	transfermx = translate(0.128,0,0);
	permx = zeroMatrix();
	permx.x[0][0] = 1; permx.x[1][1] = 1; permx.x[2][2] = 1; permx.x[3][2] = 1/imgplanedist;

	vec3 orignalpos , laterpos;
	float pixelsize = getpixelsize(30 , 1 , 600);
	float dd,sumdd=0;
	for (int i = 0; i < _newmesh.getNumVerts(); i++)
	{
		orignalpos = transformLoc( permx*projectmx , _newmesh.getVertex(i).getXYZ());
		laterpos = transformLoc( transfermx , _newmesh.getVertex(i).getXYZ());
		dd = (orignalpos - transformLoc( permx*projectmx , laterpos )).length();
		_newmesh.getVertex(i)._blur = dd;
		_mesh->getVertex(i)._blur = dd;

		//dd = dist( orignalpos , transformLoc( permx*projectmx , laterpos ) );
		//printf("%f\n",dd);
		sumdd += dd;
	}
	return true;
}