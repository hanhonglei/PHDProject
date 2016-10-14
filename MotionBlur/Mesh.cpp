#define _CRT_SECURE_NO_DEPRECATE

#include <vector>
#include "Mesh.h"
#include "Matrix.h"

#define LINE_BUFFER_SIZE 1024
// GETNUM just gets the next number from a line of input in an OBJ file
#ifndef GETNUM
#define GETNUM(lineBuffer, numBuffer, lindex, nindex, tval)  \
	nindex=0;\
	while ((lineBuffer[lindex] == ' ') || lineBuffer[lindex] == '/') lindex++;\
	while ((lineBuffer[lindex] != ' ') && (lineBuffer[lindex] != '/') && \
	(lineBuffer[lindex] != '\0') && (lineBuffer[lindex] != '\n') && (lindex != LINE_BUFFER_SIZE)) { \
	numBuffer[nindex] = lineBuffer[lindex]; \
	nindex++; \
	lindex++; \
} \
	numBuffer[nindex] = '\0'; \
	tval = atoi(numBuffer);
#endif

using namespace std;


Mesh::~Mesh() {

	// data members
	//what's going on here?
	//delete[] mptr;
	//delete[] verts;
	//delete[] vertUVs;
	//delete[] vertNs;
	//delete[] vertUVNs;
	//delete[] vertexIndex;
}




bool Mesh::ReadObjFile(char *filename)
{
	FILE* fin;
	fin = fopen(filename, "r");
	if (!fin) {
		std::cout << "Failed to open file" << std::endl;
		return false;	
	}

	// temporary input buffers
	vector<Vector3> points;//保存点的坐标
	vector<Vector3> file_normals;
	vector<Vector2> file_uvvector;//保存纹理坐标

	vector<int> vertsIndex;//保存三角形对应顶点的索引
	vector<int> normalIndex;
	vector<int> uvIndex;


	Vector3 ptmp;
	Vector3 ntmp;
	Vector2 uvtmp;
	float uv1, uv2;

	char lineBuffer[LINE_BUFFER_SIZE];
	char numBuffer[32];
	int lindex=0;
	int nindex=0;
	int ival, uvval, nval;
	ntris=0;
	nverts=0;

	// parse the data in
	while (fgets(lineBuffer, LINE_BUFFER_SIZE, fin)) {
		switch (lineBuffer[0]) {
			case 'v':
				// case vertex information
				if (lineBuffer[1] == ' ') {
					// regular vertex point
					nverts++;
					float x,y,z;
					sscanf(&lineBuffer[2], "%f %f %f", &x, &y,&z);
					
					if (x>m_max.x())
						m_max.setX(x);
					if (x<m_min.x())
						m_min.setX(x);
					
					if (y>m_max.y())
							m_max.setY(y);
					if (y<m_min.y())
						m_min.setY(y);
					
					if (z>m_max.z())
						m_max.setZ(z);
					if (z<m_min.z())
						m_min.setZ(z);
					
					ptmp = Vector3(x,y,z);
					points.push_back(ptmp);
				} else if (lineBuffer[1] == 't') {
					// texture coordinates
					sscanf(&lineBuffer[3], "%f %f", &uv1, &uv2);
					uvtmp = Vector2(uv1,uv2);
					file_uvvector.push_back(uvtmp);
				} else if (lineBuffer[1] == 'n') {
					// normal vector
					sscanf(&lineBuffer[2], "%f %f %f", &ntmp.e[0], &ntmp.e[1], &ntmp.e[2]);
					file_normals.push_back(ntmp);
				}
				break;
			case 'f':
				// case face information
				lindex = 2;//每行上的索引
				ntris++;
				for (int i=0; i < 3; i++) {

					GETNUM(lineBuffer, numBuffer, lindex, nindex, ival)

						// obj files go from 1..n, this just allows me to access the memory
						// directly by droping the index value to 0...(n-1)
						ival--;
					vertsIndex.push_back(ival);//三角形对应第i个顶点的points索引

					if (lineBuffer[lindex] == '/') {
						lindex++;
						GETNUM(lineBuffer, numBuffer, lindex, nindex, uvval)
							uvIndex.push_back(uvval-1);//三角形对应第i个顶点的纹理索引
					}

					if (lineBuffer[lindex] == '/') {
						lindex++;
						GETNUM(lineBuffer, numBuffer, lindex, nindex, nval)
							normalIndex.push_back(nval-1);//三角形对应第i个顶点的法向索引
					}
					lindex++;
				}
				break;
			case 'g':
				// not really caring about faces or materials now
				// so just making life easier, I'll ignoring it
				break;
		}
	}

	fclose(fin);


	bool useNormals = !normalIndex.empty();
	bool useUVs = !uvIndex.empty();

	if (!useNormals && !useUVs) { 
		std::cout << "Copying points" << std::endl;
		// just copy the points into the array
	}

	useNormals = false;
	useUVs = false;
	if (!useNormals&&!useUVs)
	{
		std::cout<<"this mesh can only be used by MeshTriangle class"<<std::endl;
		vertextype = 0;
		verts = new Vector3[nverts];
		shifts = new float[nverts];
		for (unsigned int i=0; i < nverts; i++)
		{
			verts[i] = points[i];
		}
	} 
	else if(useNormals&&!useUVs)
	{
		std::cout<<"this mesh can only be used by MeshTriangleN class"<<std::endl;
		vertextype = 1;
		vertNs = new VertexN[nverts];
		shifts = new float[nverts];
		for (unsigned int i=0; i < nverts; i++)
		{
			vertNs[i].vertex = points[i];
			vertNs[i].normal = file_normals[i];
		}
	}
	else if(!useNormals&&useUVs)
	{
		std::cout<<"this mesh can only be used by MeshTriangleUV class"<<std::endl;
		vertextype = 2;
		vertUVs = new VertexUV[nverts];
		shifts = new float[nverts];
		for (unsigned int i=0; i < nverts; i++)
		{
			vertUVs[i].vertex = points[i];
			vertUVs[i].uv = file_uvvector[i];
		}
	}
	else if (useNormals&&useUVs)
	{
		std::cout<<"this mesh can only be used by MeshTriangleUVN class"<<std::endl;
		vertextype = 3;
		vertUVNs = new VertexUVN[nverts];
		shifts = new float[nverts];
		for (unsigned int i=0; i < nverts; i++)
		{
			vertUVNs[i].vertex = points[i];
			vertUVNs[i].normal = file_normals[i];
			vertUVNs[i].uv = file_uvvector[i];
			//std::cout<<i<<"#"<<vertUVNs[i].vertex.x()<<","<<vertUVNs[i].vertex.y()<<","<<vertUVNs[i].vertex.z()<<std::endl;
		}
	}
	else
	{
		vertextype = -1;
		std::cerr<<"what the hell is this obj!"<<std::endl;
		return false;
	}

	vertexIndex = new int[ntris*3];//三角形数目*3
	triAreas = new float[ntris];
	triDepths = new float[ntris];
	for (unsigned int i=0; i < ntris*3; i++)
	{	
		vertexIndex[i] = vertsIndex[i];	
	}

	points.clear();
	file_normals.clear();
	file_uvvector.clear();
	vertsIndex.clear();
	normalIndex.clear();
	uvIndex.clear();

	if (useNormals) std::cout << "Used normals" << std::endl;
	if (useUVs) std::cout << "Used texture coords" << std::endl;

	return true;
}


void Mesh::meshTransform( const Matrix & mx )
{
	if (vertextype == 0)
	{	
		for(unsigned int i=0;i<nverts;i++)
		{
			verts[i] = transformLoc(mx,verts[i]);
		}
	}

	if (vertextype == 3)
	{	
		for(unsigned int i=0;i<nverts;i++)
		{
			vertUVNs[i].vertex = transformLoc(mx,vertUVNs[i].vertex);
			vertUVNs[i].normal = transformVec(mx,vertUVNs[i].normal).getNormalized();
		}
	}

}


Mesh::Mesh( const Mesh & ms )
{
	mptr = ms.mptr;
	ntris = ms.ntris;
	nverts = ms.nverts;
	vertextype = ms.vertextype;
	smptime = ms.smptime;

	verts=NULL;
	vertNs=NULL;
	vertUVs=NULL;
	vertUVNs=NULL;

	switch (vertextype)
	{
	case 0:
		verts = new Vector3[nverts];
		for (unsigned int i = 0 ; i < nverts ; i++)
			verts[i] = ms.verts[i];
		break;
	case 1:
		vertNs = new VertexN[nverts];
		for (unsigned int i = 0 ; i < nverts ; i++)
			vertNs[i] = ms.vertNs[i];
		break;
	case 2:
		vertUVs = new VertexUV[nverts];
		for (unsigned int i = 0 ; i < nverts ; i++)
			vertUVs[i] = ms.vertUVs[i];
		break;
	case 3:
		vertUVNs = new VertexUVN[nverts];
		for (unsigned int i = 0 ; i < nverts ; i++)
			vertUVNs[i] = ms.vertUVNs[i];
		break;
	default:
		std::cerr<<"ERROR Mesh(const Mesh &)"<<std::endl;
	}
	
	vertexIndex = new int[ntris*3];

	for (unsigned int i = 0 ; i < ntris*3 ; i++)
	{
		vertexIndex[i] = ms.vertexIndex[i];
	}
	
}

void Mesh::calTexCoord( const char & type )
{	
	
	// z即 纵坐标 ty
	// 通过x,y的坐标 得到角度，除以2π也就得到横坐标tx
	float temp,temp2;
	float pai = 3.14159265358979323846f;
	float x,y,z;
	float phi , theta;
	
	switch (type)
	{
	case 's':
		for ( unsigned int i=0 ; i < nverts ; i++ )
		{

			x =  vertUVNs[i].vertex.x();
			y =  vertUVNs[i].vertex.y();
			z =  vertUVNs[i].vertex.z();
			

			if ( (z  > 0.999999999f) )
			{
				phi = 0.0f;
			}
			else if( (z  < -0.999999999f) )
			{
				phi = pai;
			}
			else
			{
				phi = acos(z);
			}
			phi /= pai;

			theta  = atan2f(y, x);
			(theta < 0.f) ? theta + 2.f*pai : theta;
			theta = theta / (2.0f*pai);

			vertUVNs[i].uv.setX(theta);
			vertUVNs[i].uv.setY(phi);

			//std::cout<<theta<<std::endl;
		}
		break;

	case 'c':
		for ( unsigned int i=0 ; i < nverts ; i++ )
		{
			x = vertUVNs[i].vertex.x();
			y = vertUVNs[i].vertex.y();
			z = vertUVNs[i].vertex.z();
			temp = acos(z);
			temp2 = acos(x/sqrt(1.0f-temp*temp)) ;
			if (y<0)//属于第3,4象限
				temp2 = 2*pai - temp2;
			std::cout<<temp<<","<<temp2<<std::endl;

			vertUVNs[i].uv.setX(temp2/(2*pai));
			vertUVNs[i].uv.setY(temp/pai);
			//std::cout<<vertUVNs[i].uv.x()<<","<<vertUVNs[i].uv.y()<<std::endl;
		}

		break;
	}
}

void Mesh::calVertexShift( const Matrix curmx ,const Matrix & mx )
{
	Matrix mm ;//= mx * curmx;
	Vector3 shiftedvert;
	Vector3 shiftedvert1, shiftedvert2;
	bool complicated = true;

	if (!complicated) // if it's not complicated
	{	
		for ( int i = 0 ; i < nverts ; i++)
		{
			shiftedvert = transformLoc( mx , verts[i] );
			shifts[i] = dist(verts[i] , shiftedvert);
		}
	}
	else // if it's complicated, divide the shift into several parts to approximate the shift.
	{
		for ( int i = 0 ; i < nverts ; i++ )
		{
			shifts[i] = 0.0f; // initial the shift
			shiftedvert1 = verts[i];
			for (int parts = 0 ; parts<20 ; parts++)
			{
				float angle = 60.0;
				//mm = rotate(Vector3(1,1,1) ,60.0*parts/20 );
				mm = rotateY(60*parts/20);
				shiftedvert2 = transformLoc( mm , verts[i] );
				shifts[i] += dist(shiftedvert1, shiftedvert2);
				shiftedvert1 = shiftedvert2;
			}

		}
	}

}

//面积S=根号(p*(p-a)*(p-b)*(p-c)) 其中p=(a+b+c)/2。
float triangleArea( Vector3 p1, Vector3 p2 , Vector3 p3)
{
	float a = dist(p1,p2);
	float b = dist(p2,p3);
	float c = dist(p1,p3);
	float p = (a+b+c)/2;
	return sqrt(p*(p-a)*(p-b)*(p-c));
	
}

void Mesh::calTriangleArea()
{
	float area;
	for (int i = 0 ; i<ntris ; i++)
	{
		area = triangleArea(verts[vertexIndex[3*i]], verts[vertexIndex[3*i+1]] ,verts[vertexIndex[3*i+2]] );
		triAreas[i] = area;
		triDepths[i] = (verts[vertexIndex[3*i]].z() + verts[vertexIndex[3*i+1]].z() + verts[vertexIndex[3*i+2]].z())/3.f;
		//printf("%d , <%d,%d,%d> = %f \n",i,vertexIndex[3*i] , vertexIndex[3*i+1] ,vertexIndex[3*i+2], area );
	}
}