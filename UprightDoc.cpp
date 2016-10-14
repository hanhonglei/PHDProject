#include "stdafx.h"
#include "MeshSimpDoc.h"
#include <stdmix.h>
#include <MxTimer.h>
#include "qslim.h"
#include "MeshSimp.h"
#include "MainFrm.h"
#include "MSDM2/MSDM2.h"
#include "ViewSelectionBenchmark/ViewSelect-OUR/ppm.h"
#include "Segment/segment-image.h"
#include "SegmentParam.h"
//#include "vld.h"
#include "MeshSimpView.h"
#include "gfx/gfx.h"
//////////////////////////////////////////////////////////////////////////
wchar_t *sh;

int offset4[]={-1,0,0,1,0,-1,1,0};
int offset[]={-1,-1,-1,0,-1,1,0,-1,0,1,1,-1,1,0,1,1};
// 判断是否输出每个候选视点及其深度图 [1/28/2014 Han]
#define  TEST_IMAGE_OUTPUT

// 使用8邻域的方法对输入的图像进行平滑滤波
void LaplacianFilterImage(int w, int h, unsigned char *normalData,unsigned char *depthData, int fType=8)
{
	unsigned char *normalDataF = new unsigned char [w*h*3];
	unsigned char *depthDataF = new unsigned char [w*h];
	int pixel[3], pixelDepth = 0;
	switch (fType)
	{
	case 8:
		for(int i=0;i<w;i++){
			for(int j=0;j<h;j++)
			{
				for (int ii = 0; ii < 3; ii++)
					pixel[ii] = 0;
				pixelDepth = 0;
				int num = 0;
				int p = j*w + i;
				for (int k=0;k<8;k++) 
				{ 
					// 向x方向移动dx，向y方向移动dy，移动的结果为(nx,ny)
					int nx = i + offset[k*2], ny = j + offset[k*2+1]; 
					int pN = ny*w + nx;
					// 判断(nx,ny)是不是在图像里
					if (0 <= nx && nx < w && 0 <= ny && ny < h ) 
					{
						num++;
						pixelDepth += depthData[pN];
						for (int ii = 0; ii < 3; ii++)
							pixel[ii] += normalData[pN*3+ii];
					}
				} 
				depthDataF[p] = (pixelDepth+depthData[p])/(num+1);
				for (int ii = 0; ii < 3; ii++)
					normalDataF[p*3+ii] =  (pixel[ii]+normalData[p*3+ii])/(num+1);
			}
		}
		break;
	default:
		break;
	}
	memcpy(depthData, depthDataF, w*h*sizeof(unsigned char));
	memcpy(normalData, normalDataF, w*h*3*sizeof(unsigned char));
	delete []normalDataF;
	delete []depthDataF;
}
// ↓From Hua [11/26/2013 Han]
void dfs4(int x, int y, int w, int h, unsigned char *depthData, int **mask, int res)
{
	mask[x][y] = res;
	// 测试 [11/27/2013 Han]
	depthData[y*w + x] = 0/*res%256*/;
	for (int k=0;k<4;k++) 
	{ 
		// 向x方向移动dx，向y方向移动dy，移动的结果为(nx,ny)
		int nx = x + offset4[k*2], ny = y + offset4[k*2+1]; 
		// 判断(nx,ny)是不是在图像里
		if (0 <= nx && nx < w && 0 <= ny && ny < h && mask[nx][ny] == -1 /*&& depthData[x*h+y] == depthData[nx*h+ny]*/) 
			dfs4(nx, ny,w,h,depthData,mask,res); 
	} 
}
// 使用8邻域的话一致度更高 [11/27/2013 Han]
void dfs(int x, int y, int w, int h, unsigned char *depthData, int **mask, int res)
{
	mask[x][y] = res;
	// 测试，将当前深度值设定为所属的区域值 [11/27/2013 Han]
	depthData[y*w + x] = 0/*res%256*/;
	for (int k=0;k<8;k++) 
	{ 
		// 向x方向移动dx，向y方向移动dy，移动的结果为(nx,ny)
		int nx = x + offset[k*2], ny = y + offset[k*2+1]; 
		// 判断(nx,ny)是不是在图像里
		if (0 <= nx && nx < w && 0 <= ny && ny < h && mask[nx][ny] == -1 /*&& depthData[x*h+y] == depthData[nx*h+ny]*/) 
			dfs(nx, ny,w,h,depthData,mask,res); 
	} 
}
// 采用非递归的方法，避免栈溢出
void DFS_NR(int x, int y, int w, int h,unsigned char *depthData, int **mask, int res)  
{  
	int top = 0;  
	mask[x][y] = res;
	// 测试，将当前深度值设定为所属的区域值 [11/27/2013 Han]
	depthData[y*w + x] = 0/*res%256*/;
	int *stackX = new int [w*h];
	int *stackY = new int [w*h];
	stackX[top] = x;
	stackY[top] = y;
	while ( top != -1)  
	{  
		x = stackX[top];  
		y = stackY[top];
		int k=0;
		for (k=0;k<8;k++) 
		{ 
			// 向x方向移动dx，向y方向移动dy，移动的结果为(nx,ny)
			int nx = x + offset[k*2], ny = y + offset[k*2+1]; 
			// 判断(nx,ny)是不是在图像里
			if (0 <= nx && nx < w && 0 <= ny && ny < h && mask[nx][ny] == -1 /*&& depthData[x*h+y] == depthData[nx*h+ny]*/) 
			{  
				mask[nx][ny] = res;
				depthData[ny*w + nx] = 0/*res%256*/;
				++top;
				stackX[top] = nx;
				stackY[top] = ny;
				break ;  
			}  
		}  
		if( k == 8)  
		{  
			top -- ;  
		}  
	}  
	delete []stackX;
	delete []stackY;
} 
//////////////////////////////////////////////////////////////////////////
// 返回值：分块数量
int LinkedParts(int w, int h, unsigned char *normalData,unsigned char *depthData
	,std::vector<vec3> &patchNormal,std::vector<int> &patchArea,float &cX, float &cY) 
{
	//对mask分配空间，存储是否遍历过该像素
	int **mask = new int *[w];
	for(int i=0;i<w;i++)
		mask[i] = new int[h];
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++)
		{
			int p = j*w+i;
			if(depthData[p]==255) // 该像素为背景
				mask[i][j]=-2;
			else mask[i][j]=-1;
		}
	}
	int res=0;
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++)
		{
			if(mask[i][j]==-1) //如果未遍历该像素
			{ 
				DFS_NR/*dfs*/(i, j,w,h,depthData,mask,res); 
				res++; 
			} 
		}
	}
	//找到n个不连通的图像块A1~An
	//	计算每个图像块的面积R1~Rn
	//	在normalMap中：
	//	计算A1~An图像块的平均法线N1~Nn
	//	计算A1~An图像块的质心Ca
	patchNormal.resize(res);
	float *patchCenterX = new float [res];
	float *patchCenterY = new float [res];
	patchArea.resize(res, 0);

	memset(patchCenterX, 0.f, (res)*sizeof(float));
	memset(patchCenterY, 0.f, (res)*sizeof(float));
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++)
		{
			if(mask[i][j]!=-1&&mask[i][j]!=-2) //如果该像素非背景像素的话
			{
				int p = j*w+i;
				if(normalData[p*3] != 0 || normalData[p*3+1]!=0 || normalData[p*3+2]!=0)// 不应该是背景色
				{
					(patchArea)[mask[i][j]]++;
					patchCenterX[mask[i][j]] += i;
					patchCenterY[mask[i][j]] += j;
					(patchNormal)[mask[i][j]].x += normalData[p*3];
					(patchNormal)[mask[i][j]].y += normalData[p*3+1];
					(patchNormal)[mask[i][j]].z += normalData[p*3+2];
				}
			}
		}
	}
	cX = 0.f, cY = 0.f;
	for (int i = 0; i < res; i++)
	{
		if (patchArea[i] > 0)
		{
			(patchNormal)[i].x /= (patchArea)[i];
			(patchNormal)[i].y /= (patchArea)[i];
			(patchNormal)[i].z /= (patchArea)[i];
			patchCenterX[i] /= (patchArea)[i];
			patchCenterY[i] /= (patchArea)[i];
			cX += patchCenterX[i];
			cY += patchCenterY[i];
		}
		else 
			res--;
	}
	// 计算质心
	cX /= res;
	cY /= res;
	delete []patchCenterX;
	delete []patchCenterY;
	for(int i=0;i<w;i++)
		delete []mask[i] ;
	delete []mask;
	return res;
	//如果最终res=1，则说明只有一个封闭区域
}
// ↑FromHua
//////////////////////////////////////////////////////////////////////////
/*
	Input 候选底部法线方向V，投影深度图depthMap，模型投影法线图normalMap
	Output 支撑性得分S，稳定性得分H
	在depthMap中
		找到投影图像中心点C
		计算和C最远的像素点距离D
		计算图像区域的面积R
		删除深度大于Z的所有像素
		找到n个不连通的图像块A1~An
		计算每个图像块的面积R1~Rn
	在normalMap中：
		计算A1~An图像块的平均法线N1~Nn
		计算A1~An图像块的质心Ca
	支撑性得分：S=∑_n▒((V∙N_i ) R_i/R) 
	稳定性得分：H=|Ca-C|/D
*/
// 计算候选视点的得分,包括是否接近四分之三,是否有一个平整地面,是否有多个支撑脚 [11/11/2013 Han]
// normalData:黑色表示背景,rgb表示全局空间的法线值
// depthData:白色表示背景,越接近白色说明深度越大
int LegAndBase(int w, int h, unsigned char *normalData,unsigned char *depthData, float *baseNormal, float depthPerc, float &S, float &H)
{
	float cx = 0.f, cy = 0.f;	// 图像中心点坐标
	int area = 0;				// 图像中有效像素的数目
	int zMax = 0, zMin = 255;	// 找到图像中的最大和最小深度值
	// 找到中心点坐标和图像面积
	for(int i=0;i<w;i++)
		for(int j=0;j<h;j++)
		{
			int p = j*w + i;
			if(depthData[p]!=255) // 该像素非背景
			{ 
				cx += i;
				cy += j;
				area++;
				if (depthData[p] < zMin)
					zMin = depthData[p];
				if (depthData[p] > zMax)
					zMax = depthData[p];
			} 
		}
	cx /= area;
	cy /= area;

	// 1. 找到离中心点最远的像素距离,并将深度大于Z的像素标记为255
	float dMax = 0.f;
	float zThreshold = zMin + (zMax - zMin) * depthPerc;
	for(int i=0;i<w;i++)
		for(int j=0;j<h;j++)
		{
			int p = j*w + i;
			if(depthData[p]!=255) //该像素非背景
			{ 
				float dx = i - cx;
				float dy = j - cy;
				float d = dx*dx + dy*dy;
				if( d > dMax )
					dMax = d;
				// 删除那些深度离得较远的像素
				if (depthData[p] > zThreshold)
					depthData[p] = 255;
			} 
		}
	// 2. 计算剩余的像素图中有多少不连通的分块，这些分块即是支撑底部
	std::vector<vec3> patchNormal;
	std::vector<int> patchArea;
	float patchCX, patchCY;
	int patchN = LinkedParts(w, h, normalData, depthData, patchNormal, patchArea, patchCX, patchCY);
	
	// 3. 计算候选底面的得分
	//支撑性得分：S=∑_n▒((V∙N_i ) R_i/R) 
	//稳定性得分：H=|Ca-C|/D	
	S = 0.f;
	vec3 viewV = baseNormal;
	viewV.normalize();
	for (int i = 0; i < patchN; i++)
	{
		patchNormal[i].normalize();
		S += /*abs(patchNormal[i].dot(viewV)) * */patchArea[i] / float(area);
	}

	H = 1.0-((patchCX - cx)*(patchCX - cx) + (patchCY - cy)*(patchCY - cy))/dMax;

	patchArea.clear();
	patchNormal.clear();
	return patchN;
}

// @in candBase 候选底部的列表，每个底部视点已经经过了归一化
int* CMeshSimpDoc::CandBaseScore(const vector<vec3*> &candBase, const vec3 &bestView, int &bestBaseIdx)
{
	// 以下作为评分：
	// 1.视点舒适度（C）：选择和最优视点最接近四分之三视角的作为最终的底面,即两个向量的夹角为180*3/4=135，cos(135)=-0.7
	// 2.稳定性（S）：模型底部往往有比较平稳的底面，比如人造物品
	// 3.支撑性（H）:模型底部往往具有足部(2个或以上)，
	float depthPerc = 2/100.0; // play with this number: 2/100.0;
	int nCands = candBase.size();	
	float *C = new float[nCands];
	float *S = new float[nCands];
	float *H = new float[nCands];
	// 升序排列，依靠键值对
	std::multimap<float, int> sortedC;
	std::multimap<float, int> sortedS;
	std::multimap<float, int> sortedH;
	std::multimap<int, int> sortedScores;
	int *scores = new int[nCands];
	memset(scores, 0, nCands * sizeof(int));

	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 
	CMeshSimpView::RenderType stb = pView->showing_type;
	pView->showing_type = CMeshSimpView::SHOWING_NORMALMAP;
	glShadeModel(GL_SMOOTH); 
	float d = pView->GetViewDistance();

	GLint viewport[4];											//space for viewport data
	glGetIntegerv (GL_VIEWPORT, viewport);						//get actual viewport
	int w = viewport[2], h = viewport[3];

	unsigned char *normalData = new unsigned char [w*h*3];
	unsigned char *depthData = new unsigned char [w*h];
	float baseNormal[3];
#ifdef TEST_IMAGE_OUTPUT
	PPM *ppmWriter = new PPM;
#endif	
	for (int i = 0; i < nCands; i++)
	{
		C[i] = 1 - abs(candBase[i]->dot(bestView)+0.7f)/1.7f;

		baseNormal[0] = -candBase[i]->x;
		baseNormal[1] = -candBase[i]->y;
		baseNormal[2] = -candBase[i]->z;

		pView->eye.x = candBase[i]->x;
		pView->eye.y = candBase[i]->y;
		pView->eye.z = candBase[i]->z;
		pView->eye = pView->eye*d;
		pView->AdjustView();

		pView->DisplayModel();
		glFinish();
		//ust change the packing to ensure no overruns!
		glPixelStorei(GL_PACK_ALIGNMENT, 1);

		// OpenGL以左下角为坐标原点，而windows系统是左上角
		glReadPixels(0, 0, w, h, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, depthData);
#ifdef TEST_IMAGE_OUTPUT
		//////////////////////////////////////////////////////////////////////////
		// Test 测试 用后删除，检验输出的法线图和深度图是否正确 [11/11/2013 Han]
		ppmWriter->width = w;
		ppmWriter->height = h;
		ppmWriter->version = "P6";
		ppmWriter->data = normalData;
		wchar_t fnRe[256];
		//swprintf( fnRe,   L"OutPutFiles\\normalMap_%hs_%hs_Cand_%d_%d_%d.ppm", sh, m_sModelName, i+1, nCands, rand());
		//ppmWriter->save(fnRe);
		//for (int xxx = 0; xxx < w*h*3; xxx+=3)
		//{
		//	normalData[xxx] = normalData[xxx+1] = normalData[xxx+2] = depthData[xxx/3];
		//}

		// test 暂时不输出深度图,节约空间 [12/22/2013 Han]
		//swprintf( fnRe,   L"OutPutFiles\\depthMapBefore_%s_%hs_Cand_%d-%d_%d.ppm", sh,m_sModelName, i+1, nCands, rand());
		//ppmWriter->save(fnRe);
		//////////////////////////////////////////////////////////////////////////
#endif
		glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, normalData);

		// 平滑滤波 [11/27/2013 Han]
		LaplacianFilterImage(w, h, normalData, depthData, 8);

#ifdef TEST_IMAGE_OUTPUT
		//swprintf( fnRe,   L"OutPutFiles\\normalMapSmoothed_%s_Cand_%d_%d_%d.ppm", m_sModelName, i+1, nCands, rand());
		//ppmWriter->save(fnRe);
		for (int xxx = 0; xxx < w*h*3; xxx+=3)
		{
			normalData[xxx] = normalData[xxx+1] = normalData[xxx+2] = depthData[xxx/3];
		}
		swprintf( fnRe,   L"OutPutFiles\\depthMapBeforeSmoothed_%s_%hs_Cand_%d-%d_%d.ppm", sh, m_sModelName, i+1, nCands, rand());
		ppmWriter->save(fnRe);
#endif
		// 通过候选底部得到的深度图和法线图计算这个候选底部的得分，进而得到最终的底部
		int patchNum = LegAndBase(w, h, normalData, depthData, baseNormal, depthPerc, S[i], H[i]);	
		fprintf(outputFile, "Cand index:%d\tBase Num: %d\tC:%f;\tS:%f;\tH:%f\n", i, patchNum, C[i], S[i], H[i]);
		sortedC.insert(make_pair(C[i], i));
		sortedS.insert(make_pair(S[i], i));
		sortedH.insert(make_pair(H[i], i));

#ifdef TEST_IMAGE_OUTPUT			
		// Test 测试 用后删除，检验输出的法线图和深度图是否正确 [11/11/2013 Han]
		ppmWriter->width = w;
		ppmWriter->height = h;
		ppmWriter->version = "P6";
		ppmWriter->data = normalData;
		for (int xxx = 0; xxx < w*h*3; xxx+=3)
		{
			normalData[xxx] = normalData[xxx+1] = normalData[xxx+2] = depthData[xxx/3];
		}
		swprintf( fnRe,   L"OutPutFiles\\depthMapAfter_%s_%hs_Cand_%d-%d_%d.ppm", sh, m_sModelName, i+1, nCands, rand());
		ppmWriter->save(fnRe);
#endif
	}

	std::multimap<float, int>::iterator  ItC = sortedC.begin(), ItS = sortedS.begin(), ItH = sortedH.begin();
	int Ci = 0, Si = 0, Hi = 0;
	int Cs, Ss, Hs;
	for (int i = 0; i < nCands; i++)
	{
		scores[ItC->second] += Ci;
		scores[ItS->second] += Si;
		scores[ItH->second] += Hi;
		// 取得分到小数点后第三位，得分相同，则即为并列 [12/4/2013 Han]
		Cs = (ItC->first+0.005) * 100;
		Ss = (ItS->first+0.005) * 100;
		Hs = (ItH->first+0.005) * 100;
		ItC++, ItS++, ItH++;
		if (ItC == sortedC.end())	// 如果已经越界，则直接退出
			break;
		if (Cs != int((ItC->first+0.005) * 100))
			Ci++;
		if (Ss != int((ItS->first+0.005) * 100))
			Si++;
		if (Hs != int((ItH->first+0.005) * 100))
			Hi++;
	}
	for (int i = 0; i < nCands; i++)
	{
		sortedScores.insert(make_pair(scores[i], i));
		fprintf(outputFile, "Candindex:%d_Score:%d;\t", i, scores[i]);
	}
	std::multimap<int, int>::reverse_iterator  It = sortedScores.rbegin();
	bestBaseIdx = It->second; // 得分最高的视点即为最终的底部，然后如果其他视点和其得分相同，则选择本来在前的视点作为最终底部
	int nTBestScore = It->first;
	It++;
	for (; It != sortedScores.rend(); It ++)
	{	
		if (It->first == nTBestScore)
		{
			if(It->second < bestBaseIdx)
				bestBaseIdx = It->second;
		}
		else
			break;
	}
	fprintf(outputFile, "\n");
	delete []C;
	delete []S;
	delete []H;
	delete []normalData;
	delete []depthData;
	pView->showing_type = stb ;
#ifdef TEST_IMAGE_OUTPUT	
	ppmWriter->data = 0;
	delete ppmWriter;
#endif
	return scores;
}

//typedef struct {
//	float theta,phi,score;
//} XYZ;
//typedef struct {
//	XYZ p1,p2,p3;
//} FACET3;
//
//XYZ MidPoint(XYZ,XYZ);
//void Normalise(XYZ *);
//float CMeshSimpDoc::RefineCandBase(ViewSelectType FilterType, vec3 view, vec3 &refinedView)
//{
//	int i,j;
//	int n=0,nstart;
//	int iterations = 2;
//	FACET3 *f = NULL;
//	double theta[3] = {0.0,35.0,80.0}, phi[3] = {10.0,15.0,80.0}; // corner in polar coordinates
//	XYZ p1,p2,p3;
//
//	// Start with the vertices of the triangle
//	f = malloc(sizeof(FACET3));
//	f[0].p1.x = cos(phi[0]*DTOR) * cos(theta[0]*DTOR);
//	f[0].p1.y = cos(phi[0]*DTOR) * sin(theta[0]*DTOR);
//	f[0].p1.z = sin(phi[0]*DTOR);
//	f[0].p2.x = cos(phi[1]*DTOR) * cos(theta[1]*DTOR);
//	f[0].p2.y = cos(phi[1]*DTOR) * sin(theta[1]*DTOR);
//	f[0].p2.z = sin(phi[1]*DTOR);
//	f[0].p3.x = cos(phi[2]*DTOR) * cos(theta[2]*DTOR);
//	f[0].p3.y = cos(phi[2]*DTOR) * sin(theta[2]*DTOR);
//	f[0].p3.z = sin(phi[2]*DTOR);
//	n = 1;
//
//		for (j=0;j<nstart;j++) {
//			f = realloc(f,(n+3)*sizeof(FACET3));
//
//			// Create initially copies for the new facets 
//			f[n  ] = f[j];
//			f[n+1] = f[j];
//			f[n+2] = f[j];
//
//			// Calculate the midpoints 
//			p1 = MidPoint(f[j].p1,f[j].p2);
//			Normalise(&p1);
//			p2 = MidPoint(f[j].p2,f[j].p3);
//			Normalise(&p2);
//			p3 = MidPoint(f[j].p3,f[j].p1);
//			Normalise(&p3);
//
//			// Replace the current facet 
//			f[j].p2 = p1;
//			f[j].p3 = p3;
//
//			// Create the changed vertices in the new facets 
//			f[n  ].p1 = p1;
//			f[n  ].p3 = p2;
//			f[n+1].p1 = p3;
//			f[n+1].p2 = p2;
//			f[n+2].p1 = p1;
//			f[n+2].p2 = p2;
//			f[n+2].p3 = p3;
//			n += 3;
//		}
//}
// ***将得到的候选底部进行refine，避免离散视点带来的底部偏差问题 [1/28/2014 Han]
// 具体做法是:在这些候选底部的一定邻居范围内(小于最近的离散邻居)进行细分,然后对细分后的视点接着打分,选取其中最小的视点作为refine后的位置
// 暂未实现,考虑使用曲面拟合的方法进行 [1/28/2014 Han]
// 输出最差视点及其周围六个视点 [2/14/2014 Han]
// @vec3 *candBase:待调整位置的视点
// @int baseID:待调整视点在观察球模型中的ID
// @ViewSelectType FilterType:当前所使用的视点评分方法
float CMeshSimpDoc::Subdivide(vec3 *candBase,int baseID, ViewSelectType FilterType, bool bSmallest)
{
	bool bCur = true;
	int nF = -1;
	// 1 对最差视点的邻域进行细分,并求取细分后的视点质量及1 ring邻居的质量,将其作为控制点
	MxVertex &vCurr = m_pViewPhereMdl->vertex(baseID);		// 当前顶点 [4/18/2012 Han]
	vec3 P = vec3(vCurr.as.pos);
	MxVertexList star;
	star.reset();
	m_pViewPhereMdl->collect_vertex_star(baseID, star);
	fprintf(outputFile, "Worst point original:%f\t%f\t%f\t\t%f\n", vCurr[0], vCurr[1],vCurr[2], vCurr.view_importance);
#ifdef TEST_IMAGE_OUTPUT	
	fprintf(outputFile, "Worst ID:%d:neighbour num=%d\n", baseID, star.length());	// 观察球的邻居至少有4个,一般是6个
#endif
	float Theta = acos(P.z);
	float Phi = atan(P.y / P.x);
#ifdef TEST_IMAGE_OUTPUT			
	fprintf(outputFile, "%f\t%f\t%f\t(Theta:%f,Phi:%f)\t\t%f\n", vCurr[0], vCurr[1],vCurr[2], Theta, Phi,vCurr.view_importance);
#endif
	MxVertex extremPoint = vCurr;		// 将极值点设置为当前点
	extremPoint.view_importance = vCurr.view_importance;
	vector<vec3> ctrlPoints;
	ctrlPoints.push_back(vec3(Theta, Phi, vCurr.view_importance));
	for(int j = 0; j < star.length(); j++)
	{
		MxVertex &vThat = m_pViewPhereMdl->vertex(star(j));
		float Theta = acos(vThat[2]);
		float Phi = atan(vThat[1] / vThat[0]);
#ifdef TEST_IMAGE_OUTPUT			
		fprintf(outputFile, "%f\t%f\t%f\t(Theta:%f,Phi:%f)\t\t%f\n", vThat[0],vThat[1],vThat[2], Theta, Phi, vThat.view_importance);
#endif
		ctrlPoints.push_back(vec3(Theta, Phi, vThat.view_importance));
	}
	// a 使用邻接边中点细分
	/*
	fprintf(outputFile, "Refine edge center:\n");
	for(int j = 0; j < star.length(); j++)
	{
	vec3 p = vec3(m_pViewPhereMdl->vertex_position(star(j)));
	vec3 c = 0.5*(P+p);
	c.normalize();
	refinedV[i][0] = c.x;
	refinedV[i][1] = c.y;
	refinedV[i][2] = c.z;
	refinedV[i].view_importance = ViewQuality(c, FilterType);
	if (refinedV[i].view_importance < smallestPoint.view_importance)	// 当前点的质量最低
	smallestPoint = refinedV[i];
	Theta = acos(c.z);
	Phi = atan(c.y / c.x);
	fprintf(outputFile, "%f\t%f\t%f\t(Theta:%f,Phi:%f)\t\t%f\n", c.x, c.y, c.z, Theta, Phi,refinedV[i].view_importance);
	ctrlPoints.push_back(vec3(Theta, Phi, refinedV[i].view_importance));

	}*/
	// b 使用邻接三角形中点细分
#ifdef TEST_IMAGE_OUTPUT			
	fprintf(outputFile, "Refine face center:\n");
#endif
	const MxFaceList &faceN = m_pViewPhereMdl->neighbors(baseID);	// 一定要使用const，否则出错 [2/17/2014 Han]
	MxVertex *refinedV = new MxVertex[faceN.length()]; 
	for(int j = 0; j < faceN.length(); j++)
	{
		float x = 0, y = 0, z = 0;
		vec3 c;
		c.setZero();
		for (int mm = 0; mm < 3; mm++)
		{
			c += vec3( m_pViewPhereMdl->vertex_position(m_pViewPhereMdl->face(faceN(j))(mm)));
		}
		c = c / 3.f;
		c.normalize();
		refinedV[j][0] = c.x;
		refinedV[j][1] = c.y;
		refinedV[j][2] = c.z;
		refinedV[j].view_importance = ViewQuality(c, FilterType);
		if ((bSmallest && refinedV[j].view_importance < extremPoint.view_importance)	// 当前点的质量最低
			|| (!bSmallest && refinedV[j].view_importance > extremPoint.view_importance))
		{
			extremPoint = refinedV[j];
			extremPoint.view_importance = refinedV[j].view_importance;
			bCur = false;
			nF = j;
		}
		Theta = acos(c.z);
		Phi = atan(c.y / c.x);
#ifdef TEST_IMAGE_OUTPUT
		fprintf(outputFile, "%f\t%f\t%f\t(Theta:%f,Phi:%f)\t\t%f\n", c.x, c.y, c.z, Theta, Phi,refinedV[j].view_importance);
#endif
		ctrlPoints.push_back(vec3(Theta, Phi, refinedV[j].view_importance));
	}
	// 2 对这些控制点进行拟合,求取拟合曲面的极小值点
	fprintf(outputFile, "Worst point refined:%f\t%f\t%f\t\t%f\n", extremPoint[0], extremPoint[1],extremPoint[2],extremPoint.view_importance);
	if (bCur)		// 如果当前点是极值点
	{
		for(int j = 0; j < faceN.length(); j++)	{
			vec3 c;
			c.setZero();
			c += vec3(refinedV[j].as.pos) + P;	// 取当前点和邻域三角形中点的中心作为细分点
			c = c / 2.f;
			c.normalize();	
			float vq = ViewQuality(c, FilterType);
			if ((bSmallest && vq < extremPoint.view_importance)	// 当前点的质量最低
				|| (!bSmallest && vq > extremPoint.view_importance))
			{
				extremPoint[0] = c.x;
				extremPoint[1] = c.y;
				extremPoint[2] = c.z;
				extremPoint.view_importance = vq;
			}
			Theta = acos(c.z);
			Phi = atan(c.y / c.x);
#ifdef TEST_IMAGE_OUTPUT			
			fprintf(outputFile, "%f\t%f\t%f\t(Theta:%f,Phi:%f)\t\t%f\n", c.x, c.y, c.z, Theta, Phi,vq);
#endif
		}

	}
	else			// 如果细分点是极值点
	{
		ASSERT(nF != -1);
		vec3 c = vec3(extremPoint.as.pos);
		vec3 v = vec3(m_pViewPhereMdl->vertex_position(m_pViewPhereMdl->face(faceN(nF))(0)));
		vec3 *nvs = new vec3[6];
		for (int j = 0; j < 3; j++)
		{
			vec3 d = (c+v)/2.f;					// 中点
			d.normalize();	
			nvs[j*2] = d;

			v = vec3(m_pViewPhereMdl->vertex_position(m_pViewPhereMdl->face(faceN(nF))((j+1)%3)));
			vec3 f = (c+v)/2.f;
			f.normalize();
			f = (f+d)/2.f;
			f.normalize();
			nvs[j*2+1] = f;
		}
		for (int j = 0; j < 6; j++)
		{
			float vq = ViewQuality(nvs[j], FilterType);
			if ((bSmallest && vq < extremPoint.view_importance)	// 当前点的质量最低
				|| (!bSmallest && vq > extremPoint.view_importance))
			{
				extremPoint[0] = nvs[j].x;
				extremPoint[1] = nvs[j].y;
				extremPoint[2] = nvs[j].z;
				extremPoint.view_importance = vq;
			}
			Theta = acos(nvs[j].z);
			Phi = atan(nvs[j].y / nvs[j].x);
#ifdef TEST_IMAGE_OUTPUT			
			fprintf(outputFile, "%f\t%f\t%f\t(Theta:%f,Phi:%f)\t\t%f\n", nvs[j].x, nvs[j].y, nvs[j].z, Theta, Phi,vq);
#endif
		}
		delete []nvs;
	}
	fprintf(outputFile, "Worst point refined2:%f\t%f\t%f\t\t%f\n", extremPoint[0], extremPoint[1],extremPoint[2],extremPoint.view_importance);

	// 3 将你和曲面极小值点替换原来的最差视点作为候选底部,这样消除离散观察球带来的误差
	candBase->x = extremPoint[0];
	candBase->y = extremPoint[1];
	candBase->z = extremPoint[2];

	// test delete 输出调整后的最差视点 [7/25/2014 Han]
	ViewQuality(*candBase, FilterType, true, false);

	delete []refinedV;
	return extremPoint.view_importance;

}
// ↓计算模型竖直方向[10/9/2013 Han]
// Get best and some (3 ?) worst viewpoints from center points of clustered viewpoints
// Compair the angle between best viewpoint and worst viewpoints (include their opposite viewpoints)
// The worst viewpoint whose angle is more close to 3/4
ViewPoint* CMeshSimpDoc::BestBaseView(float k, vector<vec3*> &candBase, int &nBestBase,ViewSelectType FilterType)
{
	bool bUseClusteredBestView = false/*true*/;	// 是否使用聚类后的中心点作为最优视点，还是仍然使用直接计算的离散视点
	float minAngleCos = 0.9;		// 接近180°

	std::multimap<float, float *> sortedViews;		// 按照得分从小到大排序经过聚类的视点中心
	std::multimap<float, int> IDs;
	map<int,struct Seg>::iterator it;
	for ( it=m_pViewPhereMdl->saliencySegData.begin() ; it != m_pViewPhereMdl->saliencySegData.end(); it++ )
	{
		sortedViews.insert(make_pair(it->second.Saliency, it->second.pos));
		// 输出分块以后的观察球信息 [2/16/2014 Han]
		//fprintf(outputFile, "Seg scoring:%f:num=%d\n", it->second.Saliency, it->second.num);
		IDs.insert(make_pair(it->second.Saliency, it->second.ID));
	}
	
	// 1 选择视点作为候选底面视点,要求视点不能接近平行,而且得分都很低
	std::multimap<float, float *>::iterator  It = sortedViews.begin();
	std::multimap<float, int>::iterator ItID = IDs.begin();
	vector<int> sortedIDs;
	sortedIDs.push_back(ItID->second);
	candBase.push_back(new vec3(It->second));
	candBase[candBase.size()-1]->normalize();		// 加入候选底面视角的向量都是单位向量,便于后续的点乘运算
	for (++It, ++ItID; (It != sortedViews.end() /*&& nn < 3*/); It++, ItID++)	// 如果其他候选视点的得分比较接近，并且和已有候选视点夹角足够大话，将他们作为候选底面
	{	
		if (It->first-sortedViews.begin()->first < k)
		{
			bool bCand = true;
			vec3 *vecT = new vec3(It->second);
			vecT->normalize();
			for (int i = 0; i < candBase.size(); i++)			{
				if (abs(candBase[i]->dot(*vecT)) > minAngleCos)				{
					bCand = false;
					break;
				}
			}
			if (bCand)	{
					candBase.push_back(vecT);
					sortedIDs.push_back(ItID->second);
			}
			else
				delete vecT;
		}
		else
			break;
	}
	// ***将得到的候选底部进行refine，避免离散视点带来的底部偏差问题 [1/28/2014 Han]
	for (int i = 0; i < candBase.size(); i++)
	{
		Subdivide(candBase[i], sortedIDs[i], FilterType);
	}
	fflush(outputFile);
	// ***

	// 2 将候选底面的对向也加入候选底面中
	int nn = candBase.size();
	for (int i = 0; i < nn; i++)
		candBase.push_back(new vec3(-*candBase[i]));


	// 3 求每个候选底面的得分，按照文章要求：一、是否接近四分之三；二、是否是平整地面；三、是否有多个支撑脚
	nBestBase = 0;
	vec3 bestView = bUseClusteredBestView ? vec3(sortedViews.rbegin()->second)
		:vec3(viewpoint_candidates[m_nBestViewID]->pos[0], viewpoint_candidates[m_nBestViewID]->pos[1],viewpoint_candidates[m_nBestViewID]->pos[2]);
	bestView.normalize();
	if (m_clusterdBestView==NULL)
		m_clusterdBestView = new ViewPoint;
	// 保存经过聚类以后的最优视点
	vec3 vt = vec3(sortedViews.rbegin()->second);
	vt.normalize();
	m_clusterdBestView->pos[0] = vt.x;
	m_clusterdBestView->pos[1] = vt.y;
	m_clusterdBestView->pos[2] = vt.z;

	// 得到每个候选底部的得分，然后从中选出最终的视点
	double ts = 0.0;
	int *candScores;
	TIMING(ts, candScores = CandBaseScore(candBase, bestView, nBestBase));
	fprintf(outputFile, "Select final base time:%lfs\n", ts);

	if (m_baseViewpoint == NULL)
		m_baseViewpoint = new ViewPoint;
	m_baseViewpoint->pos[0] = candBase[nBestBase]->x;
	m_baseViewpoint->pos[1] = candBase[nBestBase]->y;
	m_baseViewpoint->pos[2] = candBase[nBestBase]->z;

	delete []candScores;
	return m_baseViewpoint;
}

void CMeshSimpDoc::Upright(ViewSelectType FilterType)
{
	double tu=get_cpu_time();
	// 1 Score viewpoints
	switch (FilterType)
	{
	case CURVATURE_ENTROPY:
		fprintf(outputFile, "-------------------------------------------------\n☆ Upright Using Curvature Entropy Scoring Method.\n");
		OnViewpointselectionCe();
		sh = L"CE";
		break;
	case SALIENCY:
		fprintf(outputFile, "-------------------------------------------------\n☆ Upright Using Mesh Saliency Scoring Method.\n");
		OnEditSaliencyvs();
		sh = L"MS";
		break;
	case MEAN_CURVATURE_ENTROPY:
		fprintf(outputFile, "-------------------------------------------------\n☆ Upright Using Our Mean Curvature + Entropy Scoring Method.\n");
		OnViewpointselectionMeancentropy();
		sh = L"MCE";
		break;
	case MAX_AREA:
		fprintf(outputFile, "-------------------------------------------------\n☆ Upright Using Max Area Scoring Method.\n");
		OnViewpointselectionMaxarea();
		sh = L"MA";
		break;
	case MEAN_CURVATURE:
		fprintf(outputFile, "-------------------------------------------------\n☆ Upright Using Mean Curvature Scoring Method.\n");
		OnViewpointselectionMeancurvature();
		sh = L"MC";
		break;
	case VIEW_ENTROPY:
		fprintf(outputFile, "-------------------------------------------------\n☆ Upright Using View Entropy Scoring Method.\n");
		OnViewpointselectionViewentropy();
		sh = L"VE";
		break;
	case BASE_CURVATURE:
		fprintf(outputFile, "-------------------------------------------------\n☆ Upright Using Base Curvature Scoring Method.\n");
		ViewpointSelectionBaseCurvature();
		sh = L"BC";
		break;
	default:
		break;
	}	
	//ASSERT(m_nBestViewID != -1);	// 有时候，模型会得不到最优视点，即视点评分失败，这种情况下，直接返回，不进行进一步计算
	if (m_nBestViewID == -1)
	{
		fprintf(outputFile, "!! There are not exisit bestviewpoints.\n");
		return;
	}

	UpdateViewsImportance();// Copy importance to viewsphere model

	// 2 Candidate viewpoints clustering
	float nSegNum = 20.f;
	int n = 0;
	double ts =0.0;
	int min_size = 1;	float k = (m_pViewPhereMdl->MaxImportanceV - m_pViewPhereMdl->MinImportanceV) / nSegNum;
	TIMING(ts, segment_image(m_pViewPhereMdl,  sigma
		, k
		, min_size , &n, false));
	fprintf(outputFile, "Clustering viewpoints Time:%lf seconds, Segment Params:sigma=%f,k=%f,min=%d---->viewpoints num:%d-Segment Number:%d\n"
		,ts, sigma, k,min_size, m_pViewPhereMdl->vert_count(),n);

	m_pViewPhereMdl->saliencySegData.clear();
	m_pViewPhereMdl->UpdateSaliencyBySegment(false, false);
	UpdateViewsImportance(false);	// Copy seged importance from model to viewsphere

	// 3 Get best and some (3 ?) worst viewpoints from clustered viewpoints
	// Compair the angle between best viewpoint and worst viewpoints (include their opposite viewpoints)
	// The worst viewpoint whose angle is more close to 3/4
	vector<vec3*> candBase;
	int nBestBase = 0;
	BestBaseView(k, candBase, nBestBase, FilterType);

	tu=get_cpu_time() - tu;
	fprintf(outputFile, "Uprighting using Total time:%lf seconds;Base Viewpoint: %f %f %f; Clustered best viewpoints:%f %f %f\n-------------------------------------------------\n"
		,tu, m_baseViewpoint->pos[0], m_baseViewpoint->pos[1],m_baseViewpoint->pos[2]
	,m_clusterdBestView->pos[0], m_clusterdBestView->pos[1], m_clusterdBestView->pos[2] );

	// 将聚类后的最优视点按照base在下的原则进行截图 [10/11/2013 Han]
	CMainFrame   *pMain=(CMainFrame *)AfxGetApp()->m_pMainWnd; 
	CMeshSimpView   *pView=(CMeshSimpView   *)pMain-> GetActiveView(); 
	pView->showing_type = pView->SHOWING_MODEL;
	pView->bShowUprightDir = true;
	void SaveRenderPic(const point *viewpos,float d, wchar_t* fn, CMeshSimpView *pView);
	wchar_t fn[256];
//#define USE_CLUSTERED_BEST_VIEW	// 是否使用cluster的最优视点的话 [10/12/2013 Han]	

#ifdef USE_CLUSTERED_BEST_VIEW
	swprintf( fn,   L"☆%s_UprightedClusteredBestView_%hs_%d.ppm ", sh, m_sModelName, rand());
	SaveRenderPic(&m_clusterdBestView->pos, pView->GetViewDistance(),  fn, pView);
#endif
	swprintf( fn,   L"☆%s_UprightedOrigBestView_%hs_%d_%d.ppm ", sh, m_sModelName, nBestBase, rand());
	SaveRenderPic(&viewpoint_candidates[m_nBestViewID]->pos, pView->GetViewDistance(),  fn, pView);

	point tmpBase = m_baseViewpoint->pos;
	int nn = candBase.size();
#ifdef TEST_IMAGE_OUTPUT
	for (int i = 0; i < nn; i++)
	{
			swprintf( fn,   L"☆%s_CandBaseView_%hs_%d-%d_%d.ppm ", sh,m_sModelName,i+1,nn,  rand());
			m_baseViewpoint->pos[0] = candBase[i]->x;
			m_baseViewpoint->pos[1] = candBase[i]->y;
			m_baseViewpoint->pos[2] = candBase[i]->z;

			//应该使用bestview进行绘制
			SaveRenderPic(&viewpoint_candidates[m_nBestViewID]->pos, pView->GetViewDistance(),  fn, pView);
	}	
	bool tt = pView->bShowUprightDir;
	ViewPoint* ttt = m_baseViewpoint;
	pView->bShowUprightDir = false;
	m_baseViewpoint = NULL;
	swprintf( fn,   L"☆%s_BestView_%hs_%d.ppm ", sh,m_sModelName, rand());
	//应该使用bestview进行绘制
	SaveRenderPic(&viewpoint_candidates[m_nBestViewID]->pos, pView->GetViewDistance(),  fn, pView);
	pView->bShowUprightDir = tt;
	m_baseViewpoint = ttt;
#endif

	m_baseViewpoint->pos = tmpBase;
	fprintf(outputFile, "There are %d dandidate base viewpoints:\n", candBase.size());
	for (int i = 0; i < nn; i++)
	{
		if (nBestBase == i)
			fprintf(outputFile, "The Best Base-->");
		fprintf(outputFile, "%f, %f, %f\n"
			,candBase[i]->x, candBase[i]->y, candBase[i]->z);
		fflush(outputFile);
		delete candBase[i];
	}
	candBase.clear();
	pView->AdjustView();
	pView->Invalidate(true);
}
// ↑计算模型竖直方向[10/9/2013 Han]
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// ↓使用的各种视点评分方法
void CMeshSimpDoc::OnUprightCurshannon()
{
	Upright(CURVATURE_ENTROPY);
}

// 使用Our的curvature+shannon方法进行视点选择
void CMeshSimpDoc::OnUprightMcs()
{
	Upright(MEAN_CURVATURE_ENTROPY);
}


void CMeshSimpDoc::OnUprightMeancur()
{	
	Upright(MEAN_CURVATURE);	
}


void CMeshSimpDoc::OnUprightMeshsaliency()
{
	Upright(SALIENCY);
}


void CMeshSimpDoc::OnUprightMiniarea()
{
	Upright(MAX_AREA);
}


void CMeshSimpDoc::OnUprightViewshannon()
{
	Upright(VIEW_ENTROPY);
}

//////////////////////////////////////////////////////////////////////////
// Wang
void CMeshSimpDoc::OnUprightBasecur()
{
	Upright(BASE_CURVATURE);
}

void CMeshSimpDoc::OnUprightAll()
{
	// ↓使用的各种视点评分方法进行向上方向确定

	OnUprightMiniarea();
	OnUprightMeancur();
	OnUprightViewshannon();
	//OnUprightCurshannon();
	OnUprightMeshsaliency();
	// 使用Wang的方法 [3/24/2014 Han]
	//OnUprightBasecur();
	//OnUprightMcs();
}
// ↑使用的各种视点评分方法
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// ↓和观察球绘制,更新有关
void CMeshSimpDoc::UpdateViewsImportance(bool bFromView2Model)
{
	if (m_pViewPhereMdl == NULL)
		return;
	vector<ViewPoint*>::iterator iter;
	int n = 0;
	float sss = 0.0f;
	if (bFromView2Model)
	{
		m_pViewPhereMdl->MinImportanceV = FLT_MAX, m_pViewPhereMdl->MaxImportanceV = 0;
		for (iter = viewpoint_candidates.begin(); iter != viewpoint_candidates.end(); ++iter)
		{
			sss = (*iter)->importance;
			m_pViewPhereMdl->vertex(n).view_importance = sss;
			if (sss > m_pViewPhereMdl->MaxImportanceV)
				m_pViewPhereMdl->MaxImportanceV = sss;
			if (sss < m_pViewPhereMdl->MinImportanceV)
				m_pViewPhereMdl->MinImportanceV = sss;
			n++;
		}
		m_pViewPhereMdl->IdentityVertexImportance(/*m_pQslimMesh->MinImportance, m_pQslimMesh->MaxImportance*/);
		m_pViewPhereMdl->UpdateFaceImportanceByVert();
		m_pViewPhereMdl->VertexColor(3/*, m_pQslimMesh->MaxImportance, m_pQslimMesh->MinImportance*/);
	}
	else
	{
		for (iter = viewpoint_candidates.begin(); iter != viewpoint_candidates.end(); ++iter)
		{
			sss = m_pViewPhereMdl->vertex(n).view_importance;
			(*iter)->importance = sss;
			n++;
		}
		UpdateSphereColorByImportance();
	}
}
void CMeshSimpDoc::DrawArrow(vec3 posA, vec3 posB, float *color)
{
	if (color == NULL)
	{
		float f[3] = {1.0f, 0.f, 0.f};
		color = f;
	}
	vec3 vec = posB - posA;
	vec3 t0, t1, t2;
	//NxNormalToTangents(vec, t1, t2);
	vec3 c1 = vec.cross(vec3(0.0, 0.0, 1.0)); 
	vec3 c2 = vec.cross(vec3(0.0, 1.0, 0.0)); 
	if (c1.length() > c2.length())
		t1 = c1;	
	else
		t1 = c2;	
	t1.normalize();
	t2 = vec.cross(t1); 
	t2.normalize();

	t0 = posB - posA;
	t0.normalize();
	float ar = 0.05f;
	vec3 lobe1  = posB - t0*ar + t1 * ar;
	vec3 lobe2  = posB - t0*ar - t1 * ar;
	vec3 lobe3  = posB - t0*ar + t2 * ar;
	vec3 lobe4  = posB - t0*ar - t2 * ar;

	vec3 v3ArrowShape[] = {
		posA, posB,
		posB, lobe1,
		posB, lobe2,
		posB, lobe3,
		posB, lobe4,
	};

	glPushAttrib(GL_COLOR_BUFFER_BIT|GL_POINT_BIT|GL_POLYGON_BIT|GL_LINE_BIT|GL_DEPTH_BUFFER_BIT);
	//glDisable(GL_DEPTH_TEST);
	glEnable(GL_LINE_SMOOTH);
	glDisable(GL_LIGHTING);
	glLineWidth(3.0f);
	glColor3fv(color);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, sizeof(vec3), &v3ArrowShape[0].x);
	glDrawArrays(GL_LINES, 0, sizeof(v3ArrowShape)/sizeof(vec3));
	glDisableClientState(GL_VERTEX_ARRAY);
	glColor4f(1.0f,1.0f,1.0f,1.0f);

	glPopAttrib();
}

void CMeshSimpDoc::DrawUprightDir()
{
	//if (m_nWorstViewID == -1)
	//	return;
	if (m_baseViewpoint == NULL)
		return;
	float len = GetModelRadius();
	vec3 eyeP = 0.9f * len * m_baseViewpoint->pos;

	DrawArrow(eyeP, (-eyeP));

	// 再绘制最优视点方向
	//float c[3] = {0.f, 0.f, 1.f};
	//DrawArrow(vec3(viewpoint_candidates[m_nBestViewID]->pos[0], viewpoint_candidates[m_nBestViewID]->pos[1],viewpoint_candidates[m_nBestViewID]->pos[2]), vec3(0,0,0), c);
}
// ↑和观察球绘制,更新有关
//////////////////////////////////////////////////////////////////////////
