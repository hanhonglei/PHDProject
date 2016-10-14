#include "MSDM2.h"
#include<set>
#include<map>
#include<list>
#include<vector>
#include<stack>

double MSDM2::mini_radius = 0.002;
double MSDM2::radius_step = 0.001;
std::vector<vec3> MSDM2::TabPoint1;
std::vector<vec3> MSDM2::TabPoint2 ;

/*
MSDM2 algrithm steps:
1. A scale-dependent curvature is computed on vertices
from both meshes .
2. Each vertex of the distorted mesh Md is matched with
its corresponding 3D point and curvature value from the
reference mesh Mr, using fast projection and barycentric
interpolation.
3. For each vertex of Md a local distortion measure is
computed as a difference of Gaussian-weighted statistics
computed over a local spherical neighborhood of radius
hi.
4. Steps (1-3) are repeated for multiple scales hi, leading to
several distortion maps.
5. The final distortion map is constructed by adding the local
distortion maps at all scales (see section 4.4).
6. The global multiscale distortion score is then obtained by
combining the local values using a Minkowski pooling
*/


double MSDM2::MSDM2_computation(MxStdModel *pSrcMesh, MxStdModel *pSimpMesh, bool bSym)
{
	double MSDM2Value = 0.0;
		double MSDM2Value0;
		double MSDM2Value1;
	//Curvature_ComponentPtr component_ptr_curvature = findOrCreateComponentForViewer<Curvature_ComponentPtr, Curvature_Component>(viewer, polyhedron_ptr);

	if(!bSym)	// 非对称
	{
		//////////////processing here/////////////////////////////////////////////
		double maxdim=getMaxDim(pSrcMesh);
		ProcessMSDM2_Multires(pSrcMesh,pSimpMesh,MSDM_SCALE,maxdim,MSDM2Value);
	}
	else //symmetric
	{
		//////////////processing here/////////////////////////////////////////////
		// 得到包围盒最大边长
		double maxdim0=getMaxDim(pSrcMesh);
		double maxdim1=getMaxDim(pSimpMesh);


		ProcessMSDM2_Multires(pSimpMesh,pSrcMesh,MSDM_SCALE,maxdim1,MSDM2Value1);

		ProcessMSDM2_Multires(pSrcMesh,pSimpMesh,MSDM_SCALE,maxdim0,MSDM2Value0);

		MSDM2Value=(MSDM2Value0+MSDM2Value1)/2;

	}
	char pszError[1024];
	sprintf(pszError, "msdm2 error: %lf(sys);%lf(s2o);%lf(o2s)\n", MSDM2Value, MSDM2Value0, MSDM2Value1);
	MessageBox(NULL, LPCSTR(pszError), NULL, MB_ICONEXCLAMATION);
	//return;


	pSimpMesh->VertexColor(1);
	return MSDM2Value;
}
//
//void MSDM2::DistanceToColorMap()
//{
//	// active viewer
//	if (mw->activeMdiChild() != 0)
//	{
//#if defined (_MSC_VER) || ((defined (__linux__) || defined (__APPLE__)) && (CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(3,8,0)))
//		QApplication::setOverrideCursor(Qt::WaitCursor);
//
//		Viewer* viewer = (Viewer *)mw->activeMdiChild();
//		PolyhedronPtr polyhedron_ptr = viewer->getScenePtr()->get_polyhedron();
//
//		MSDM2_ComponentPtr component_ptr = findOrCreateComponentForViewer<MSDM2_ComponentPtr, MSDM2>(viewer, polyhedron_ptr);
//		component_ptr->ComputeMaxMin(polyhedron_ptr,1);	
//		component_ptr->ConstructColorMap(polyhedron_ptr,1);
//
//		viewer->recreateListsAndUpdateGL();
//
//		QApplication::restoreOverrideCursor();
//#else
//		QMessageBox::information(mw, APPLICATION, tr("Warning: to use this component under LINUX or MAC, you must install CGAL 3.8."));
//#endif
//	}	
//}


// 找到pSimMesh中每个顶点最近的pSrcMesh中的顶点
void MSDM2::Matching_Multires_Init( MxStdModel *pDesMesh, MxStdModel *pSrcMesh/*, MxFace * _TabMatchedFacet*/)
{

	// constructs AABB tree
	// 不进行aabb查找最近点，利用边折叠的简化序列，只寻找对应折叠点附近的最近点和最近平面 [4/19/2012 Han]
	//AABB_Tree tree(m_PolyOriginal->facets_begin(),m_PolyOriginal->facets_end());
	//tree.accelerate_distance_queries();

	//Searching for the closest point and facet for each vertex

	//int ind=0;
	//for(Vertex_iterator	pVertex	= m_PolyDegrad->vertices_begin();
	//	pVertex	!= m_PolyDegrad->vertices_end();
	//	pVertex++)
	//{
	//	pVertex->MSDM2_Local=0;
	//	// computes closest point and primitive id
	//	Point_and_primitive_id pp = tree.closest_point_and_primitive(pVertex->point());
	//	Point3d Nearest=pp.first;
	//	Facet_iterator f_Nearest = pp.second; // closest primitive id

	//	pVertex->match=Nearest;
	//	_TabMatchedFacet[ind]=*f_Nearest;
	//	ind++;		

	//}

	for (MxVertexID vID = 0; vID < pDesMesh->vert_count(); vID++)
	{
		if (!pDesMesh->vertex_is_valid(vID))
			continue;

		MxVertex &v = pDesMesh->vertex(vID);
		v.MSDM2_Local[0]=0;v.MSDM2_Local[1]=0;v.MSDM2_Local[2]=0;
		v.MSDM2_Local_sum = 0;
		FindNearest(pDesMesh, pSrcMesh, vID, vID);
	}
}

void MSDM2::Matching_Multires_Init( MxStdModel *pDesMesh, MxStdModel *pSrcMesh, MxVertexID fromVID, MxVertexID toVID)
{

	FindNearest(pDesMesh, pSrcMesh, toVID, pDesMesh->vertex(fromVID).nearestVertexID);
	if (pDesMesh->vertex(fromVID).nearestVertexID != pDesMesh->vertex(toVID).nearestVertexID)
		FindNearest(pDesMesh, pSrcMesh, toVID, pDesMesh->vertex(toVID).nearestVertexID);
}	

// 更新最近点
void MSDM2::FindNearest( MxStdModel *pDesMesh, MxStdModel *pSrcMesh, MxVertexID vID, MxVertexID nearestVID)
{
	if (!pDesMesh->vertex_is_valid(vID) || pSrcMesh == NULL || pDesMesh == NULL)
		return;

	MxVertex &v = pDesMesh->vertex(vID);
	vec3 vPos = vec3(v.as.pos);
	if(!pSrcMesh->vertex_is_valid(nearestVID))
		nearestVID = v.to;	
	const MxFaceList& nFL = pSrcMesh->neighbors(nearestVID);
	float disV = 1024;
	float disF = 1024;
	vec3 T[3];
	for(MxFaceID fSID = 0; fSID < nFL.length(); fSID++)
	{
		const MxFace &fS = pSrcMesh->face(nFL[fSID]);
		for(MxVertexID i = 0; i < 3; i++)
		{
			T[i] = vec3(pSrcMesh->vertex_position(fS.v[i]));
			float dV = (T[i] - vPos).length2();
			if (dV < disV)
			{
				disV = dV;
				v.nearestVertexID = fS.v[i];
			}				
		}
		float dF = vPos.Dist2ToTriangle(T);
		if (dF < disF)
		{
			disF = dF;
			v.nearestFaceID = nFL[fSID];
		}
	}
	// 保留匹配的最近点和最近平面信息
	v.match = vec3(pSrcMesh->vertex_position(v.nearestVertexID));
	//_TabMatchedFacet[vID] = pSrcMesh->face(v.nearestFaceID);
}

void MSDM2::Matching_Multires_Update( MxStdModel *pDesMesh, MxStdModel *pSrcMesh, int nScale)
{
	nScale == nScale < 0 ? 0:nScale;
	assert(nScale <= MSDM_SCALE);
	for(MxVertexID vID = 0; vID < pDesMesh->vert_count(); vID++)
	{
		if (!pDesMesh->vertex_is_valid(vID))
			continue;

		MxVertex &pVertex = pDesMesh->vertex(vID);
		pVertex.tag(vID);
		Matching_Multires_Update(pDesMesh, pSrcMesh, nScale, vID);
	}
}

void MSDM2::Matching_Multires_Update( MxStdModel *pDesMesh, MxStdModel *pSrcMesh, int nScale, int vID)
{
	if (pDesMesh == NULL || pSrcMesh == NULL)
		return;
	
	int ind=0;
	if (!pDesMesh->vertex_is_valid(vID))
		return;

	MxVertex &pVertex = pDesMesh->vertex(vID);

	MxFace &f_Nearest = pSrcMesh->face(pVertex.nearestFaceID);
	//Facet* f_Nearest=&_TabMatchedFacet[ind];

	//for debug
	//pVertex->point()=Nearest;

	///calculation of the nearest point curvature value using vertices of the Nearest triangle	
	//we use linear interpolation using barycentric coordinates
	//Point3d x1=f_Nearest->halfedge()->vertex()->point();
	//Point3d x2=f_Nearest->halfedge()->next()->vertex()->point();
	//Point3d x3=f_Nearest->halfedge()->next()->next()->vertex()->point();
	vec3 x1=vec3(pSrcMesh->vertex_position(f_Nearest.v[0]));
	vec3 x2=vec3(pSrcMesh->vertex_position(f_Nearest.v[1]));
	vec3 x3=vec3(pSrcMesh->vertex_position(f_Nearest.v[2]));

	double l1=sqrt((x3-x2).dot(x3-x2));
	double l2=sqrt((x1-x3).dot(x1-x3));
	double l3=sqrt((x1-x2).dot(x1-x2));

	//Vector v1=f_Nearest->halfedge()->vertex()->point()-pVertex->point();
	//Vector v2=f_Nearest->halfedge()->next()->vertex()->point()-pVertex->point();
	//Vector v3=f_Nearest->halfedge()->next()->next()->vertex()->point()-pVertex->point();
	vec3 v1=x1-vec3(pVertex.as.pos);
	vec3 v2=x2-vec3(pVertex.as.pos);
	vec3 v3=x3-vec3(pVertex.as.pos);

	double t1=sqrt(v1.dot(v1));
	double t2=sqrt(v2.dot(v2));
	double t3=sqrt(v3.dot(v3));

	double p1=(l1+t2+t3)/2;
	double p2=(t1+l2+t3)/2;
	double p3=(t1+t2+l3)/2;

	double A1=(p1*(p1-l1)*(p1-t3)*(p1-t2));
	double A2=(p2*(p2-l2)*(p2-t3)*(p2-t1));
	double A3=(p3*(p3-l3)*(p3-t1)*(p3-t2));

	if(A1>0) A1=sqrt(A1); else	A1=0;
	if(A2>0) A2=sqrt(A2); else	A2=0;
	if(A3>0) A3=sqrt(A3); else	A3=0;

	//double c1=f_Nearest->halfedge()->vertex()->KmaxCurv;
	//double c2=f_Nearest->halfedge()->next()->vertex()->KmaxCurv;
	//double c3=f_Nearest->halfedge()->next()->next()->vertex()->KmaxCurv;
	double c1=pSrcMesh->vertex(f_Nearest.v[0]).KmaxCurv[nScale];
	double c2=pSrcMesh->vertex(f_Nearest.v[1]).KmaxCurv[nScale];
	double c3=pSrcMesh->vertex(f_Nearest.v[2]).KmaxCurv[nScale];

	if((A1+A2+A3)>0)
		pVertex.curvmatch[nScale]=(A1*c1+A2*c2+A3*c3)/(A1+A2+A3);
	else
		pVertex.curvmatch[nScale]=(c1+c2+c3)/3;		

}

void MSDM2::ComputeStatistics(MxStdModel *pMesh, MxVertexID verID,  double Param,
	std::vector<double> & TabDistance1,std::vector<double>& TabDistance2,std::vector<vec3> &TabPoint1,std::vector<vec3> &TabPoint2,double radius, double dim, int nScale)
{
	MxVertex &pVertex = pMesh->vertex(verID);

	double moyenne1=0;
	double variance1=0;

	double moyenne2=0;
	double variance2=0;

	double covariance=0;

	////gaussian normalisation

	double variance = radius/2;
	double *tab_wi1=new double[TabPoint1.size()];
	double *tab_wi2=new double[TabPoint1.size()];
	double SommeDistance1=0;
	double SommeDistance2=0;
	double SommeWi1=0;
	double SommeWi2=0;
	for(unsigned int i=0;i<TabPoint1.size();i++) // MT
	{
		vec3 DistancePt1=TabPoint1[i]-vec3(pVertex.as.pos);
		double distPt1=sqrt(DistancePt1.dot(DistancePt1));
		double wi1=1/variance/sqrt(2*3.141592)*exp(-(distPt1*distPt1)/2/variance/variance);

		tab_wi1[i]=wi1;
		SommeWi1+=wi1;
		SommeDistance1+=TabDistance1[i]*wi1;

		vec3 DistancePt2=TabPoint2[i]-pVertex.match;
		double distPt2=sqrt(DistancePt2.dot(DistancePt2));

		double wi2=1/variance/sqrt(2*3.141592)*exp(-(distPt2*distPt2)/2/variance/variance);
		tab_wi2[i]=wi2;
		SommeWi2+=wi2;
		SommeDistance2+=TabDistance2[i]*wi2;
	}
	moyenne1=SommeDistance1/(double)SommeWi1;
	moyenne2=SommeDistance2/(double)SommeWi2;


	for(unsigned int i=0;i<TabPoint1.size();i++) // MT
	{
		variance1+=tab_wi1[i]*pow(TabDistance1[i]-moyenne1,2);
		variance2+=tab_wi2[i]*pow(TabDistance2[i]-moyenne2,2);
		covariance+=tab_wi2[i]*(TabDistance1[i]-moyenne1)*(TabDistance2[i]-moyenne2);
	}

	variance1=variance1/SommeWi1;
	variance2=variance2/SommeWi2;
	variance1=sqrt(variance1);
	variance2=sqrt(variance2);

	covariance=covariance/SommeWi1;

	///we then compute the MSDM2_Local value
	//double C1=1; // MT
	//double C2=1; // MT
	//double C3=0.5; // MT
#ifdef _MSC_VER
	double fact1=(fabs(moyenne1-moyenne2))/(max(moyenne1,moyenne2)+1);
	double fact2=(fabs(variance1-variance2))/(max(variance1,variance2)+1);

#else
	double fact1=(fabs(moyenne1-moyenne2))/(std::max(moyenne1,moyenne2)+1);
	double fact2=(fabs(variance1-variance2))/(std::max(variance1,variance2)+1);
#endif
	double fact3=(fabs(variance1*variance2-covariance))/(variance1*variance2+1);

	pVertex.MSDM2_Local[nScale]+=pow((pow(fact1,1)+pow(fact2,1)+Param*pow(fact3,1))/(2.+Param),1./1.);

	delete []tab_wi1;
	delete []tab_wi2;
}

double MSDM2::ProcessMSDM2_Multires(MxStdModel *pSrcMesh, MxStdModel *pDesMesh,int NbLevel, double maxdim,double & FastMSDM)
{	 
	double somme_MSDM2=0;
	int NbVertex=0;	

	//Facet * TabMatchedFacet=new Facet[pDesMesh->vert_count()];

	Matching_Multires_Init(pDesMesh, pSrcMesh/*,TabMatchedFacet*/);

	int nV = 0;

	double RadiusCurvature=mini_radius;
	for (int i=0;i<NbLevel;i++)
	{
		//m_PolyDegrad->set_index_vertices();
		pDesMesh->set_index_vertices(); 

		//component_ptr_curvature->principal_curvature(m_PolyOriginal,true,RadiusCurvature*maxdim);
		//component_ptr_curvature->principal_curvature(m_PolyDegrad,true,RadiusCurvature*maxdim);
		// 计算不同半径下的模型曲率
		pSrcMesh->principal_curvature(true,RadiusCurvature*maxdim, i);
		pDesMesh->principal_curvature(true,RadiusCurvature*maxdim, i);


		// 调整每个顶点的曲率
		KmaxKmean(pSrcMesh,maxdim, i);
		KmaxKmean(pDesMesh,maxdim, i);

		Matching_Multires_Update(pDesMesh, pSrcMesh, i);

		for(nV = 0;	nV < pDesMesh->vert_count(); nV++
			/*Vertex_iterator	pVertex	=	m_PolyDegrad->vertices_begin();
			pVertex	!= m_PolyDegrad->vertices_end();
			pVertex++*/)
		{
			if (!pDesMesh->vertex_is_valid(nV))
				continue;

			std::vector<double>  TabDistance1;std::vector<double> TabDistance2;
			TabPoint1.clear();
			TabPoint2.clear();

			// 计算每个顶点的MSDM2信息
			//ProcessMSDM2_per_vertex(pVertex,RadiusCurvature*5*maxdim,TabDistance1,TabDistance2,TabPoint1,TabPoint2);
			//ComputeStatistics((&(*pVertex)), 0.5,TabDistance1,TabDistance2,TabPoint1,TabPoint2,RadiusCurvature*5*maxdim,maxdim);	
			ProcessMSDM2_per_vertex(pDesMesh, nV,RadiusCurvature*5*maxdim,TabDistance1,TabDistance2,TabPoint1,TabPoint2, i);
			ComputeStatistics(pDesMesh, nV, 0.5,TabDistance1,TabDistance2,TabPoint1,TabPoint2,RadiusCurvature*5*maxdim,maxdim, i);	
			// han
			pDesMesh->vertex(nV).MSDM2_Local_sum += pDesMesh->vertex(nV).MSDM2_Local[i];
		}
		RadiusCurvature+=radius_step;
	}

	for(nV = 0;	nV < pDesMesh->vert_count(); nV++
		/*Vertex_iterator	pVertex	=	m_PolyDegrad->vertices_begin();
		pVertex	!= m_PolyDegrad->vertices_end();
		pVertex++*/)
	{
		if (!pDesMesh->vertex_is_valid(nV))
			continue;

		MxVertex &pVertex = pDesMesh->vertex(nV);

		//somme_MSDM2+=pow(pVertex.MSDM2_Local/NbLevel,3);	
		somme_MSDM2+=pow(pVertex.MSDM2_Local_sum/NbLevel,3);	
		NbVertex++;
	}


	FastMSDM=pow(somme_MSDM2/(double)NbVertex,0.33333);

	//delete [] TabMatchedFacet;
	//pDesMesh->IsDistanceComputed=true;
	return 0;
}

void MSDM2::KmaxKmean(MxStdModel *pMesh, double coef, int nScale)
{
	MxVertexID nV;
	for(nV = 0;	nV < pMesh->vert_count(); nV++)
	{
		if (!pMesh->vertex_is_valid(nV))
			continue;

		KmaxKmean(pMesh, coef, nScale, nV);
	}
}

void MSDM2::KmaxKmean(MxStdModel *pMesh, double coef, int nScale, MxVertexID nV)
{
	nScale == nScale < 0 ? 0:nScale;
	assert(nScale <= MSDM_SCALE);
	if (!pMesh->vertex_is_valid(nV))
		return;

	MxVertex &pVertex =  pMesh->vertex(nV);
	double kmax=pVertex.KmaxCurv[nScale]*coef;
	double kmin=pVertex.KminCurv[nScale]*coef;

	pVertex.KmaxCurv[nScale]=(kmax+kmin)/2.;
}
//
//void MSDM2::ConstructColorMap(PolyhedronPtr P, int MetricOrHausdorff)
//{
//	if(MetricOrHausdorff==1)
//	{
//		if(P->IsDistanceComputed==true)
//		{
//
//			double R;
//			int indiceLut;
//			Vertex_iterator pVertex = NULL;
//
//			if(MaxMSDM2>MinMSDM2)
//			{
//				for (pVertex = P->vertices_begin();
//					pVertex != P->vertices_end();
//					pVertex++)
//				{
//
//
//					R=(pVertex->MSDM2_Local-MinMSDM2)/(MaxMSDM2-MinMSDM2)*255;
//
//					indiceLut=floor(R);
//					pVertex->color(LUT_CourbureClust[3*indiceLut],LUT_CourbureClust[3*indiceLut+1],LUT_CourbureClust[3*indiceLut+2]);
//
//				}
//			}
//		}
//	}
//	else
//	{
//		if(P->IsDistanceComputed==true)
//		{
//
//			double R;
//			int indiceLut;
//			Vertex_iterator pVertex = NULL;
//
//			if(MaxMSDM2>MinMSDM2)
//			{
//				for (pVertex = P->vertices_begin();
//					pVertex != P->vertices_end();
//					pVertex++)
//				{
//
//					float d=sqrt((pVertex->point()-pVertex->match)*(pVertex->point()-pVertex->match));
//					R=(d-MinMSDM2)/(MaxMSDM2-MinMSDM2)*255;
//
//					indiceLut=floor(R);
//					pVertex->color(LUT_CourbureClust[3*indiceLut],LUT_CourbureClust[3*indiceLut+1],LUT_CourbureClust[3*indiceLut+2]);
//
//				}
//			}
//		}
//
//	}
//}



	//void MSDM2::ComputeMaxMin(PolyhedronPtr P, int MetricOrHausdorff)
	//{
	//	if(MetricOrHausdorff==1)
	//	{
	//		if(P->IsDistanceComputed==true)
	//		{
	//			MinMSDM2=1000;
	//			MaxMSDM2=0;
	//			for(Vertex_iterator	pVertex	=	P->vertices_begin();pVertex!= P->vertices_end();pVertex++)
	//			{
	//				if(pVertex->MSDM2_Local>MaxMSDM2)
	//					MaxMSDM2=pVertex->MSDM2_Local;
	//				if(pVertex->MSDM2_Local<MinMSDM2)
	//					MinMSDM2=pVertex->MSDM2_Local;

	//			}
	//		}
	//	}
	//	else
	//	{
	//		if(P->IsDistanceComputed==true)
	//		{
	//			MinMSDM2=1000;
	//			MaxMSDM2=0;
	//			for(Vertex_iterator	pVertex	=	P->vertices_begin();pVertex!= P->vertices_end();pVertex++)
	//			{
	//				float d=sqrt((pVertex->point()-pVertex->match)*(pVertex->point()-pVertex->match));
	//				if(d>MaxMSDM2)
	//					MaxMSDM2=d;
	//				if(d<MinMSDM2)
	//					MinMSDM2=d;

	//			}
	//		}

	//	}
	//

	//}

double MSDM2::getMaxDim(MxStdModel *pMesh)
{
	float size[3];
	pMesh->GetBoundSize(size);
	float max = 0;
	for (int i = 0; i < 3; i ++)
		max = max < size[i]? size[i] : max;

	return max;

}


bool sphere_clip_vector_MSDM2(vec3 &O, double r,const vec3 &P, vec3 &V)
{

	vec3 W = P - O ;
	double a = (V.dot(V));
	double b = 2.0 * V.dot(W) ;
	double c = (W.dot(W)) - r*r ;
	double delta = b*b - 4*a*c ;



	if( a==0)
		return true ;

	if(delta < 0) {
		// Should not happen, but happens sometimes (numerical precision)

		return true ;
	}
	double t = (- b + std::sqrt(delta)) / (2.0 * a) ;
	if(t < 0.0) {

		// Should not happen, but happens sometimes (numerical precision)
		return true ;
	}
	if(t >= 1.0) {
		// Inside the sphere
		return false ;
	}

	if(t==0)
	{

		t=0.01;
	}

	V=V*t;

	return true ;
}

	
void MSDM2::ProcessMSDM2_per_vertex(MxStdModel *pMesh, MxVertexID verID, double radius,std::vector<double> & TabDistance1
		,std::vector<double>& TabDistance2,std::vector<vec3> &TabPoint1,std::vector<vec3> &TabPoint2, int nScale)
{
	std::set<int> vertices ;

	std::stack<MxVertexID> S ;
	MxVertex &pVertex = pMesh->vertex(verID);

	vec3 O = vec3(pMesh->vertex_position(verID));
	//Point3d O = pVertex->point() ;

	//S.push(pVertex) ;	
	S.push(verID);
	vertices.insert(pVertex.tag()) ;


	// 当前模型的curvature和顶点信息
	TabDistance1.push_back(pVertex.KmaxCurv[nScale]);
	TabPoint1.push_back(/*pVertex->point()*/O);

	// 对应顶点的curvature和顶点信息
	TabDistance2.push_back(pVertex.curvmatch[nScale]);
	TabPoint2.push_back(pVertex.match);

	int NbSommetInSphere=0;
	//double SommeDistance=0; // MT

	while(!S.empty())
	{
		//Vertex_iterator v = S.top() ;
		MxVertexID v = S.top();
		// 当前顶点 [4/18/2012 Han]
		MxVertex &vCurr = pMesh->vertex(v);

		S.pop() ;
		//Point3d P = v->point() ;
		vec3 P = vec3(vCurr.as.pos);

		MxVertexList star;
		star.reset();
		pMesh->collect_vertex_star(v, star);
		//Halfedge_around_vertex_circulator h = v->vertex_begin();
		//Halfedge_around_vertex_circulator pHalfedgeStart = h;
		//CGAL_For_all(h,pHalfedgeStart)
		//MxEdgeList h;
		//pMesh->edge_around_vertex_circular(v, h);
		//for(uint j=0; j<h.length(); j++)
		for(int j = 0; j < star.length(); j++)
		{
			//const MxVertex &ve1 = pMesh->vertex(h[j].v1);
			//const MxVertex &ve2 = pMesh->vertex(h[j].v2);
			const MxVertex &ve1 = pMesh->vertex(v);
			const MxVertex &ve2 = pMesh->vertex(star(j));

			//Point3d p1 = h->vertex()->point();
			//Point3d p2 = h->opposite()->vertex()->point();
			vec3 p1 = vec3(ve1.as.pos);
			vec3 p2 = vec3(ve2.as.pos);;

			//Point3d p1m = h->vertex()->match;
			//Point3d p2m = h->opposite()->vertex()->match;
			vec3 p1m = ve1.match;
			vec3 p2m = ve2.match;

			vec3 V = (p2-p1);
			vec3 Vm = (p2m-p1m);

			if(v==verID || V.dot (P - O) > 0.0) 
			{
				double len_old = std::sqrt(V.dot(V));
				bool isect = sphere_clip_vector_MSDM2(O, radius, P, V) ;
				double len_edge = std::sqrt(V.dot(V));

				NbSommetInSphere++;

				double WeightedCurv1,WeightedCurv2;
				vec3 WeightedP1,WeightedP2;

				bool IsAlreadyIntegrated=false;
				if(!isect) 
				{
					//Vertex_iterator w=h->opposite()->vertex();
					//MxVertex &w=pMesh->vertex(h[j].v2);
					MxVertex &w=pMesh->vertex(star(j));

					if(vertices.find(w.tag()) == vertices.end())
					{
						vertices.insert(w.tag()) ;
						//S.push(h[j].v2) ;
						S.push(star(j));
					}
					else
						IsAlreadyIntegrated=true;
				}

				if (IsAlreadyIntegrated==false)
				{
					if(len_old!=0)
					{
						if(isect)
						{
							WeightedCurv1=(1-len_edge/len_old)*/*h->vertex()*/vCurr.KmaxCurv[nScale]+len_edge/len_old*/*h->opposite()->vertex()*/ve2.KmaxCurv[nScale];
							WeightedP1=p1+V;

							WeightedCurv2=(1-len_edge/len_old)*/*h->vertex()*/vCurr.curvmatch[nScale]+len_edge/len_old*/*h->opposite()->vertex()*/ve2.curvmatch[nScale];
							WeightedP2=p1m+(len_edge/len_old)*Vm;
						}
						else
						{
							WeightedCurv1=/*h->opposite()->vertex()*/ve2.KmaxCurv[nScale];
							WeightedCurv2=/*h->opposite()->vertex()*/ve2.curvmatch[nScale];
							WeightedP1=p2;
							WeightedP2=p2m;
						}
					}
					else
					{
						WeightedCurv1=/*h->opposite()->vertex()*/ve2.KmaxCurv[nScale];
						WeightedCurv2=/*h->opposite()->vertex()*/ve2.curvmatch[nScale];
						WeightedP1=p2;
						WeightedP2=p2m;
					}

					TabDistance1.push_back(WeightedCurv1);
					TabPoint1.push_back(WeightedP1);

					TabDistance2.push_back(WeightedCurv2);
					TabPoint2.push_back(WeightedP2);
				}
			}

		}

	}
}

MSDM2::MSDM2() 
{
}



///////
// han

void MSDM2::MSDM2SimpInit(MxStdModel *pSrcMesh, MxStdModel *pSimpMesh)
{
	MatchingSimpInit(pSimpMesh);
	if(pSrcMesh != NULL)
			MatchingSimpInit(pSrcMesh);
}

// 找到pSimMesh中每个顶点最近的pSrcMesh中的顶点
void MSDM2::MatchingSimpInit( MxStdModel *pMesh)
{
	for (MxVertexID vID = 0; vID < pMesh->vert_count(); vID++)
	{
		if (!pMesh->vertex_is_valid(vID))
			continue;

		MxVertex &v = pMesh->vertex(vID);
		const MxFaceList& nFL = pMesh->neighbors(vID);
		// init
		v.tag(vID);
		v.from = vID;
		v.to = vID;
		v.nearestFaceID = nFL(0);
		v.nearestVertexID = vID;
		v.MSDM2_Local_sum = 0;
		v.match = vec3(v.as.pos);
	}	

	double maxdim=getMaxDim(pMesh);
	double RadiusCurvature=mini_radius;	// 需要调整
	for (int i=0;i<MSDM_SCALE;i++)
	{
		// 计算不同半径下的模型曲率
		pMesh->principal_curvature(true,RadiusCurvature*maxdim, i);

		// 调整每个顶点的曲率
		KmaxKmean(pMesh,maxdim, i);

		//Matching_Multires_Update(pDesMesh, pSrcMesh, i);

		for(int nV = 0;	nV < pMesh->vert_count(); nV++)
		{
			if (!pMesh->vertex_is_valid(nV))
				continue;

			MxVertex &v = pMesh->vertex(nV);
			v.MSDM2_Local[i]=0;
			v.curvmatch[i] = v.KmaxCurv[i];

		}
		RadiusCurvature+=radius_step;
	}
}

// 更新当前顶点的所有信息，并计算由此引起的msdm误差
// 注意：由于两个点折叠为一个点，除了需要更新这个点本身的信息之外，这个点的相邻点以及相邻点周围半径范围内的所有点的信息都需要更新
void MSDM2::UpdateMSDM2(MxStdModel *pSrcMesh, MxStdModel *pSimpMesh, MxVertexID fromVID, MxVertexID toVID, double & FastMSDM ,bool bGlobalMsdm, bool bSym)
{
	// 首先更新from，to指针 [4/26/2012 Han]
	pSimpMesh->vertex(fromVID).to = toVID;
	pSimpMesh->vertex(toVID).from = fromVID;

	double somme_MSDM2=0;
	int NbVertex=0;	
	double maxdim=getMaxDim(pSimpMesh);
	//Facet * TabMatchedFacet=new Facet[pDesMesh->vert_count()];

	vec3 O = vec3(pSimpMesh->vertex_position(toVID));
	Matching_Multires_Init(pSimpMesh, pSrcMesh, fromVID, toVID);

	double RadiusCurvature=mini_radius;
	for (int i=0;i<MSDM_SCALE;i++)
	{
		// 遍历所有离toVID距离在radius内的顶点，并更新他们的信息
		std::set<MxVertexID> verts ;
		std::stack<MxVertexID> S ;
		S.push(toVID) ;
		verts.insert(toVID) ;
		while(!S.empty())
		{
			MxVertexID vID = S.top();
			MxVertex &vCurr = pSimpMesh->vertex(vID);
			S.pop() ;

			vec3 P = vec3(vCurr.as.pos);

			MxVertexList star;
			star.reset();
			pSimpMesh->collect_vertex_star(vID, star);
			for(int j = 0; j < star.length(); j++)
			{
				vec3 p1 = vec3(vCurr.as.pos);
				vec3 p2 = vec3(pSimpMesh->vertex_position(star(j)));;

				vec3 V = (p2-p1);
				//if (vID==toVID || V.dot(P - O) > 0.0)
				{
					bool isect = pSimpMesh->sphere_clip_vector(O, RadiusCurvature/*han 增加为直径*/*2, P, V);
					//bool isect = (p2 - O).length() <= RadiusCurvature;

					if (pSrcMesh != NULL)
						pSrcMesh->principal_curvature(true,RadiusCurvature*maxdim, i, vID);
					pSimpMesh->principal_curvature(true,RadiusCurvature*maxdim, i, vID);


					// 调整每个顶点的曲率
					if (pSrcMesh != NULL)
						KmaxKmean(pSrcMesh,maxdim, i, vID);
					KmaxKmean(pSimpMesh,maxdim, i, vID);

					Matching_Multires_Update(pSimpMesh, pSrcMesh, i, vID);

					if (!isect) {
						MxVertexID w = star(j);
						if (verts.find(w) == verts.end()) {
							verts.insert(w) ;
							S.push(w) ;
						}
					}
				}
			}
		}
		set<MxVertexID>::iterator it;

		for ( it=verts.begin() ; it != verts.end(); it++ )
		{
			MxVertexID nV = *it;
			if (!pSimpMesh->vertex_is_valid(nV))
				continue;

			std::vector<double>  TabDistance1;std::vector<double> TabDistance2;
			TabPoint1.clear();
			TabPoint2.clear();

			// 计算每个顶点的MSDM2信息
			//ProcessMSDM2_per_vertex(pVertex,RadiusCurvature*5*maxdim,TabDistance1,TabDistance2,TabPoint1,TabPoint2);
			//ComputeStatistics((&(*pVertex)), 0.5,TabDistance1,TabDistance2,TabPoint1,TabPoint2,RadiusCurvature*5*maxdim,maxdim);	
			ProcessMSDM2_per_vertex(pSimpMesh, nV,RadiusCurvature*5*maxdim,TabDistance1,TabDistance2,TabPoint1,TabPoint2, i);
			ComputeStatistics(pSimpMesh, nV, 0.5,TabDistance1,TabDistance2,TabPoint1,TabPoint2,RadiusCurvature*5*maxdim,maxdim, i);	
			// han
			pSimpMesh->vertex(nV).MSDM2_Local_sum += pSimpMesh->vertex(nV).MSDM2_Local[i];
		}

		RadiusCurvature+=radius_step;
	}

	if (!bGlobalMsdm)		// 无需得到总体msdm信息
		return;
	double minM, maxM;
	GetMSDMInfo(pSimpMesh, FastMSDM, minM, maxM);
}

void MSDM2::GetMSDMInfo(MxStdModel* pMesh, double & FastMSDM ,double &MinLocalMSDM, double &MaxLocalMSDM)
{	
	MinLocalMSDM=1000;
	MaxLocalMSDM=0;
	double somme_MSDM2 = 0;
	int NbVertex = 0;

	for(MxVertexID nV = 0;	nV < pMesh->vert_count(); nV++
		/*Vertex_iterator	pVertex	=	m_PolyDegrad->vertices_begin();
		pVertex	!= m_PolyDegrad->vertices_end();
		pVertex++*/)
	{
		if (!pMesh->vertex_is_valid(nV))
			continue;

		MxVertex &pVertex = pMesh->vertex(nV);

		if(pVertex.MSDM2_Local_sum>MaxLocalMSDM)
			MaxLocalMSDM=pVertex.MSDM2_Local_sum;
		if(pVertex.MSDM2_Local_sum<MinLocalMSDM)
			MinLocalMSDM=pVertex.MSDM2_Local_sum;

		//somme_MSDM2+=pow(pVertex.MSDM2_Local/NbLevel,3);	
		somme_MSDM2+=pow(pVertex.MSDM2_Local_sum/MSDM_SCALE,3);	
		NbVertex++;
	}
	
	FastMSDM=pow(somme_MSDM2/(double)NbVertex,0.33333);
}