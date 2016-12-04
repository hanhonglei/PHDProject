/************************************************************************

  Surface simplification using quadric error metrics

  Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.
  
  $Id: MxQSlim.cxx,v 1.42.2.2 2004/07/01 18:47:32 garland Exp $

 ************************************************************************/

#include "stdmix.h"
#include "MxQSlim.h"
#include "MxGeom3D.h"
#include "MxVector.h"

typedef MxQuadric3 Quadric;

MxQSlim::MxQSlim(MxStdModel& _m)
    : MxStdSlim(&_m),
      quadrics(_m.vert_count())
{
    // Externally visible variables
    object_transform = NULL;

	m_bUpdateMSDM = false;
}

void MxQSlim::initialize()
{
    collect_quadrics();
    if( boundary_weight > 0.0 )
	constrain_boundaries();
    if( object_transform )
	transform_quadrics(*object_transform);

    is_initialized = true;
}
//
//
//float getpixelsize(float fovy, float planedist , float imageHeight)
//{
//	return tan(fovy/2.0f/180.0f*3.141592653)*planedist/imageHeight*2.0f;
//}
//
//// 计算每个顶点的模糊值 [11/25/2011 Han]
//bool CalcBlurWeights(/*Vec3 eye,Vec3 gaze,Vec3 up*/)
//{
//	//Mesh m;
//	Matrix curmx = translate(0,0,0);
//	Matrix mx1 = translate(1,1,1);
//	Matrix mx2 = rotate( Vec3(1,1,1) , 60);
//	double areasum = 0 , perareasum = 0;
//	Matrix projectmx , transfermx , permx;
//	float imgplanedist=0.5;
//	Vec3 eye = Vec3(0,0,2);
//	Vec3 gaze = Vec3(0,0,-1);
//	Vec3 up = Vec3(0,1,0);
//
//	//m.ReadObjFile("Wei_Bunny01(0.00081).obj");
//	//m.calTriangleArea();
//	//Vec3 maxscale = m.m_max - m.m_min;
//	//float f_scale = max(max( maxscale.x() , maxscale.y()) , maxscale.z());
//	//m.meshTransform(scale(1.f/f_scale , 1.f/f_scale , 1.f/f_scale));
//	//Vec3 diag = m.m_max - m.m_min;
//	//float f_scale = 1.0 / diag.length() * 1.5;
//	//Vec3 trans = -f_scale * (m.m_min + 0.5 * (m.m_max - m.m_min));
//	//m.meshTransform(translate(trans.x , trans.y , trans.z));
//	//m.meshTransform(scale(f_scale ,f_scale , f_scale));
//	//m.calTriangleArea();
//	//std::ofstream fin("d:\\Wei_Bunny01(0.00081).blr"); 
//
//	//
//	projectmx = viewMatrix( eye , gaze , up );
//	transfermx = translate(0.128,0,0);
//	permx = zeroMatrix();
//	permx.x[0][0] = 1; permx.x[1][1] = 1; permx.x[2][2] = 1; permx.x[3][2] = 1/imgplanedist;
//
//	Vec3 orignalpos , laterpos;
//	float pixelsize = getpixelsize(30 , 1 , 600);
//	float dd,sumdd=0;
//	for (int i = 0; i < _newmesh.getNumVerts(); i++)
//	{
//		orignalpos = transformLoc( permx*projectmx , _newmesh.getVertex(i).getXYZ());
//		laterpos = transformLoc( transfermx , _newmesh.getVertex(i).getXYZ());
//		dd = (orignalpos - transformLoc( permx*projectmx , laterpos )).length();
//		_newmesh.getVertex(i)._blur = dd;
//		_mesh->getVertex(i)._blur = dd;
//
//		//dd = dist( orignalpos , transformLoc( permx*projectmx , laterpos ) );
//		//printf("%f\n",dd);
//		sumdd += dd;
//	}
//	return true;
//}

/*
其他3篇利用QSlim中的quadric调节进行简化的权重设置方法，都是将Q矩阵乘以一个权重weight

** Perceptually Guided Polygon Reduction
weight = 1.0/(1.0+pow(I(v), α)
其中I(v)表示顶点的重要程度，α用来调节相对权重，文章使用2

**Mesh Saliency
			{ λl(v)	if l(v) >= α
weight =	| 
			{ l(v)		if l(v) < α
其中λ = 100， α表示前30个重要度,l(v)表示顶点重要度

**Visibility-Guided Simplification
weight = V(t)^2
V(t)表示三角形法线和视线的点乘加权平均的结果。每个三角形的误差是每个视线方向的加权平均投影距离。
这表示当某个边处于低可见区域的时候，越有可能被简化掉。文章使用了原始的quadric矩阵来得到新顶点位置。
*/
void MxQSlim::collect_quadrics()
{
    uint j;

    for(j=0; j<quadrics.length(); j++)
	quadrics(j).clear();

	for(MxFaceID i=0; i<m->face_count(); i++)
    {
	MxFace& f = m->face(i);

	Vec3 v1(m->vertex(f(0)));
	Vec3 v2(m->vertex(f(1)));
	Vec3 v3(m->vertex(f(2)));

	Vec4 p = (weighting_policy==MX_WEIGHT_RAWNORMALS) ?
		    triangle_raw_plane<Vec3,Vec4>(v1, v2, v3):
		    triangle_plane<Vec3,Vec4>(v1, v2, v3);
	Quadric Q(p[X], p[Y], p[Z], p[W], m->compute_face_area(i));

	switch( weighting_policy )
	{
	case MX_WEIGHT_ANGLE:
	    for(j=0; j<3; j++)
	    {
		Quadric Q_j = Q;
		Q_j *= m->compute_corner_angle(i, j);
		quadrics(f[j]) += Q_j;
	    break;
		}
	case MX_WEIGHT_AREA:
	case MX_WEIGHT_AREA_AVG:
	    Q *= Q.area();
	    // no break: fallthrough
	default:
	    quadrics(f[0]) += Q;
	    quadrics(f[1]) += Q;
	    quadrics(f[2]) += Q;
	    break;
	}
    }


	// 为每个顶点添加权重：从文件载入motion blr数据或者直接计算每个顶点的旋转半径 [3/30/2012 Han]
#define USE_VERTEX_BLUR
//#define USE_VERTEX_RADIUS


	// 1.利用顶点权重来调整简化顺序:在motionblur中,权重表示顶点的速度.速度越快,越容易简化 [3/30/2012 Han]
	// 2.利用sPower来调整权重对简化的影响:0表示权重不产生影响,1表示产生最大影响
	double fMin = FLT_MAX, fMax = 0;
	float fTmpMin = FLT_MAX, fTmpMax = 0;
	// test, add blur weights [11/28/2011 Han]
	float fMinGard = 0.1, fMaxGard = 10.0;


#ifdef USE_VERTEX_BLUR
	if (m->pVertexMBlur.size() == 0 || m->pVertexMBlur.size() < m->vert_count()){
		return;
	}
	// 计算顶点的最大最小motion weight [12/16/2011 Han]
	for (MxVertexID i = 0; i < m->vert_count(); i++)
	{
		float tmp = m->pVertexMBlur[i];
		if (tmp < fMin)
			fMin = tmp;
		if (tmp > fMax)
			fMax = tmp;
	}
#else
#ifdef USE_VERTEX_RADIUS
	// 下面直接采取的是顶点离模型中心的距离作为简化参数，影响简化结果
	// 即:离中心的距离越远,越容易简化
	m->pVertexMBlur.clear();
	Vec3 vMax, vMin, vCent;
	vMax[0] = -1024;
	vMax[1] = -1024;
	vMax[2] = -1024;
	vMin[0] = 1024;
	vMin[1] = 1024;
	vMin[2] = 1024;
	for (MxVertexID i = 0; i < m->vert_count(); i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if(m->vertex(i).as.pos[j] < vMin[j])
				vMin[j] = m->vertex(i).as.pos[j];
			if(m->vertex(i).as.pos[j] > vMax[j])
				vMax[j] = m->vertex(i).as.pos[j];
		}
	}
	// 计算模型的包围盒中心，按照顶点到中心的平面距离作为旋转时的速度参数 [12/16/2011 Han]
	vCent = (vMax + vMin) / 2.0;

	// 计算每个顶点到中心的平面距离并保存 [12/16/2011 Han]
	for (MxVertexID i = 0; i < m->vert_count(); i++)
	{
		Vec3 vAdjust;
		for (int j = 0; j < 3; j++)
			vAdjust[j] = m->vertex(i).as.pos[j] - vCent[j];

		//float tmp = /*1- 100**/(m->vertex(i).as.pos[0]*m->vertex(i).as.pos[0]+m->vertex(i).as.pos[2]*m->vertex(i).as.pos[2]);
		float tmp = pow(vAdjust[0], 2) + pow(vAdjust[2], 2);
		m->pVertexMBlur.push_back(tmp);
		if (tmp < fMin)
			fMin = tmp;
		if (tmp > fMax)
			fMax = tmp;
	}
#endif
#endif
	if (m->pVertexMBlur.size() == 0 || m->pVertexMBlur.size() < m->vert_count()){
		return;
	}
	// 将在最大和最小权重保存于最后两位 [3/30/2012 Han]
	m->pVertexMBlur.push_back(fMax);
	m->pVertexMBlur.push_back(fMin);


	for (MxVertexID i = 0; i < m->vert_count(); i++)
	{
		double tmp = m->pVertexMBlur[i];/*1- 100**//*(m->vertex(f(j))[2]*m->vertex(f(j))[2]+m->vertex(f(j))[0]*m->vertex(f(j))[0]);*/
		// 直接乘一个权重太过集中 [11/28/2011 Han]
		//quadrics(f[j]) *= abs(tmp);
		//if(abs(tmp)<abs(fMin)) 
		//	tmp = fMin;
		//tmp = fMinGard+ (fMaxGard-fMinGard)*(abs(tmp) - abs(fMin));
		/*tmp = fMaxGard - (fMinGard+ (fMaxGard-fMinGard)*(tmp - fMin)/(fMax-fMin));*/
		//tmp = 1.0/(1.0+1000*tmp/*pow((tmp - fMin)/(fMax-fMin), 2)*/);
		//if (m->vertex(i).as.pos[2] < 0)
		//	tmp = 0.0001;
		//else
		//	tmp = 1;
		//  [12/16/2011 Han]
		// 暂时删除，测试其他方法 [3/30/2012 Han]
		tmp = /*fMaxGard -*/ 1.0 / (1.0 + (fMinGard+ (fMaxGard-fMinGard)*(tmp - fMin)/(fMax-fMin)));

		// 计算weight影响Q的方法:
		// 1.将weight调整到0~1之间
		// 2.利用幂指数来调整期对Q的影响
		// test
		//tmp = 1.0 + 2.0*(1.0 / (1.0 + (tmp - fMin)/(fMax-fMin)) - 0.5);
		//if (tmp > 0.6)
		//	tmp *= 100;
		
		quadrics(i) *= pow(float(tmp), m->fWeightPower*7);
		if (tmp < fTmpMin)
			fTmpMin = tmp;
		if (tmp > fTmpMax)
			fTmpMax = tmp;
	}
}
/*
其他3篇利用QSlim中的quadric调节进行简化的权重设置方法，都是将Q矩阵乘以一个权重weight

** Perceptually Guided Polygon Reduction
weight = 1.0/(1.0+pow(I(v), α)
其中I(v)表示顶点的重要程度，α用来调节相对权重，文章使用2

**Mesh Saliency
			{ λl(v)	if l(v) >= α
weight =	| 
			{ l(v)		if l(v) < α
其中λ = 100， α表示前30个重要度,l(v)表示顶点重要度

**Visibility-Guided Simplification
weight = V(t)^2
V(t)表示三角形法线和视线的点乘加权平均的结果。每个三角形的误差是每个视线方向的加权平均投影距离。
这表示当某个边处于低可见区域的时候，越有可能被简化掉。文章使用了原始的quadric矩阵来得到新顶点位置。
*/

//////////////////////////////////////////////////////////////////////////
// 依靠距离因素调整LOD级别 [9/28/2011 Han Honglei]
// 屏幕投影为1个像素的时候，反投影到模型所在位置d后的长度Δd = 0.0254*d/(2*r*l),r是屏幕分辨率（单位是dpi），l是近裁剪面的距离
double MxQSlim::AdjustLOD(float d, float l, double fPixelInMeter)
{
	// calc 1 pixel unprojected to mesh's position, this is the largest error can be tolerated by viewer [3/27/2012 Han]
	//float QError = pow(0.0254*d/(2*r*l), 2);

	double fUnprojPixel = fPixelInMeter*d/l;
	double QError = pow(fUnprojPixel, 2);
	DecimateByError(QError);
	return fUnprojPixel;
}

void MxQSlim::transform_quadrics(const Mat4& P)
{
    for(uint j=0; j<quadrics.length(); j++)
	quadrics(j).transform(P);
}

void MxQSlim::discontinuity_constraint(MxVertexID i, MxVertexID j,
				       const MxFaceList& faces)
{
    for(uint f=0; f<faces.length(); f++)
    {
	Vec3 org(m->vertex(i)), dest(m->vertex(j));
	Vec3 e = dest - org;

	Vec3 n;
	m->compute_face_normal(faces[f], n);

	Vec3 n2 = e ^ n;
	unitize(n2);

	MxQuadric3 Q(n2, -(n2*org));
	Q *= boundary_weight;

	if( weighting_policy == MX_WEIGHT_AREA ||
	    weighting_policy == MX_WEIGHT_AREA_AVG )
	{
	    Q.set_area(norm2(e));
	    Q *= Q.area();
	}

	quadrics(i) += Q;
	quadrics(j) += Q;
    }
}

void MxQSlim::constrain_boundaries()
{
    MxVertexList star;
    MxFaceList faces;

    for(MxVertexID i=0; i<m->vert_count(); i++)
    {
	star.reset();
	m->collect_vertex_star(i, star);

	for(uint j=0; j<star.length(); j++)
	    if( i < star(j) )
	    {
		faces.reset();
		m->collect_edge_neighbors(i, star(j), faces);
		if( faces.length() == 1 )
		    discontinuity_constraint(i, star(j), faces);
	    }
    }
}




MxEdgeQSlim::MxEdgeQSlim(MxStdModel& _m)
  : MxQSlim(_m),
    edge_links(_m.vert_count())
{
    contraction_callback = NULL;
}

MxEdgeQSlim::~MxEdgeQSlim()
{
    // Delete everything remaining in the heap
    for(uint i=0; i<heap.size(); i++)
	delete ((MxQSlimEdge *)heap.item(i));
}

///////////////////////////////////////////////////////////////////////////
//
// IMPORTANT NOTE:  These check_xxx functions assume that the local
//                  neighborhoods have been marked so that each face
//                  is marked with the number of involved vertices it has.
//

double MxEdgeQSlim::check_local_compactness(uint v1, uint/*v2*/,
					    const float *vnew)
{
    const MxFaceList& N1 = m->neighbors(v1);
    double c_min = 1.0;

    for(uint i=0; i<N1.length(); i++)
	if( m->face_mark(N1[i]) == 1 )
	{
	    const MxFace& f = m->face(N1[i]);
	    Vec3 f_after[3];
	    for(uint j=0; j<3; j++)
		f_after[j] = (f[j]==v1)?Vec3(vnew):Vec3(m->vertex(f[j]));

	    double c=triangle_compactness(f_after[0], f_after[1], f_after[2]);

	    if( c < c_min ) c_min = c;
	}

    return c_min;
}

double MxEdgeQSlim::check_local_inversion(uint v1,uint/*v2*/,const float *vnew)
{
    double Nmin = 1.0;
    const MxFaceList& N1 = m->neighbors(v1);

    for(uint i=0; i<N1.length(); i++)
	if( m->face_mark(N1[i]) == 1 )
	{
	    const MxFace& f = m->face(N1[i]);
	    Vec3 n_before;
	    m->compute_face_normal(N1[i], n_before);

	    Vec3 f_after[3];
	    for(uint j=0; j<3; j++)
		f_after[j] = (f[j]==v1)?Vec3(vnew):Vec3(m->vertex(f[j]));

	    double delta = n_before *
		triangle_normal(f_after[0], f_after[1], f_after[2]);

	    if( delta < Nmin ) Nmin = delta;
	}

    return Nmin;
}

uint MxEdgeQSlim::check_local_validity(uint v1, uint /*v2*/, const float *vnew)

{
    const MxFaceList& N1 = m->neighbors(v1);
    uint nfailed = 0;
    uint i;

    for(i=0; i<N1.length(); i++)
	if( m->face_mark(N1[i]) == 1 )
	{
	    MxFace& f = m->face(N1[i]);
	    uint k = f.find_vertex(v1);
	    uint x = f[(k+1)%3];
	    uint y = f[(k+2)%3];

	    float d_yx[3], d_vx[3], d_vnew[3], f_n[3], n[3];
	    mxv_sub(d_yx, m->vertex(y),  m->vertex(x), 3);   // d_yx = y-x
	    mxv_sub(d_vx, m->vertex(v1), m->vertex(x), 3);   // d_vx = v-x
	    mxv_sub(d_vnew, vnew, m->vertex(x), 3);          // d_vnew = vnew-x

	    mxv_cross3(f_n, d_yx, d_vx);
	    mxv_cross3(n, f_n, d_yx);     // n = ((y-x)^(v-x))^(y-x)
	    mxv_unitize(n, 3);

	    // assert( mxv_dot(d_vx, n, 3) > -FEQ_EPS );
	    if(mxv_dot(d_vnew,n,3)<local_validity_threshold*mxv_dot(d_vx,n,3))
		nfailed++;
	}

    return nfailed;
}

uint MxEdgeQSlim::check_local_degree(uint v1, uint v2, const float *)
{
    const MxFaceList& N1 = m->neighbors(v1);
    const MxFaceList& N2 = m->neighbors(v2);
    uint i;
    uint degree = 0;

    // Compute the degree of the vertex after contraction
    //
    for(i=0; i<N1.length(); i++)
	if( m->face_mark(N1[i]) == 1 )
	    degree++;

    for(i=0; i<N2.length(); i++)
	if( m->face_mark(N2[i]) == 1 )
	    degree++;


    if( degree > vertex_degree_limit )
	return degree - vertex_degree_limit;
    else
	return 0;
}

void MxEdgeQSlim::apply_mesh_penalties(MxQSlimEdge *info)
{
    uint i;

    const MxFaceList& N1 = m->neighbors(info->v1);
    const MxFaceList& N2 = m->neighbors(info->v2);

    // Set up the face marks as the check_xxx() functions expect.
    //
    for(i=0; i<N2.length(); i++) m->face_mark(N2[i], 0);
    for(i=0; i<N1.length(); i++) m->face_mark(N1[i], 1);
    for(i=0; i<N2.length(); i++) m->face_mark(N2[i], m->face_mark(N2[i])+1);

    double base_error = info->heap_key();
    double bias = 0.0;
    
    // Check for excess over degree bounds.
    //
    uint max_degree = MAX(N1.length(), N2.length());
    if( max_degree > vertex_degree_limit )
	bias += (max_degree-vertex_degree_limit) * meshing_penalty * 0.001;

#if ALTERNATE_DEGREE_BIAS
    // ??BUG:  This code was supposed to be a slight improvement over
    //         the earlier version above.  But it performs worse.
    //         Should check into why sometime.
    //
    uint degree_excess = check_local_degree(info->v1, info->v2, info->vnew);
    if( degree_excess )
 	bias += degree_excess * meshing_penalty;
#endif

    // Local validity checks
    //
    uint nfailed = check_local_validity(info->v1, info->v2, info->vnew);
    nfailed += check_local_validity(info->v2, info->v1, info->vnew);
    if( nfailed )
	bias += nfailed*meshing_penalty;

    if( compactness_ratio > 0.0 )
    {
	double c1_min=check_local_compactness(info->v1, info->v2, info->vnew);
	double c2_min=check_local_compactness(info->v2, info->v1, info->vnew);
	double c_min = MIN(c1_min, c2_min);

	// !!BUG: There's a small problem with this: it ignores the scale
	//        of the errors when adding the bias.  For instance, enabling
	//        this on the cow produces bad results.  I also tried
	//        += (base_error + FEQ_EPS) * (2-c_min), but that works
	//        poorly on the flat planar thing.  A better solution is
	//        clearly needed.
	//
	//  NOTE: The prior heuristic was
	//        if( ratio*cmin_before > cmin_after ) apply penalty;
	//
	if( c_min < compactness_ratio )
	    bias += (1-c_min);
    }

#if USE_OLD_INVERSION_CHECK
    double Nmin1 = check_local_inversion(info->v1, info->v2, info->vnew);
    double Nmin2 = check_local_inversion(info->v2, info->v1, info->vnew);
    if( MIN(Nmin1, Nmin2) < 0.0 )
	bias += meshing_penalty;
#endif

    info->heap_key(base_error - bias);
}

// 此处寻找顶点V1和V2折叠以后的最佳位置，并且计算误差，将误差最为后续简化的依据 [3/30/2012 Han]
void MxEdgeQSlim::compute_target_placement(MxQSlimEdge *info)
{
    MxVertexID i=info->v1, j=info->v2;

    const Quadric &Qi=quadrics(i), &Qj=quadrics(j);

    Quadric Q = Qi;  Q += Qj;
    double e_min;

    if( placement_policy==MX_PLACE_OPTIMAL &&
	Q.optimize(&info->vnew[X], &info->vnew[Y], &info->vnew[Z]) )
    {
	e_min = Q(info->vnew);
    }
    else
    {
	Vec3 vi(m->vertex(i)), vj(m->vertex(j));	
	Vec3 best;

	if( placement_policy>=MX_PLACE_LINE && Q.optimize(best, vi, vj) )
	    e_min = Q(best);
	else
	{
	    double ei=Q(vi), ej=Q(vj);

	    if( ei < ej ) { e_min = ei; best = vi; }
	    else          { e_min = ej; best = vj; }

	    if( placement_policy>=MX_PLACE_ENDORMID )
	    {
		Vec3 mid = (vi+vj)/2.0;
		double e_mid = Q(mid);

		if( e_mid < e_min ) { e_min = e_mid; best = mid; }
	    }
	}

	info->vnew[X] = best[X];
	info->vnew[Y] = best[Y];
	info->vnew[Z] = best[Z];
    }

    if( weighting_policy == MX_WEIGHT_AREA_AVG )
 	e_min /= Q.area();

    info->heap_key(-e_min);
}

void MxEdgeQSlim::finalize_edge_update(MxQSlimEdge *info)
{
    if( meshing_penalty > 1.0 )
	apply_mesh_penalties(info);

    if( info->is_in_heap() )
	heap.update(info);
    else
	heap.insert(info);
}


void MxEdgeQSlim::compute_edge_info(MxQSlimEdge *info)
{
    compute_target_placement(info);

    finalize_edge_update(info);
}

void MxEdgeQSlim::create_edge(MxVertexID i, MxVertexID j)
{
    MxQSlimEdge *info = new MxQSlimEdge;

    edge_links(i).add(info);
    edge_links(j).add(info);

    info->v1 = i;
    info->v2 = j;

    compute_edge_info(info);
}

void MxEdgeQSlim::collect_edges()
{
    MxVertexList star;

    for(MxVertexID i=0; i<m->vert_count(); i++)
    {
	star.reset();
	m->collect_vertex_star(i, star);

	for(uint j=0; j<star.length(); j++)
	    if( i < star(j) )  // Only add particular edge once
		create_edge(i, star(j));
    }
}

void MxEdgeQSlim::initialize()
{
    MxQSlim::initialize();
    collect_edges();
}

void MxEdgeQSlim::initialize(const MxEdge *edges, uint count)
{
    MxQSlim::initialize();
    for(uint i=0; i<count; i++)
	create_edge(edges[i].v1, edges[i].v2);
}

void MxEdgeQSlim::update_pre_contract(const MxPairContraction& conx)
{
    MxVertexID v1=conx.v1, v2=conx.v2;
    uint i, j;

    star.reset();
    //
    // Before, I was gathering the vertex "star" using:
    //      m->collect_vertex_star(v1, star);
    // This doesn't work when we initially begin with a subset of
    // the total edges.  Instead, we need to collect the "star"
    // from the edge links maintained at v1.
    //
    for(i=0; i<edge_links(v1).length(); i++)
	star.add(edge_links(v1)[i]->opposite_vertex(v1));

    for(i=0; i<edge_links(v2).length(); i++)
    {
	MxQSlimEdge *e = edge_links(v2)(i);
	MxVertexID u = (e->v1==v2)?e->v2:e->v1;
	SanityCheck( e->v1==v2 || e->v2==v2 );
	SanityCheck( u!=v2 );

	if( u==v1 || varray_find(star, u) )
	{
	    // This is a useless link --- kill it
	    bool found = varray_find(edge_links(u), e, &j);
	    assert( found );
	    edge_links(u).remove(j);
	    heap.remove(e);
	    if( u!=v1 ) delete e; // (v1,v2) will be deleted later
	}
	else
	{
	    // Relink this to v1
	    e->v1 = v1;
	    e->v2 = u;
	    edge_links(v1).add(e);
	}
    }

    edge_links(v2).reset();
}

void MxEdgeQSlim::update_post_contract(const MxPairContraction& conx)
{
}

void MxEdgeQSlim::apply_contraction(const MxPairContraction& conx, bool bUpdateMSDM)
{
    //
    // Pre-contraction update
    valid_verts--;
    valid_faces -= conx.dead_faces.length();
    quadrics(conx.v1) += quadrics(conx.v2);

    update_pre_contract(conx);

    m->apply_contraction(conx);

    update_post_contract(conx);

    // Must update edge info here so that the meshing penalties
    // will be computed with respect to the new mesh rather than the old
    for(uint i=0; i<edge_links(conx.v1).length(); i++)
	compute_edge_info(edge_links(conx.v1)[i]);

	// update vertex's normal [5/28/2012 Han]
	float n[3];
	m->compute_vertex_normal(conx.v1, n);
	m->normal(conx.v1).set(n);


	// 完成边折叠以后更新msdm信息 [4/26/2012 Han]
	if (bUpdateMSDM)
		m->UpdateAfterSimp(conx.v1, conx.v2);	
}

void MxEdgeQSlim::update_pre_expand(const MxPairContraction&)
{
}

void MxEdgeQSlim::update_post_expand(const MxPairContraction& conx)
{
    MxVertexID v1=conx.v1, v2=conx.v2;
    uint i;

    star.reset(); star2.reset();
    PRECAUTION(edge_links(conx.v2).reset());
    m->collect_vertex_star(conx.v1, star);
    m->collect_vertex_star(conx.v2, star2);

    i = 0;
    while( i<edge_links(v1).length() )
    {
	MxQSlimEdge *e = edge_links(v1)(i);
	MxVertexID u = (e->v1==v1)?e->v2:e->v1;
	SanityCheck( e->v1==v1 || e->v2==v1 );
	SanityCheck( u!=v1 && u!=v2 );

	bool v1_linked = varray_find(star, u);
	bool v2_linked = varray_find(star2, u);

	if( v1_linked )
	{
	    if( v2_linked )  create_edge(v2, u);
	    i++;
	}
	else
	{
	    // !! BUG: I expected this to be true, but it's not.
	    //         Need to find out why, and whether it's my
	    //         expectation or the code that's wrong.
	    // SanityCheck(v2_linked);
	    e->v1 = v2;
	    e->v2 = u;
	    edge_links(v2).add(e);
	    edge_links(v1).remove(i);
	}

	compute_edge_info(e);
    }

    if( varray_find(star, v2) )
	// ?? BUG: Is it legitimate for there not to be an edge here ??
	create_edge(v1, v2);
}


void MxEdgeQSlim::apply_expansion(const MxPairContraction& conx)
{
    update_pre_expand(conx);

    m->apply_expansion(conx);

    //
    // Post-expansion update
    valid_verts++;
    valid_faces += conx.dead_faces.length();
    quadrics(conx.v1) -= quadrics(conx.v2);

    update_post_expand(conx);
}

bool MxEdgeQSlim::decimate(uint target)
{
    MxPairContraction local_conx;

    while( valid_faces > target )
    {
	MxQSlimEdge *info = (MxQSlimEdge *)heap.extract();
	if( !info ) { return false; }

	MxVertexID v1=info->v1, v2=info->v2;

	if( m->vertex_is_valid(v1) && m->vertex_is_valid(v2) )
	{
	    MxPairContraction& conx = local_conx;

	    m->compute_contraction(v1, v2, &conx, info->vnew);

	    if( will_join_only && conx.dead_faces.length()>0 ) continue;

	    if( contraction_callback )
		(*contraction_callback)(conx, -info->heap_key());
	    
	    apply_contraction(conx, m_bUpdateMSDM);
	}

	delete info;
    }

    return true;
}

// QError is the quadric distance of collapsed vertex to it's original position
bool MxEdgeQSlim::DecimateByError(double QError)
{
	MxPairContraction local_conx;
	double QECurr = 0;

	while( valid_faces > 0 )
	{
		// 求出目前哈希表中最小误差 [4/1/2012 Han]
		MxQSlimEdge *info = (MxQSlimEdge *)heap.top();
		Quadric Q = quadrics(info->v1);
		Quadric Q2 = quadrics(info->v2);
		Q += Q2;
		// 如果求平均距离的话，发现误差会急剧，导致简化过于剧烈
		QECurr = Q(info->vnew)/* / Q.nf*/;		// Get current error, average distance to all original planes
		
		////////////
		if (QECurr > QError)
			return true;
		
		decimate(valid_faces -1);

		//// test [3/29/2012 Han]
		//MxVertex v1Vert = m->vertex(info->v1);
		//info = (MxQSlimEdge *)heap.extract();
		//if( !info ) { return false; }

		//MxVertexID v1=info->v1, v2=info->v2;

		//if( m->vertex_is_valid(v1) && m->vertex_is_valid(v2) )
		//{
		//	MxPairContraction& conx = local_conx;

		//	m->compute_contraction(v1, v2, &conx, info->vnew);

		//	if( will_join_only && conx.dead_faces.length()>0 ) continue;

		//	if( contraction_callback )
		//		(*contraction_callback)(conx, -info->heap_key());

		//	apply_contraction(conx);

		//	// test [3/29/2012 Han]
		//	//MxVertex v1new = m->vertex(info->v1);
		//	//QECurr = pow(v1Vert[0] - v1new[0], 2)+ pow(v1Vert[1] - v1new[1], 2)+pow(v1Vert[2] - v1new[2], 2);
		//	//Q = quadrics(info->v1);
		//	//// 在此处计算边折叠带来的新点到原来点所在平面的误差 [3/29/2012 Han]
		//	//Quadric Q = quadrics(info->v1);
		//	//QECurr = Q(info->vnew);
		//}

		//delete info;
	}
	return false;
}


void MxFaceQSlim::compute_face_info(MxFaceID f)
{
    tri_info& info = f_info(f);
    info.f = f;

    MxVertexID i = m->face(f)(0);
    MxVertexID j = m->face(f)(1);
    MxVertexID k = m->face(f)(2);

    const Quadric& Qi = quadrics(i);
    const Quadric& Qj = quadrics(j);
    const Quadric& Qk = quadrics(k);

    Quadric Q = Qi;
    Q += Qj;
    Q += Qk;

    if( placement_policy == MX_PLACE_OPTIMAL &&
	Q.optimize(&info.vnew[X], &info.vnew[Y], &info.vnew[Z]) )
    {
      info.heap_key(-Q(info.vnew));
    }
    else
    {
      Vec3 v1(m->vertex(i)), v2(m->vertex(j)), v3(m->vertex(k));
      double e1=Q(v1), e2=Q(v2), e3=Q(v3);

      Vec3 best;
      double e_min;

      if( e1<=e2 && e1<=e3 ) { e_min=e1; best=v1; }
      else if( e2<=e1 && e2<=e3 ) { e_min=e2; best=v2; }
      else { e_min=e3; best=v3; }

      info.vnew[X] = best[X];
      info.vnew[Y] = best[Y];
      info.vnew[Z] = best[Z];
      info.heap_key(-e_min);
    }

    if( weighting_policy == MX_WEIGHT_AREA_AVG )
	info.heap_key(info.heap_key() / Q.area());

    if( info.is_in_heap() )
	heap.update(&info);
    else
	heap.insert(&info);
}


MxFaceQSlim::MxFaceQSlim(MxStdModel& _m)
  : MxQSlim(_m), f_info(_m.face_count())
{
}

void MxFaceQSlim::initialize()
{
    MxQSlim::initialize();

    for(MxFaceID f=0; f<m->face_count(); f++)
      compute_face_info(f);
}

// 进行边折叠 [11/30/2011 Han]
bool MxFaceQSlim::decimate(uint target)
{
    unsigned int i;

    MxFaceList changed;

    while( valid_faces > target )
    {
	tri_info *info = (tri_info *)heap.extract();
	if( !info ) { return false; }

	MxFaceID f = info->f;
	MxVertexID v1 = m->face(f)(0),
	  v2 = m->face(f)(1),
	  v3 = m->face(f)(2);

	if( m->face_is_valid(f) )
	{
	    //
	    // Perform the actual contractions
	    m->contract(v1, v2, v3, info->vnew, changed);

	    quadrics(v1) += quadrics(v2);  	// update quadric of v1
	    quadrics(v1) += quadrics(v3);

	    //
	    // Update valid counts
	    valid_verts -= 2;
	    for(i=0; i<changed.length(); i++)
		if( !m->face_is_valid(changed(i)) ) valid_faces--;

	    for(i=0; i<changed.length(); i++)
		if( m->face_is_valid(changed(i)) )
		    compute_face_info(changed(i));
		else
		    heap.remove(&f_info(changed(i)));
	}
    }

    return true;
}

// QError is the quadric distance of collapsed vertex to it's original position
bool MxFaceQSlim::DecimateByError(double QError)
{
	unsigned int i;

	MxFaceList changed;

	//while( valid_faces > target )
	//{
	//	tri_info *info = (tri_info *)heap.extract();
	//	if( !info ) { return false; }
	//	info->

	//	MxFaceID f = info->f;
	//	MxVertexID v1 = m->face(f)(0),
	//		v2 = m->face(f)(1),
	//		v3 = m->face(f)(2);

	//	if( m->face_is_valid(f) )
	//	{
	//		//
	//		// Perform the actual contractions
	//		m->contract(v1, v2, v3, info->vnew, changed);

	//		quadrics(v1) += quadrics(v2);  	// update quadric of v1
	//		quadrics(v1) += quadrics(v3);

	//		//
	//		// Update valid counts
	//		valid_verts -= 2;
	//		for(i=0; i<changed.length(); i++)
	//			if( !m->face_is_valid(changed(i)) ) valid_faces--;

	//		for(i=0; i<changed.length(); i++)
	//			if( m->face_is_valid(changed(i)) )
	//				compute_face_info(changed(i));
	//			else
	//				heap.remove(&f_info(changed(i)));
	//	}
	//}

	return true;
}
