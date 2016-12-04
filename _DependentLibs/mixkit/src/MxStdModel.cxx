/************************************************************************

  MxStdModel

  Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.
  
  $Id: MxStdModel.cxx,v 1.42 2000/11/20 20:36:38 garland Exp $

 ************************************************************************/

#include "stdmix.h"
#include "MxStdModel.h"
#include "MxVector.h"
#include "C:\_Hanhonglei\Projects\MeshSimp\jmspmesh\vec3.h"
#include "C:\_Hanhonglei\Projects\MeshSimp\MSDM2\MSDM2.h"

#include<set>
#include<map>
#include<list>
#include<vector>
#include<stack>
//#include <gl\gl.h>
//#include <gl\glu.h>
#include "C:\_Hanhonglei\Projects\MeshSimp\ViewSelectionBenchmark\Common\mesh.h"
#include "C:\_Hanhonglei\Projects\MeshSimp\ViewSelectionBenchmark/ViewSelect-OUR/ppm.h"

const static  double m_pi = 3.14159265359;

MxPairContraction& MxPairContraction::operator=(const MxPairContraction& c)
{
    v1 = c.v1;
    v2 = c.v2;
    mxv_set(dv1, c.dv1, 3);
    mxv_set(dv2, c.dv2, 3);

    delta_faces.reset();
    dead_faces.reset();

    for(uint i=0; i<c.delta_faces.length(); i++)
	delta_faces.add(c.delta_faces[i]);
    for(uint j=0; j<c.dead_faces.length(); j++)
	dead_faces.add(c.dead_faces[j]);

    delta_pivot = c.delta_pivot;

    return *this;
}

MxStdModel::~MxStdModel()
{
    for(uint i=0; i<face_links.length(); i++)  delete face_links[i];
	if(pFaceInfo != NULL)
		delete[] pFaceInfo;
	//if (m_pSrcMesh != NULL)
	//	delete m_pSrcMesh;
	//m_pSrcMesh = NULL;
}

MxVertexID MxStdModel::alloc_vertex(float x, float y, float z)
{
    MxVertexID id = MxBlockModel::alloc_vertex(x,y,z);
    v_data.add();
    v_data(id).tag = 0x0;
    v_data(id).user_tag = 0x0;
    vertex_mark_valid(id);

    face_links.add(new MxFaceList);
    unsigned int l = face_links.last_id();
    SanityCheck( l == id );
    SanityCheck( neighbors(id).length() == 0 );

    return id;
}

void MxStdModel::free_vertex(MxVertexID v)
{
    delete face_links[v];
    face_links.remove(v);
    v_data.remove(v);
}

MxFaceID MxStdModel::alloc_face(MxVertexID v1, MxVertexID v2, MxVertexID v3)
{
    MxFaceID id = MxBlockModel::alloc_face(v1,v2,v3);
    f_data.add();
    f_data(id).tag = 0x0;
    f_data(id).user_tag = 0x0;
    face_mark_valid(id);

    return id;
}

void MxStdModel::free_face(MxFaceID f)
{
    f_data.remove(f);
}

void MxStdModel::init_face(MxFaceID id)
{
    neighbors(face(id).v[0]).add(id);
    neighbors(face(id).v[1]).add(id);
    neighbors(face(id).v[2]).add(id);
}

MxStdModel *MxStdModel::clone()
{
    MxStdModel *m = new MxStdModel(vert_count(), face_count());
    // ??BUG: Current flags/marks won't be copied.  Is this the
    //        behavior we want?
    MxBlockModel::clone(m);

	// 克隆完以后将克隆后的模型指针设定为原始模型 [4/26/2012 Han]
	m_pSrcMesh = m;
    
    return m;
}

void MxStdModel::mark_neighborhood(MxVertexID vid, unsigned short mark)
{
    AssertBound( vid < vert_count() ); 

    for(unsigned int i=0; i<neighbors(vid).length(); i++)
    {
	unsigned int f = neighbors(vid)(i);
	fmark(f, mark);
    }
}

void MxStdModel::collect_unmarked_neighbors(MxVertexID vid,MxFaceList& faces)
{
    AssertBound( vid < vert_count() ); 

    for(unsigned int i=0; i<neighbors(vid).length(); i++)
    {
	unsigned int fid = neighbors(vid)(i);
	if( !fmark(fid) )
	{
	    faces.add(fid);
	    fmark(fid, 1);
	}
    }
}

void MxStdModel::mark_neighborhood_delta(MxVertexID vid, short delta)
{
    AssertBound( vid < vert_count() );
    for(uint i=0; i<neighbors(vid).length(); i++)
    {
	uint f = neighbors(vid)(i);
	fmark(f, fmark(f)+delta);
    }
}

void MxStdModel::partition_marked_neighbors(MxVertexID v, unsigned short pivot,
					    MxFaceList& lo, MxFaceList& hi)
{
    AssertBound( v < vert_count() );
    for(uint i=0; i<neighbors(v).length(); i++)
    {
	uint f = neighbors(v)(i);
	if( fmark(f) )
	{
	    if( fmark(f) < pivot )  lo.add(f);
	    else  hi.add(f);
	    fmark(f, 0);
	}
    }
}

void MxStdModel::mark_corners(const MxFaceList& faces, unsigned short mark)
{
    for(uint i=0; i<faces.length(); i++)
	for(uint j=0; j<3; j++)
	    vmark(face(faces(i))(j), mark);
}

void MxStdModel::collect_unmarked_corners(const MxFaceList& faces,
					  MxVertexList& verts)
{
    for(uint i=0; i<faces.length(); i++)
	for(uint j=0; j<3; j++)
	{
	    MxVertexID v = face(faces(i))(j);
	    if( !vmark(v) )
	    {
		verts.add(v);
		vmark(v, 1);
	    }
	}
}

void MxStdModel::collect_edge_neighbors(unsigned int v1, unsigned int v2,
					MxFaceList& faces)
{
    mark_neighborhood(v1, 1);
    mark_neighborhood(v2, 0);
    collect_unmarked_neighbors(v1, faces);
}

void MxStdModel::collect_vertex_star(unsigned int v, MxVertexList& verts)
{
    const MxFaceList& N = neighbors(v);

    mark_corners(N, 0);
    vmark(v, 1); // Don't want to include v in the star
    collect_unmarked_corners(N, verts);
}

void MxStdModel::collect_neighborhood(MxVertexID v, int depth,
				      MxFaceList& faces)
{
    // TODO: This method is somewhat inefficient.  It will repeatedly touch
    // all the faces within the collected region at each iteration.  For now,
    // this is acceptable.  But ultimately it will need to be fixed.

    int i;

    faces.reset();

    // Initially copy immediate neighbors of v
    for(i=0; i<neighbors(v).length(); i++)
	faces.add(neighbors(v)[i]);

    for(; depth>0; depth--)
    {
	// Unmark the neighbors of all vertices in region
	for(i=0; i<faces.length(); i++)
	{
	    mark_neighborhood(face(faces[i])[0], 0);
	    mark_neighborhood(face(faces[i])[1], 0);
	    mark_neighborhood(face(faces[i])[2], 0);
	}

	// Mark all faces already accumulated
	for(i=0; i<faces.length(); i++)
	    fmark(faces[i], 1);

	// Collect all unmarked faces
	uint limit = faces.length();
	for(i=0; i<limit; i++)
	{
	    collect_unmarked_neighbors(face(faces[i])[0], faces);
	    collect_unmarked_neighbors(face(faces[i])[1], faces);
	    collect_unmarked_neighbors(face(faces[i])[2], faces);
	}
    }
}


void MxStdModel::compute_vertex_normal(MxVertexID v, float *n)
{
    MxFaceList& star = neighbors(v);
    mxv_set(n, 0.0f, 3);

    unsigned int i;
    for(i=0; i<star.length(); i++)
    {
	float fn[3];

	// Weight normals uniformly
	compute_face_normal(star(i), fn, false);
	//
	// Weight normals by angle around vertex
	// 		uint c = face(star[i]).find_vertex(v);
	// 		compute_face_normal(star[i], fn);
	// 		mxv_scale(fn, compute_corner_angle(star[i], c), 3);

	mxv_addinto(n, fn, 3);
    }
    if( i>0 )
	mxv_unitize(n, 3);
}


void MxStdModel::synthesize_normals(uint start)
{
    float n[3];

    if( normal_binding() == MX_PERFACE )
    {
	for(MxFaceID f=start; f<face_count(); f++)
	{
	    compute_face_normal(f, n);
	    add_normal(n[X], n[Y], n[Z]);
	}
    }
    else if( normal_binding() == MX_PERVERTEX )
    {
	for(MxVertexID v=start; v<vert_count(); v++)
	{
	    compute_vertex_normal(v, n);
	    add_normal(n[X], n[Y], n[Z]);
	}
    }
    else
	mxmsg_signal(MXMSG_WARN, "Unsupported normal binding.",
		     "MxStdModel::synthesize_normals");
}



void MxStdModel::remap_vertex(unsigned int from, unsigned int to)
{
    AssertBound( from < vert_count() ); 
    AssertBound( to < vert_count() ); 
    SanityCheck( vertex_is_valid(from) );
    SanityCheck( vertex_is_valid(to) );
    
    for(unsigned int i=0; i<neighbors(from).length(); i++)
	face(neighbors(from)(i)).remap_vertex(from, to);

    mark_neighborhood(from, 0);
    mark_neighborhood(to, 1);
    collect_unmarked_neighbors(from, neighbors(to));

    vertex_mark_invalid(from);
    neighbors(from).reset();   // remove links in old vertex
}

unsigned int MxStdModel::split_edge(unsigned int a, unsigned int b)
{
    float *v1 = vertex(a), *v2 = vertex(b);

    return split_edge(a, b,
                      (v1[X] + v2[X])/2.0f,
                      (v1[Y] + v2[Y])/2.0f,
                      (v1[Z] + v2[Z])/2.0f);
}

static
void remove_neighbor(MxFaceList& faces, unsigned int f)
{
    unsigned int j;
    if( varray_find(faces, f, &j) )
	faces.remove(j);
}

unsigned int MxStdModel::split_edge(unsigned int v1, unsigned int v2,
			    float x, float y, float z)
{
    AssertBound( v1 < vert_count() );   AssertBound( v2 < vert_count() );
    SanityCheck( vertex_is_valid(v1) ); SanityCheck( vertex_is_valid(v2) );
    SanityCheck( v1 != v2 );

    MxFaceList faces;
    collect_edge_neighbors(v1, v2, faces);
    SanityCheck( faces.length() > 0 );

    unsigned int vnew = add_vertex(x,y,z);

    for(unsigned int i=0; i<faces.length(); i++)
    {
	unsigned int f = faces(i);
	unsigned int v3 = face(f).opposite_vertex(v1, v2);
	SanityCheck( v3!=v1 && v3!=v2 );
	SanityCheck( vertex_is_valid(v3) );

	// in f, remap v2-->vnew
	face(f).remap_vertex(v2, vnew);
	neighbors(vnew).add(f);

	// remove f from neighbors(v2)
	remove_neighbor(neighbors(v2), f);

	// assure orientation is consistent
	if( face(f).is_inorder(vnew, v3) )
	    add_face(vnew, v2, v3);
	else
	    add_face(vnew, v3, v2);
    }

    return vnew;
}

void MxStdModel::flip_edge(unsigned int v1, unsigned int v2)
{
    MxFaceList faces;  collect_edge_neighbors(v1, v2, faces);
    if( faces.length() != 2 ) return;

    unsigned int f1 = faces(0);
    unsigned int f2 = faces(1);
    unsigned int v3 = face(f1).opposite_vertex(v1, v2);
    unsigned int v4 = face(f2).opposite_vertex(v1, v2);

    // ?? Should we check for convexity or assume thats been taken care of?

    remove_neighbor(neighbors(v1), f2);
    remove_neighbor(neighbors(v2), f1);
    neighbors(v3).add(f2);
    neighbors(v4).add(f1);

    face(f1).remap_vertex(v2, v4);
    face(f2).remap_vertex(v1, v3);
}

void MxStdModel::split_face4(unsigned int f, unsigned int *newverts)
{
    unsigned int v0 = face(f).v[0];
    unsigned int v1 = face(f).v[1];
    unsigned int v2 = face(f).v[2];

    unsigned int pivot = split_edge(v0, v1);
    unsigned int new1 = split_edge(v1, v2);
    unsigned int new2 = split_edge(v0, v2);

    if( newverts )
    {
	newverts[0] = pivot;
	newverts[1] = new1;
	newverts[2] = new2;
    }

    flip_edge(pivot, v2);
}

void MxStdModel::compact_vertices()
{
    MxVertexID oldID;
    MxVertexID newID = 0;

    for(oldID=0; oldID<vert_count(); oldID++)
    {
	if( vertex_is_valid(oldID) )
	{
	    if( newID != oldID )
	    {
		vertex(newID) = vertex(oldID);
		if( normal_binding() == MX_PERVERTEX )
		    normal(newID)=normal(oldID);
		if( color_binding() == MX_PERVERTEX )
		    color(newID)=color(oldID);
		if( texcoord_binding() == MX_PERVERTEX )
		    texcoord(newID) = texcoord(oldID);

		// Because we'll be freeing the link lists for the
		// old vertices, we actually have to swap values instead
		// of the simple copying in the block above.
		//
		MxFaceList *t = face_links(newID);
		face_links(newID) = face_links(oldID);
		face_links(oldID) = t;

		vertex_mark_valid(newID);

		for(unsigned int i=0; i<neighbors(newID).length(); i++)
		    face(neighbors(newID)(i)).remap_vertex(oldID, newID);
	    }
	    newID++;
	}
    }

    for(oldID = newID; newID < vert_count(); )
	remove_vertex(oldID);
}

void MxStdModel::unlink_face(MxFaceID fid)
{
    MxFace& f = face(fid);
    face_mark_invalid(fid);

    unsigned int j; int found=0;

    if( varray_find(neighbors(f(0)), fid, &j) )
    {found++; neighbors(f(0)).remove(j);}
    if( varray_find(neighbors(f(1)), fid, &j) )
    { found++; neighbors(f(1)).remove(j); }
    if( varray_find(neighbors(f(2)), fid, &j) )
    { found++; neighbors(f(2)).remove(j); }
    SanityCheck( found > 0 );
    SanityCheck( !varray_find(neighbors(f(0)), fid, &j) );
    SanityCheck( !varray_find(neighbors(f(1)), fid, &j) );
    SanityCheck( !varray_find(neighbors(f(2)), fid, &j) );
}

void MxStdModel::remove_degeneracy(MxFaceList& faces)
{
    for(unsigned int i=0; i<faces.length(); i++)
    {
	SanityCheck( face_is_valid(faces(i)) );
	MxFace& f = face(faces(i));

	if( f(0)==f(1) || f(1)==f(2) || f(0)==f(2) )
	    unlink_face(faces(i));
    }
}

void MxStdModel::compute_contraction(MxVertexID v1, MxVertexID v2,
				     MxPairContraction *conx,
				     const float *vnew)
{
    conx->v1 = v1;
    conx->v2 = v2;

    if( vnew )
    {
	mxv_sub(conx->dv1, vnew, vertex(v1), 3);
	mxv_sub(conx->dv2, vnew, vertex(v2), 3);
    }
    else
    {
	conx->dv1[X] = conx->dv1[Y] = conx->dv1[Z] = 0.0;
	conx->dv2[X] = conx->dv2[Y] = conx->dv2[Z] = 0.0;
    }

    conx->delta_faces.reset();
    conx->dead_faces.reset();


    // Mark the neighborhood of (v1,v2) such that each face is
    // tagged with the number of times the vertices v1,v2 occur
    // in it.  Possible values are 1 or 2.
    //
    mark_neighborhood(v2, 0);
    mark_neighborhood(v1, 1);
    mark_neighborhood_delta(v2, 1);


    // Now partition the neighborhood of (v1,v2) into those faces
    // which degenerate during contraction and those which are merely
    // reshaped.
    //
    partition_marked_neighbors(v1, 2, conx->delta_faces, conx->dead_faces);
    conx->delta_pivot = conx->delta_faces.length();
    partition_marked_neighbors(v2, 2, conx->delta_faces, conx->dead_faces);
}

void MxStdModel::apply_contraction(const MxPairContraction& conx)
{
    MxVertexID v1=conx.v1, v2=conx.v2;

    // Move v1 to new position
    mxv_addinto(vertex(v1), conx.dv1, 3);

	// update V1's view selection importance [5/24/2012 Han]
	vertex(conx.v1).view_importance = (vertex(conx.v1).view_importance +vertex(conx.v2).view_importance)/2;


    uint i;
    //
    // Remove dead faces
    for(i=0; i<conx.dead_faces.length(); i++)
	unlink_face(conx.dead_faces(i));

    //
    // Update changed faces
    for(i=conx.delta_pivot; i<conx.delta_faces.length(); i++)
    {
	MxFaceID fid = conx.delta_faces(i);
	face(fid).remap_vertex(v2, v1);
	neighbors(v1).add(fid);
    }

    //
    // !!HACK: This is really only a temporary solution to the problem
    if( normal_binding() == MX_PERFACE )
    {
	float n[3];
	for(i=0; i<conx.delta_faces.length(); i++)
	{
	    compute_face_normal(conx.delta_faces[i], n);
	    normal(conx.delta_faces[i]) = MxNormal(n);
	}
    }

    //
    // Kill v2
    vertex_mark_invalid(v2);
    neighbors(v2).reset();
}

void MxStdModel::apply_expansion(const MxPairExpansion& conx)
{
    MxVertexID v1=conx.v1, v2=conx.v2;

    mxv_sub(vertex(v2), vertex(v1), conx.dv2, 3);
    mxv_subfrom(vertex(v1), conx.dv1, 3);

    uint i,j;
    for(i=0; i<conx.dead_faces.length(); i++)
    {
	MxFaceID fid = conx.dead_faces(i);
	face_mark_valid(fid);
	neighbors(face(fid)(0)).add(fid);
	neighbors(face(fid)(1)).add(fid);
	neighbors(face(fid)(2)).add(fid);
    }

    for(i=conx.delta_pivot; i<conx.delta_faces.length(); i++)
    {
	MxFaceID fid = conx.delta_faces(i);
	face(fid).remap_vertex(v1, v2);
	neighbors(v2).add(fid);
	bool found = varray_find(neighbors(v1), fid, &j);
	SanityCheck( found );
	neighbors(v1).remove(j);
    }

    //
    // !!HACK: This is really only a temporary solution to the problem
    if( normal_binding() == MX_PERFACE )
    {
	float n[3];
	for(i=0; i<conx.delta_faces.length(); i++)
	{
	    compute_face_normal(conx.delta_faces[i], n);
	    normal(conx.delta_faces[i]) = MxNormal(n);
	}

	for(i=0; i<conx.dead_faces.length(); i++)
	{
	    compute_face_normal(conx.dead_faces[i], n);
	    normal(conx.dead_faces[i]) = MxNormal(n);
	}
    }

    vertex_mark_valid(v2);
}


void MxStdModel::contract(MxVertexID v1, MxVertexID v2,
			  const float *vnew, MxPairContraction *conx)
{
    compute_contraction(v1, v2, conx);
    mxv_sub(conx->dv1, vnew, vertex(v1), 3);
    mxv_sub(conx->dv2, vnew, vertex(v2), 3);
    apply_contraction(*conx);
}

void MxStdModel::compute_contraction(MxFaceID fid, MxFaceContraction *conx)
{
    const MxFace& f = face(fid);

    conx->f = fid;
    conx->dv1[X] = conx->dv1[Y] = conx->dv1[Z] = 0.0;
    conx->dv2[X] = conx->dv2[Y] = conx->dv2[Z] = 0.0;
    conx->dv3[X] = conx->dv3[Y] = conx->dv3[Z] = 0.0;

    conx->delta_faces.reset();
    conx->dead_faces.reset();


    PARANOID( mark_neighborhood(f[0], 0) );
    mark_neighborhood(f[1], 0);
    mark_neighborhood(f[2], 0);

    mark_neighborhood(f[0], 1);
    mark_neighborhood_delta(f[1], +1);
    mark_neighborhood_delta(f[2], +1);

    fmark(fid, 0);		// don't include f in dead_faces

    partition_marked_neighbors(f[0], 2, conx->delta_faces, conx->dead_faces);
    partition_marked_neighbors(f[1], 2, conx->delta_faces, conx->dead_faces);
    partition_marked_neighbors(f[2], 2, conx->delta_faces, conx->dead_faces);
}

void MxStdModel::contract(MxVertexID v1, MxVertexID v2, MxVertexID v3,
			  const float *vnew,
			  MxFaceList& changed)
{
    mark_neighborhood(v1, 0);
    mark_neighborhood(v2, 0);
    mark_neighborhood(v3, 0);
    changed.reset();
    collect_unmarked_neighbors(v1, changed);
    collect_unmarked_neighbors(v2, changed);
    collect_unmarked_neighbors(v3, changed);

    // Move v1 to vnew
    vertex(v1)(0) = vnew[X];
    vertex(v1)(1) = vnew[Y];
    vertex(v1)(2) = vnew[Z];

    // Replace occurrences of v2 & v3 with v1
    remap_vertex(v2, v1);
    remap_vertex(v3, v1);

    remove_degeneracy(changed);

    //
    // !!HACK: Only a temporary solution
    if( normal_binding() == MX_PERFACE )
    {
	float n[3];
	for(uint i=0; i<changed.length(); i++)
	    if( face_is_valid(changed[i]) )
	    {
		compute_face_normal(changed[i], n);
		normal(changed[i]) = MxNormal(n);
	    }
    }
}

void MxStdModel::contract(MxVertexID v1, const MxVertexList& rest,
			  const float *vnew, MxFaceList& changed)
{
    uint i;

    // Collect all effected faces
    mark_neighborhood(v1, 0);
    for(i=0; i<rest.length(); i++)
	mark_neighborhood(rest(i), 0);

    changed.reset();

    collect_unmarked_neighbors(v1, changed);
    for(i=0; i<rest.length(); i++)
	collect_unmarked_neighbors(rest(i), changed);

    // Move v1 to vnew
    vertex(v1)(0) = vnew[X];
    vertex(v1)(1) = vnew[Y];
    vertex(v1)(2) = vnew[Z];

    // Replace all occurrences of vi with v1
    for(i=0; i<rest.length(); i++)
	remap_vertex(rest(i), v1);

    remove_degeneracy(changed);
}

MxVertexID MxStdModel::resolve_proxies(MxVertexID v)
{
    while( vertex_is_proxy(v) )
	v = vertex(v).as.proxy.parent;
    return v;
}

float *MxStdModel::vertex_position(MxVertexID v)
{
    return vertex(resolve_proxies(v));
}

MxVertexID& MxStdModel::vertex_proxy_parent(MxVertexID v)
{
    SanityCheck( vertex_is_proxy(v) );
    return vertex(v).as.proxy.parent;
}

MxVertexID MxStdModel::add_proxy_vertex(MxVertexID parent)
{
    MxVertexID v = alloc_vertex(0, 0, 0); // position will be ignored

    vertex_mark_proxy(v);
    vertex_proxy_parent(v) = parent;

    return v;
}

// msdm2 [4/18/2012 Han]
void MxStdModel::set_index_vertices()
{
	int index = 0;
	for (MxVertexID id= 0; id < vert_count(); id++)
		vertex(id).tag(id);
}


// from mepp [4/15/2012 Han]
//**********************************************
// compute v x v^t
//**********************************************
void vector_times_transpose_mult(double pVector[3],
                                               double ppMatrix[3][3],
                                               double coeff)
{
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      ppMatrix[i][j] = coeff * pVector[i] * pVector[j];
}

//**********************************************
// fix sine
//**********************************************
double fix_sine(double sine)
{
  if (sine >= 1)
    return m_pi/2;
  else
    if (sine <= -1)
      return -m_pi/2;
    else
      return std::asin(sine);
}

//**********************************************
// add two matrices
//**********************************************
void add(double pMatrix[3][3],
                       double pMatrixSum[3][3])
{
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      pMatrixSum[i][j] += pMatrix[i][j];
}

//**********************************************
// principal curvature for a vertex
//**********************************************
void MxStdModel::principal_curvature_per_vert(MxVertexID verID, double ppMatrix_sum[3][3])
{
	double area=0;
	MxVertex &pVertex = vertex(verID);
  // iterate over all edges
  //Halfedge_around_vertex_circulator pHalfedge = pVertex.vertex_begin();
  //Halfedge_around_vertex_circulator pHalfedgeStart = pHalfedge;
  //CGAL_For_all(pHalfedge,pHalfedgeStart)
	MxVertexList star;
	star.reset();
	collect_vertex_star(verID, star);
	for(uint j=0; j<star.length(); j++)
  {

    // build edge vector and comput its norm
	  //Point3d p1 = pHalfedge->vertex()->point();
	  //Point3d p2 = pHalfedge->opposite()->vertex()->point();
	  MxVertex &vOpposite = vertex(star(j));

	  //Point3d p1 = h->vertex()->point();
	  //Point3d p2 = h->opposite()->vertex()->point();
	  vec3 p1 = vec3(pVertex.as.pos);
	  vec3 p2 = vec3(vOpposite.as.pos);

	  vec3 edge = (p1-p2);
		double len_edge = std::sqrt(edge.dot(edge));
		if (len_edge == 0) // avoid divide by zero
		continue;

    // compute (signed) angle between two incident faces, if exists
    //Facet_handle pFacet1 = pHalfedge->facet();
    //Facet_handle pFacet2 = pHalfedge->opposite()->facet();
    //CGAL_assertion(pFacet1 != pFacet2);
    //if (pFacet1 == NULL || pFacet2 == NULL)
    //  continue; // border edge

	MxFaceList edgeN;
	edgeN.reset();
	collect_edge_neighbors(verID, star(j), edgeN);
	float n1[3], n2[3];
	if (edgeN.length() < 2)
		continue;
	//这里还是有问题，简化以后，某些情况下会出现一个边有4个相邻面片的情况，而实际应该只有2个。
	// 目前的处理是对其求取平均Normal [7/7/2012 Han]
	else if (edgeN.length() > 2)
	{
		float nt[3];
		for (int iii = 0; iii < 3; iii++)
		{
			n1[iii] = 0.0f;
			n2[iii] = 0.0f;
		}

		for(int ttt = 0; ttt < edgeN.length(); ttt++)
		{
			compute_face_normal(edgeN(ttt), nt);
			if (face(edgeN(ttt)).is_inorder(verID, star(j)))
			{
				area += compute_face_area(edgeN(ttt));	//  [7/17/2012 Han]
				for (int iii = 0; iii < 3; iii++)
					n1[iii] += nt[iii];
			}
			else
				for (int iii = 0; iii < 3; iii++)
					n2[iii] += nt[iii];
		}
		float a = float(sqrt(n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2])); if (a!=0.0f) {n1[0]/=a; n1[1]/=a; n1[2]/=a;};
		a = float(sqrt(n2[0]*n2[0] + n2[1]*n2[1] + n2[2]*n2[2])); if (a!=0.0f) {n2[0]/=a; n2[1]/=a; n2[2]/=a;};
	}
	else
	{
		// 暂时不进行等于2检测 [7/8/2012 Han]
		//assert(edgeN.length()==2);
		MxFaceID pFacet1 = edgeN(0);
		MxFaceID pFacet2 = edgeN(1);
		assert(pFacet1 != pFacet2);

		// 边（v1， v2）的左手边应该是pFacet1，右手边是pFacet2
		if (!face(pFacet1).is_inorder(verID, star(j)))
		{
			MxFaceID tmp = pFacet1;
			pFacet1 = pFacet2;
			pFacet2 = tmp;
		}

		//area+=AreaFacetTriangle(pFacet1);
		area += compute_face_area(pFacet1);

		compute_face_normal(pFacet1, n1);
		compute_face_normal(pFacet2, n2);
	}
	//Vector normal1 = pFacet1->normal();
	//Vector normal2 = pFacet2->normal();
	vec3 normal1(n1);
	vec3 normal2(n2);

    double sine = ((normal1.cross(normal2)).dot(edge))/len_edge;
    double beta = fix_sine(sine);

    // compute edge * edge^t * coeff, and add it to current matrix
    double pVector_edge[3] = {edge.x,edge.y,edge.z};
    double ppMatrix[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    vector_times_transpose_mult(pVector_edge,ppMatrix,beta/len_edge);
    add(ppMatrix,ppMatrix_sum);
  }
	if (area == 0)
		area = 1;
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
			ppMatrix_sum[i][j]/=area;
}

// Return all border info if the vertex is on an edge of the mesh.
bool MxStdModel::IsBorderEdges(MxVertexID v1, MxVertexID v2)
{
	// So go through the list of all neighboring vertices, and see how many
	// triangles this vertex has in common w/ each neighboring vertex.  Normally
	// there will be two triangles in common, but if there is only one, then this 
	// vertex is on an edge.
	MxFaceList faces;
	faces.reset();
	collect_edge_neighbors(v1, v2, faces);
	if( faces.length() == 1 )
		return true;
	else
		return false;					

}

bool MxStdModel::sphere_clip_vector(vec3 &O, double r,const vec3 &P, vec3 &V)
{

	vec3 W = P - O ;
	double a = (V.dot(V));
	double b = 2.0 * V.dot(W) ;
	double c = (W.dot(W)) - r*r ;
	double delta = b*b - 4*a*c ;
	if (delta < 0) {
		// Should not happen, but happens sometimes (numerical precision)
		return true ;
	}
	double t = (- b + ::sqrt(delta)) / (2.0 * a) ;
	if (t < 0.0) {
		// Should not happen, but happens sometimes (numerical precision)
		return true ;
	}
	if (t >= 1.0) {
		// Inside the sphere
		return false ;
	}

	V=V*t;

	return true ;
}

void MxStdModel::geodes_principal_curvature_per_vert(MxVertexID verID, double ppMatrix_sum[3][3], double radius)
{ 

	MxVertex* pVertex = &vertex(verID);
	
	std::set<MxVertexID> verts ;
	vec3 O = vec3(vertex_position(verID));
	//std::stack<MxVertex*> S ;
	std::stack<MxVertexID> S ;

	S.push(verID) ;
	verts.insert(verID) ;
	int iter=0;
	while(!S.empty())
	{
		//Vertex* v = S.top() ;
		MxVertexID v = S.top();
		// 当前顶点 [4/18/2012 Han]
		MxVertex &vCurr = vertex(v);

		S.pop() ;
		//Point3d P = v->point() ;
		vec3 P = vec3(vCurr.as.pos);

		// halfedge around vertex指的是以点v为顶点的所有边，每个边都包括v [4/23/2012 Han]
		//Halfedge_around_vertex_circulator h = v->vertex_begin();
		//Halfedge_around_vertex_circulator pHalfedgeStart = h;
		//CGAL_For_all(h,pHalfedgeStart)
		MxVertexList star;
		star.reset();
		collect_vertex_star(v, star);
		//MxEdgeList h;
		//edge_around_vertex_circular(v, h);
		//for(uint j=0; j<h.length(); j++)
		for(int j = 0; j < star.length(); j++)
		{
			//Point3d p1 = h->vertex()->point();
			//Point3d p2 = h->opposite()->vertex()->point();
			//vec3 p1 = vec3(vertex_position(h[j].v2));
			//vec3 p2 = vec3(vertex_position(h[j].v1));;
			vec3 p1 = vec3(vCurr.as.pos);
			vec3 p2 = vec3(vertex_position(star(j)));

			vec3 V = (p2-p1);
			if (v==verID || V.dot(P - O) > 0.0)
			{
				//double len_old = std::sqrt(V*V);
				bool isect = sphere_clip_vector(O, radius, P, V) ;

				if (/*!h->is_border_edge()*//*!IsBorderEdges(h[j].v1, h[j].v2)*/!IsBorderEdges(v, star(j)))
				{
					double len_edge = std::sqrt(V.dot(V));
					// compute (signed) angle between two incident faces, if exists
					MxFaceList edgeN;
					edgeN.reset();
					//collect_edge_neighbors(h[j].v1, h[j].v2, edgeN);
					collect_edge_neighbors(v, star(j), edgeN);

					float n1[3], n2[3];

				//这里还是有问题，简化以后，某些情况下会出现一个边有4个相邻面片的情况，而实际应该只有2个。
					// 目前的处理是对其求取平均Normal [7/7/2012 Han]
					if (edgeN.length() > 2)
					{
						float nt[3];
						for (int iii = 0; iii < 3; iii++)
						{
							n1[iii] = 0.0f;
							n2[iii] = 0.0f;
						}
						
						for(int ttt = 0; ttt < edgeN.length(); ttt++)
						{
							compute_face_normal(edgeN(ttt), nt);
							if (face(edgeN(ttt)).is_inorder(v, star(j)))
								for (int iii = 0; iii < 3; iii++)
									n1[iii] += nt[iii];
							else
								for (int iii = 0; iii < 3; iii++)
									n2[iii] += nt[iii];
						}
						float a = float(sqrt(n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2])); if (a!=0.0f) {n1[0]/=a; n1[1]/=a; n1[2]/=a;};
						a = float(sqrt(n2[0]*n2[0] + n2[1]*n2[1] + n2[2]*n2[2])); if (a!=0.0f) {n2[0]/=a; n2[1]/=a; n2[2]/=a;};
					}
					else
					{
						// 暂时不进行等于2检测 [7/8/2012 Han]
						//assert(edgeN.length()==2);

						//Facet_handle pFacet1 = h->facet();
						//Facet_handle pFacet2 = h->opposite()->facet();
						MxFaceID pFacet1 = edgeN(0);
						MxFaceID pFacet2 = edgeN(1);
						assert(pFacet1 != pFacet2);
						// 边（v1， v2）的左手边应该是pFacet1，右手边是pFacet2
						//if (!face(pFacet1).is_inorder(h[j].v1, h[j].v2))
						if (!face(pFacet1).is_inorder(v, star(j)))
						{
							MxFaceID tmp = pFacet1;
							pFacet1 = pFacet2;
							pFacet2 = tmp;
						}

						compute_face_normal(pFacet1, n1);
						compute_face_normal(pFacet2, n2);
					}
					vec3 normal1(n1);
					vec3 normal2(n2);

					double sine = ((normal1.cross(normal2)).dot(V))/len_edge;
					double beta = fix_sine(sine);
					// compute edge * edge^t * coeff, and add it to current matrix
					double pVector_edge[3] = {V.x,V.y,V.z};
					double ppMatrix[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
					vector_times_transpose_mult(pVector_edge,ppMatrix,beta/len_edge);
					add(ppMatrix,ppMatrix_sum);
				}
				if (!isect) {
					//Vertex_iterator w=h->opposite()->vertex();
					//MxVertexID w = h[j].v1;
					MxVertexID w = star(j);
					//if (verts.find(&(*w)) == verts.end()) {
					//    verts.insert(&(*w)) ;
					//    S.push(&(*w)) ;
					//}
					if (verts.find(w) == verts.end()) {
						verts.insert(w) ;
						S.push(w) ;
					}

				}
			}
		}
		iter++;
	}

	double area=m_pi*radius*radius;
	for (int i=0;i<3;i++)
		for (int j=0;j<3;j++)
			ppMatrix_sum[i][j]/=area;
}


#define ABS_GUY  fabs
#define SQRT_GUY sqrt
#define MACH_EPS 2e-16


#define NB_VP_PRIS 10

/*********************************************************************
*      Eigenvalues and eigenvectors of a double symmetric matrix     *
*                      by JACOBI'S METHOD                            *
* ------------------------------------------------------------------ *
* INPUTS:                                                            *
*           A(N,N) : double symmetric matrix                         *
*           N      : NMAX of matrix A                                *
*           ER     : desired pr閏ision                               *
*           IM     : maximum number of iterations                    *
* OUTPUTS:                                                           *
*           AD(N,N): eigenvalues stored in ascending order on main   *
*                    diagonal                                        *
*           U(N,N) : eigenvectors given in lines                     *
*                                                                    *
*                          C++ version from FORTRAN by J-P Moreau    *
* ------------------------------------------------------------------ *
* NOTA: index 0 not used here.                                       *
*********************************************************************/
int ValPro(int N,double **A,double ER,int IM,double**U,double **AD)  {
	double XNA,XND,ST,C,S,P,PR,Q,B,TE;
	int I,IT,J,K,NI,NJ=0;
	NI=-1;
	for (I=1; I<N+1; I++)
		for (J=1; J<N+1; J++) {
			U[I][J]=0.0;
			AD[I][J]=0.0;
		}
		IT=0;
		XNA=0.0;
		XND=0.0;
		ST=0.0;
		for (I=1; I<N+1; I++) {
			U[I][I]=1.0;
			XND=XND+A[I][I]*A[I][I];
			if (I==N) goto s10;
			for (J=I+1; J<N+1; J++) {
				ST=ST+(A[I][J]*A[I][J]);
				U[I][J]=0.0;
				U[J][I]=0.0;
			}
		}
s10: XNA=XND+2.0*ST;
		do {
			ST=0.0;
			for (I=1; I<N; I++) {
				for (J=I+1; J<N+1; J++) {
					TE=ABS_GUY(A[I][J]);
					if (TE > ST) {
						ST=TE;
						NI=I;
						NJ=J;
					}
				}
			}
			if (NI==-1)
				return -1;
			B=A[NI][NI]-A[NJ][NJ];
			if (ABS_GUY(B) > MACH_EPS) goto s15;
			C=1.0/SQRT_GUY(2.0);

			if (A[NI][NJ]>0) S=C;
			else if (A[NI][NJ]<0) S=-C;
			else S=0;
			// S=C*DSIGN(0.1D+01,A[NI][NJ]) fortran
			goto s20;
s15: Q=ABS_GUY(B);
			if (B<0) P=-2.0*A[NI][NJ];
			else if (B>0) P=2.0*A[NI][NJ];
			else P=0.0;
			// P=2.0*A[NI][NJ]*DSIGN(0.1D+01,B)  fortran
			ST=SQRT_GUY(P*P+Q*Q);
			C=SQRT_GUY((1.0+Q/ST)/2.0);
			S=P/(2.0*ST*C);
s20: for (K=1; K<N+1; K++) {
			ST=U[K][NI];
			U[K][NI]=C*ST+S*U[K][NJ];
			U[K][NJ]=C*U[K][NJ]-S*ST;
	 }
	 if (NI==1) goto s30;
	 for (K=1; K<NI; K++) {
		 ST=A[K][NI];
		 A[K][NI]=C*ST+S*A[K][NJ];
		 A[K][NJ]=C*A[K][NJ]-S*ST;
	 }
s30: if (NJ==NI+1) goto s40;
	 for (K=NI+1; K<NJ; K++) {
		 ST=A[NI][K];
		 A[NI][K]=C*ST+S*A[K][NJ];
		 A[K][NJ]=C*A[K][NJ]-S*ST;
	 }
s40: if (NJ==N) goto s50;
	 for (K=NJ+1; K<N+1; K++) {
		 ST=A[NI][K];
		 A[NI][K]=C*ST+S*A[NJ][K];
		 A[NJ][K]=C*A[NJ][K]-S*ST;
	 }
s50: XND=XND+2.0*A[NI][NJ]*A[NI][NJ];
	 ST=A[NI][NI];
	 A[NI][NI]=C*C*ST+2.0*S*C*A[NI][NJ]+S*S*A[NJ][NJ];
	 A[NJ][NJ]=C*C*A[NJ][NJ]+S*S*ST-2.0*S*C*A[NI][NJ];
	 A[NI][NJ]=0.0;
	 IT=IT+1;
	 PR=ABS_GUY(1.0-XNA/XND);
	 if (PR<ER) goto s60;
		} while (IT <= IM);

s60:
		for (I=1; I<N+1; I++)
			for (J=1; J<N+1; J++)
				AD[I][J]=A[I][J];
		return 0;
} // valpro

/*********************************************************************
* Given the eigenvalues D and eigenvectors V as output from VALPRO,  *
* this routine sorts the eigenvalues into ascending order, and       *
* rearranges the lines of V accordingly. The method used is straight *
* insertion.                                                         *
*********************************************************************/
void EigSrt(double **D,double **V,int N) {
	double P;
	int I,J,K;
	for (I=1; I<N; I++) {
		K=I;
		P=D[I][I];
		for (J=I+1; J<N+1; J++)
			if (D[J][J] >= P)  {
				K=J;
				P=D[J][J];
			}

			if (K!=I)  {
				D[K][K]=D[I][I];
				D[I][I]=P;
				for (J=1; J<N+1; J++)  {
					P=V[J][I];
					V[J][I]=V[J][K];
					V[J][K]=P;
				}
			}
	}
}
/*
void MATREAD(double A[NMAX][NMAX],int N) {
	int I,J; double temp;
    for (I=1; I<N+1; I++)
		for (J=1; J<N+1; J++) {
			fscanf(fp1,"%lg",&temp);
			A[I][J] = temp;
        }
}
void MATPRINT(char *titre,double A[NMAX][NMAX],int N) {
	int I, J;
    printf("\n\n  %s  \n\n",titre);
    for (I=1; I<N+1; I++)
		for (J=1; J<N+1; J++)  {
			printf("%13.6f",A[I][J]);
			if (J==N) printf("\n");
		}
}
int WriteHead (FILE *fp, char * string){
	if (string == NULL) return (-2);
	if (fprintf (fp,"\n%s\n%s\n%s\n\n", Separator, string, Separator) <= 0)
		return (-1);
	return 0;
}
int WriteEnd (FILE *fp){
	if (fprintf (fp,"\n\n%s\n", Separator) <= 0) return (-1);
	 return 0;
}  */

// 得到某个顶点的主曲率
void MxStdModel::principal_curvature(bool IsGeod,double radius, int nScale, MxVertexID vID)
{
	nScale == nScale < 0 ? 0:nScale;
	assert(nScale <= MSDM_SCALE);

	if(!vertex_is_valid(vID) )
		return;

	MxVertex &pVertex = vertex(vID);
	bool NoValPro=false;

	double ppMatrix_sum[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	//double eigenvalues[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
	//double eigenvectors[3][3] = {{0,0,0},{0,0,0},{0,0,0}};

	if (IsGeod==true)//geodesic neighborhood
		geodes_principal_curvature_per_vert(vID,ppMatrix_sum,radius);
	else//1-ring neighborhood
		principal_curvature_per_vert(vID,ppMatrix_sum);

	//Eigen values/vectors
	double **CovMat=(double**)malloc((4)*sizeof(double*));
	double **VectPro=(double**)malloc((4)*sizeof(double*));
	double **Valpro=(double**)malloc((4)*sizeof(double*));
	for (int i=0;i<(4);i++)
	{
		CovMat[i]=(double*)malloc((4)*sizeof(double));
		VectPro[i]=(double*)malloc((4)*sizeof(double));
		Valpro[i]=(double*)malloc((4)*sizeof(double));
	}

	CovMat[1][1]=ppMatrix_sum[0][0];
	CovMat[1][2]=ppMatrix_sum[0][1];
	CovMat[1][3]=ppMatrix_sum[0][2];
	CovMat[2][1]=ppMatrix_sum[1][0];
	CovMat[2][2]=ppMatrix_sum[1][1];
	CovMat[2][3]=ppMatrix_sum[1][2];
	CovMat[3][1]=ppMatrix_sum[2][0];
	CovMat[3][2]=ppMatrix_sum[2][1];
	CovMat[3][3]=ppMatrix_sum[2][2];


	if (ppMatrix_sum[0][1]==0 && ppMatrix_sum[0][2]==0 &&
		ppMatrix_sum[1][0]==0 && ppMatrix_sum[1][2]==0 &&
		ppMatrix_sum[2][1]==0 && ppMatrix_sum[2][0]==0)
	{
		for (int i=1;i<4;i++)
			for (int j=1;j<4;j++)
			{
				Valpro[i][j]=CovMat[i][j];
				if (i==j)
					VectPro[i][j]=1;
				else
					VectPro[i][j]=0;
			}
	}
	else
	{
		if (ValPro(3,CovMat,1e-15,10000.,VectPro,Valpro)==-1)
		{
			NoValPro=true;

		}
	}

	if(!NoValPro)
	{
		//  Call the Jacovi subroutine
		for (int u=0;u<4;u++)
			for (int v=0;v<4;v++)
			{
				Valpro[u][v]=fabs(Valpro[u][v]);

			}
			EigSrt(Valpro,VectPro,3);
			vec3 VKmaxCurv(VectPro[1][2],VectPro[2][2],VectPro[3][2]);
			vec3 VKminCurv(VectPro[1][1],VectPro[2][1],VectPro[3][1]);

			// MT
			/*eigenvalues[0][0]=Valpro[1][1];
			eigenvalues[0][1]=Valpro[1][2];
			eigenvalues[0][2]=Valpro[1][3];
			eigenvalues[1][0]=Valpro[2][1];
			eigenvalues[1][1]=Valpro[2][2];
			eigenvalues[1][2]=Valpro[2][3];
			eigenvalues[2][0]=Valpro[3][1];
			eigenvalues[2][1]=Valpro[3][2];
			eigenvalues[2][2]=Valpro[3][3];

			eigenvectors[0][0]=VectPro[1][1];
			eigenvectors[0][1]=VectPro[1][2];
			eigenvectors[0][2]=VectPro[1][3];
			eigenvectors[1][0]=VectPro[2][1];
			eigenvectors[1][1]=VectPro[2][2];
			eigenvectors[1][2]=VectPro[2][3];
			eigenvectors[2][0]=VectPro[3][1];
			eigenvectors[2][1]=VectPro[3][2];
			eigenvectors[2][2]=VectPro[3][3];*/

			pVertex.VKmaxCurv[nScale]=VKmaxCurv;
			pVertex.VKminCurv[nScale]=VKminCurv;

			pVertex.KmaxCurv[nScale]=Valpro[1][1];
			pVertex.KminCurv[nScale]=Valpro[2][2];
	}
	else
	{
		pVertex.VKmaxCurv[nScale]=vec3(0, 0, 0)/*CGAL::NULL_VECTOR*/;
		pVertex.VKminCurv[nScale]=vec3(0, 0, 0)/*CGAL::NULL_VECTOR*/;
		pVertex.KmaxCurv[nScale]=0;
		pVertex.KminCurv[nScale]=0;
	}
	// 原文3，此处应为4 [5/2/2012 Han]
	for (int i=0;i</*(3)*/4;i++)
	{
		free(CovMat[i]);
		free(VectPro[i]);
		free(Valpro[i]);
	}
	free(CovMat);
	free(VectPro);
	free(Valpro);

#ifdef _MSC_VER
	MinNrmMinCurvature[nScale]=min(MinNrmMinCurvature[nScale],pVertex.KminCurv[nScale]);
	MaxNrmMinCurvature[nScale]=max(MaxNrmMinCurvature[nScale],pVertex.KminCurv[nScale]);

	MinNrmMaxCurvature[nScale]=min(MinNrmMaxCurvature[nScale],pVertex.KmaxCurv[nScale]);
	MaxNrmMaxCurvature[nScale]=max(MaxNrmMaxCurvature[nScale],pVertex.KmaxCurv[nScale]);
#else
	MinNrmMinCurvature=CGAL::min(MinNrmMinCurvature,pVertex->KminCurv[nScale]);
	MaxNrmMinCurvature=CGAL::max(MaxNrmMinCurvature,pVertex->KminCurv[nScale]);

	MinNrmMaxCurvature=CGAL::min(MinNrmMaxCurvature,pVertex->KmaxCurv[nScale]);
	MaxNrmMaxCurvature=CGAL::max(MaxNrmMaxCurvature,pVertex->KmaxCurv[nScale]);
#endif
}

void MxStdModel::principal_curvature(bool IsGeod,double radius, int nScale)
{
	MinNrmMinCurvature[nScale]=100000;
	MaxNrmMinCurvature[nScale]=-100000;

	MinNrmMaxCurvature[nScale]=100000;
	MaxNrmMaxCurvature[nScale]=-100000;

	for(MxVertexID vID=0; vID<vert_count(); vID++)
		principal_curvature(IsGeod, radius, nScale, vID);
}

void MxStdModel::ComputeMaxMin()
{
	MinMSDM2=FLT_MAX;
	MaxMSDM2=0;
	MinImportanceV = FLT_MAX;
	MaxImportanceV = 0;
	for(int i = 0; i < MSDM_SCALE; i++)
	{
		MinCurv[i] = FLT_MAX, MaxCurv[i] = 0;
	}

	for(MxVertexID vID=0; vID<vert_count(); vID++)
	{
		if(!vertex_is_valid(vID) )
			continue;
		MxVertex &pVertex = vertex(vID);
		if(pVertex.MSDM2_Local_sum>MaxMSDM2)
			MaxMSDM2=pVertex.MSDM2_Local_sum;
		if(pVertex.MSDM2_Local_sum<MinMSDM2)
			MinMSDM2=pVertex.MSDM2_Local_sum;

		if(pVertex.view_importance>MaxImportanceV)
			MaxImportanceV=pVertex.view_importance;
		if(pVertex.view_importance<MinImportanceV)
			MinImportanceV=pVertex.view_importance;

		for(int i = 0; i < MSDM_SCALE; i++)
		{
			if(pVertex.KmaxCurv[i]>MaxCurv[i])
				MaxCurv[i]=pVertex.KmaxCurv[i];
			if(pVertex.KmaxCurv[i]<MinCurv[i])
				MinCurv[i]=pVertex.KmaxCurv[i];
		}
	}
}


void MxStdModel::ComputeMaxMin(double *maxV, double *minV, int type)
{
	*maxV = 0, *minV = FLT_MAX;

	ComputeMaxMin();

	switch (type)
	{
	case 1:			// color by msdm
		*maxV = MaxMSDM2;
		*minV = MinMSDM2;
		break;
	case 2:			// color by max curvature
		*maxV = MaxCurv[0];
		*minV = MinCurv[0];
		break;
	case 3:			// color by importance ( ei. mesh saliency)
		*maxV = MaxImportanceV;
		*minV = MinImportanceV;
		break;
	default:
		break;
	}

}

//type:	case 1:			// color by msdm
//		case 2:			// color by max curvature
//		case 3:			// color by importance ( ei. mesh saliency)
void MxStdModel::VertexColor(int type)
{
	double maxV = 1.0, minV = 0.0;
	ComputeMaxMin(&maxV, &minV, type);
	VertexColor(type, maxV, minV);
}


//type:	case 1:			// color by msdm
//		case 2:			// color by max curvature
//		case 3:			// color by importance ( ei. mesh saliency)
void MxStdModel::VertexColor(int type, double maxV, double minV)
{
	double len = 1.0/(maxV - minV);

	//v_infor_colors.room_for(vert_count());

	double R = 0.0;
	int indiceLut = 0;

	for(MxVertexID vID=0; vID<vert_count(); vID++)
	{
		if(!vertex_is_valid(vID) )
			continue;
		double v = 0;
		MxVertex &pVertex = vertex(vID);

		switch (type)
		{
		case 1:			// color by msdm
			v = pVertex.MSDM2_Local_sum;
			break;
		case 2:			// color by max curvature
			v = pVertex.KmaxCurv[0];
			break;
		case 3:			// color by importance ( ei. mesh saliency)
			v = pVertex.view_importance;
			break;
		default:
			break;
		}
		R = (v - minV) * len * 255;

		if(R>255)
			R=255;
		indiceLut=floor(R);/*
		v_infor_colors[vID][0] = LUT_CourbureClust[3*indiceLut];
		v_infor_colors[vID][1] = LUT_CourbureClust[3*indiceLut+1];
		v_infor_colors[vID][2] = LUT_CourbureClust[3*indiceLut+2];*/
		pVertex.colorIndex = indiceLut;
	}

	// 计算面片的颜色级别，面片importance不一定和vertex一致，但暂时借用vertex的max 和min  [7/25/2012 Han]
	// face的颜色使用face的最大最小值来计算 [3/7/2013 Han]
	len = 1.0/(MaxImportanceF - MinImportanceF);
	for (MxFaceID f = 0; f < face_count(); f++)
	{
		double impor = 0;
		MxFace &pFace = face(f);

		switch (type)
		{
		case 1:			// color by msdm
			//v = pVertex.MSDM2_Local_sum;
			break;
		case 2:			// color by max curvature
			//v = pVertex.KmaxCurv[0];
			break;
		case 3:			// color by importance ( ei. mesh saliency)
			impor = pFace.view_importance;
			break;
		default:
			break;
		}
		R = (impor - MinImportanceF) * len * 255;

		if(R>255)
			R=255;
		indiceLut=floor(R);/*
		v_infor_colors[vID][0] = LUT_CourbureClust[3*indiceLut];
		v_infor_colors[vID][1] = LUT_CourbureClust[3*indiceLut+1];
		v_infor_colors[vID][2] = LUT_CourbureClust[3*indiceLut+2];*/
		pFace.colorIndex = indiceLut;
	}
}

void MxStdModel::UpdateAfterSimp(MxVertexID toVID, MxVertexID fromVID)
{
	double gMsdm = 0.0;
	MSDM2::UpdateMSDM2(m_pSrcMesh, this,fromVID, toVID, gMsdm, false, false);
}


bool MxStdModel::IsOccluded(float p[3], GLfloat *bufferZ)
{
	GLdouble errorAdjust = 0.0005;		// 用于调整深度缓存误差带来的问题，比如相邻顶点占据同一个像素带来的误差
	//GLint viewport[4];											//space for viewport data
	//GLdouble mvmatrix[16], projmatrix[16];  //space for transform matricex
	GLdouble winx, winy, winz;							//space for returned projected coords
	GLdouble flareZ;												//here we will store the transformed flare Z

	// Now we will ask OGL to project some geometry for us using the gluProject function.
	// Practically we ask OGL to guess where a point in space will be projected in our current viewport,
	// using arbitrary viewport and transform matrices we pass to the function.
	// If we pass to the function the current matrices  (retrievede with the glGet funcs)
	// we will have the real position on screen where the dot will be drawn.
	// The interesting part is that we also get a Z value back, this means that 
	// reading the REAL buffer for Z values we can discover if the flare is in front or
	// if it's occluded by some objects.

	//glGetIntegerv (GL_VIEWPORT, viewport);						//get actual viewport
	//glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);			//get actual model view matrix
	//glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);	//get actual projiection matrix

	// this asks OGL to guess the 2d position of a 3d point inside the viewport
	gluProject(p[0], p[1], p[2], mvmatrix, projmatrix, viewport, &winx, &winy, &winz);

	if (winx < 0 || winx > viewport[2]-1 || winy < 0 || winy > viewport[3]-1)
		return true;

	flareZ = winz;

	// test [5/7/2012 Han]
	flareZ -= errorAdjust;

	// we read back one pixel from th depth buffer (exactly where our flare should be drawn)
	//glReadPixels(winx, winy,1,1,GL_DEPTH_COMPONENT, GL_FLOAT, &bufferZ);

	// if the buffer Z is lower than our flare guessed Z then don't draw 
	// this means there is something in front of our flare
	//int l, r, u, d;
	//l = int(winx);
	//l = l < 0 ? 0:l;
	//r = int(winx+1);
	//r = r > viewport[2] - 1 ? viewport[2] - 1 : r;
	//u = int(winy+1);
	//u = u > viewport[3] - 1 ? viewport[3] - 1 : u;
	//d = int(winy);
	//d = d < 0 ? 0 : d;
	//if (  bufferZ[u*viewport[2] + l] < flareZ
	//	&&bufferZ[u*viewport[2] + r] < flareZ
	//	&&bufferZ[d*viewport[2] + l] < flareZ
	//	&&bufferZ[d*viewport[2] + r] < flareZ)
	//	return true;
	//else
	//	return false;

	if (bufferZ[int(winy)*viewport[2] + int(winx)] < flareZ)
		return true;
	else
		return false;
}
//////////////////////////////////////////////////////////////////////////
// type =	1:saliency importance
//			2:saliency importance
//			3:view entropy using mean curvature
//			4:view entropy using saliency
//			5:segment using saliency [7/10/2012 Han]
//			6:view entropy using segmentation of saliency
double MxStdModel::GetVisibleImportance(int &visibleN, ViewSelectType type, bool bVisibleFace)
{
	double importance = 0.0;
	visibleN = 0;
	bool *pVis;

	vec3 eye_normal = eye_pos;
	eye_normal.normalize();				// 归一化视点位置

	switch(type)
	{
	case SALIENCY:
	case SALIENCY_SMOOTHED:
	case MEAN_CURVATURE:
		if (bVisibleFace)
		{
			//////////////////////////////////////////////////////////////////////////
			// 利用三角形三个顶点以及中点的深度和当前像素深度比较来决定三角形是否可见
			pVis = GetVisibleFaces(visibleN);

			for (unsigned int i=0; i<face_count(); i++)
				if (pVis[i])
					importance += face(i).view_importance;
		}
		// 利用顶点是否可见来计算重要度 [6/8/2012 Han]
		else
		{
			pVis = GetVisibleVerts(visibleN);

			for (unsigned int i=0; i<vert_count(); i++)
			{
				if (pVis[i])
				{
					importance += vertex(i).view_importance ;
				}
			}
		}	
	break;

	case SEMANTIC_DRIVEN:
		if (bVisibleFace)
			pVis = GetVisibleFaces(visibleN);
		// 利用顶点是否可见来计算重要度 [6/8/2012 Han]
		else
			pVis = GetVisibleVerts(visibleN);

		if (InitSegData(bVisibleFace))
			importance = SegmentImportance(visibleN, pVis, bVisibleFace);
		break;
	case MAX_AREA:
		pVis = GetVisibleFaces(visibleN);
		for (unsigned int i=0; i<face_count(); i++)
		{
			if (pVis[i])
			{
				float area = compute_face_area(i);
				// 计算投影到视平面的面积二维三角形的实际面积 [8/18/2012 Han]
				float fn[3];
				memcpy(fn, face(i).normal, 3*sizeof(float));
				vec3 v = fn;
				v.normalize();
				area =abs(area * v.dot(eye_normal));
				importance += area;
			}
		}
		break;
	case BASE_CURVATURE:
//////////////////////////////////////////////////////////////////////////
// Wang
//	1）对包围球面上采样的每个视点，只处理模型中离视点近的部分面片。
//		比如说，所有模型面片对于这个视点而言的深度值范围为depth，则处理的面片的深度值是位于前面depth/5的。
//		对它们用平均曲率的视点检测方法进行检测。这样，那些近处有较大平坦区域的视点，就作为候选视点。（也可加上你所称的镜面对称点）。
//	2）对候选视点，用你说的支撑性和稳定性进行进一步的检测（针对整个模型），筛选出一部分视点。
//	3）对剩下的视点，进行整个模型的视点选择时的度量（依然用平均曲率的方法），将得分最差的作为摆正方向。
//	以上是用的depth/5，也可depth/4, depth/6, 等等。

		// 使用基于图像的方法 [3/23/2014 Han]
		GLint viewport[4];											//space for viewport data
		glGetIntegerv (GL_VIEWPORT, viewport);						//get actual viewport
		int	w = viewport[2], h = viewport[3];

		unsigned char *curvData = new unsigned char [w*h];
		unsigned char *depthData = new unsigned char [w*h];
		int zMax = 0, zMin = 255;	// 找到图像中的最大和最小深度值
		int area = 0;
		//ust change the packing to ensure no overruns!
		glPixelStorei(GL_PACK_ALIGNMENT, 1);

		// OpenGL以左下角为坐标原点，而windows系统是左上角
		glReadPixels(0, 0, w, h, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, depthData);

		glReadPixels(0, 0, w, h, GL_RED, GL_UNSIGNED_BYTE, curvData);

		// 找到中心点坐标和图像面积
		for(int i=0;i<w;i++)
		{
			for(int j=0;j<h;j++)
			{
				int p = j*w + i;
				if(depthData[p]!=255) // 该像素非背景
				{ 
					if (depthData[p] < zMin)
						zMin = depthData[p];
					if (depthData[p] > zMax)
						zMax = depthData[p];
				} 
			}
		}
		float zThreshold = zMin + (zMax - zMin) * 0.1;//play with this number
		for(int i=0;i<w;i++)
		{
			for(int j=0;j<h;j++)
			{
				int p = j*w + i;
				if(depthData[p]!=255 && depthData[p] < zThreshold) //该像素非背景,并且离视点足够近
				{
					importance += curvData[p];
					area++;
				}
				else
					curvData[p] = 255;
			}
		}
		importance /= area;		// 是否计算平均曲率

		// output rendered images [3/23/2014 Han]
		//PPM *ppmWriter = new PPM;
		//ppmWriter->width = w;
		//ppmWriter->height = h;
		//ppmWriter->version = "P6";
		//ppmWriter->data = new unsigned char [w*h*3];;
		//for (int xxx = 0; xxx < w*h*3; xxx+=3)
		//{
		//	ppmWriter->data[xxx] = ppmWriter->data[xxx+1] = ppmWriter->data[xxx+2] = curvData[xxx/3];
		//}
		//wchar_t fnRe[512];
		//swprintf( fnRe,   L"OutPutFiles\\ViewCurvature_%.2lf_average_%.2lf_%d.ppm",importance*area,importance, rand());
		//ppmWriter->save(fnRe);
		//delete ppmWriter;

		delete []curvData;
		delete []depthData;
		break;
		//default:
		//	return 0.0;
	}

	if (visibleN > 0)
		delete []pVis;

	// test using mean importance
	//if (visibleN != 0)
	//	return importance / visibleN;
	//else
	//	return importance;
	// test using visible primitave
	//return visibleN;
	return importance;
}

double MxStdModel::GetVisibleImportance(int &visibleN, ViewSelectType type, bool bVisibleFace, float &E, float &current_standard_deviation)
{
	double importance = 0.0;
	visibleN = 0;
	bool *pVis;
	//////////////////////////////////////////////////////////////////////////
	// 利用三角形三个顶点以及中点的深度和当前像素深度比较来决定三角形是否可见
	if (bVisibleFace || type == VIEW_ENTROPY)
		pVis = GetVisibleFaces(visibleN);
	// 利用顶点是否可见来计算重要度 [6/8/2012 Han]
	else
		pVis = GetVisibleVerts(visibleN);
	vec3 eye_normal = eye_pos;
	eye_normal.normalize();				// 归一化视点位置
	std::vector<float> areaVis;
	areaVis.resize(visibleN, 0.0f);
	float total_area = 0.f;
	int nVis = 0;

	switch(type)
	{
	case SALIENCY_ENTROPY:
	case CURVATURE_ENTROPY:
		importance = compute_shannon_entropy_II(visibleN, pVis, type, bVisibleFace, E, current_standard_deviation);
		break;

	case SEGMENT_SEMANTIC_ENTROPY:
		if (InitSegData(bVisibleFace))
			importance = Segment_shannon_entropy_II(visibleN, pVis, bVisibleFace,E,current_standard_deviation);
		break;

	case MEAN_CURVATURE_ENTROPY:
		importance = compute_maxvisible_shannon(visibleN, pVis,bVisibleFace, E, current_standard_deviation);
		break;

	case VIEW_ENTROPY:
		for (unsigned int i=0; i<face_count(); i++)
		{
			if (pVis[i])
			{
				float area = compute_face_area(i);
				// 计算投影到视平面的面积二维三角形的实际面积 [8/18/2012 Han]
				float fn[3];
				memcpy(fn, face(i).normal, 3*sizeof(float));
				vec3 v = fn;
				v.normalize();
				area =abs(area * v.dot(eye_normal));
				total_area += area;	// 投影面积和重要度乘积作为重要度，这样考虑了三角形的大小区别
				areaVis[nVis++] = area;
			}
		}
		float pi;
		E = 0.f;
		for (vector<float>::size_type i=0; i<areaVis.size(); i++)
		{
			if (areaVis[i] != 0)
			{
				pi = areaVis[i]/total_area;
				E += -pi*log(pi)/*/log(2.0f)*/;
			}
		}
		importance = E;		
		break;

	default:
		if (visibleN > 0)
			delete []pVis;
		return 0.0;
	}

	if (visibleN > 0)
		delete []pVis;

	//if (visibleN != 0)
	//	return importance / visibleN;
	//else
	//	return importance;
	// test using visible primitave
	//return visibleN;
	return importance;
}

float MxStdModel::compute_shannon_entropy_II(int &visibleN, bool *pVis, int type, bool	bVisibleFace, float &E, float &current_standard_deviation)
{
	//curvature_counted_low_bound = 0;
	//curvature_counted_up_bound = 255;

	//mean_curvature_amplified_times = 1000;
	//current_standard_deviation = 0;
	//adaptive_box_size = 1;
	//int number_of_histogram_intervals = 512;


	// curvatures for computing entropy
	std::vector<float> curvatures;
	curvatures.resize(visibleN, 0.0f);
	int visV = 0;
	// use 'as' for abbreviation
	//int as = adaptive_box_size;
	//int as = 1;
	// statistic of curvature
	if (bVisibleFace)
	{
		for (unsigned int i=0; i<face_count(); i++)
		{
			if (pVis[i])
			{
				// 统一都设置到view_importance中 [7/12/2012 Han]
				//for (int m = 0; m < 3; m++)
				//	curvatures[visV] += (type == 3)?vertex(face(i).v[m]).KmaxCurv[0]:vertex(face(i).v[m]).view_importance;
				curvatures[visV] = face(i).view_importance;

				visV += 1;
			}
		}
	}
	else
	{
		for (unsigned int i=0; i<vert_count(); i++)
		{
			if (pVis[i])
			{
				//curvatures[visV] = (type == 3)?vertex(i).KmaxCurv[0]:vertex(i).view_importance;
				curvatures[visV] = vertex(i).view_importance;
				visV += 1;
			}
		}
	}


	float maxc=0,minc=0;
	for (unsigned int i=0; i<curvatures.size(); i++)
	{
		if (curvatures[i] > maxc)
		{
			maxc = curvatures[i];
		}
		if (curvatures[i] < minc)
		{
			minc = curvatures[i];
		}
	}

	// 使用默认的32阶直方图方式 [10/12/2013 Han]
	//const unsigned int iHISTO_WIDTH = 512;
	std::vector<unsigned int> histogram;
	float hs = (maxc - minc)/HISTO_WIDTH; //histogram step width
	histogram.resize(HISTO_WIDTH,0);
	for (unsigned int i=0; i<curvatures.size(); i++)
	{
		unsigned int idx = int((curvatures[i] - minc)/hs);
		idx = idx >= HISTO_WIDTH? (HISTO_WIDTH-1):idx;
		if (idx < histogram.size())
		{
			histogram[idx] += 1;
		}
	}

	// computing standard deviation
	float avg_mc=0; // average mean curvature,
	for (vector<int>::size_type i=0; i<histogram.size(); i++)
	{
		avg_mc += histogram[i];
	}
	avg_mc/=histogram.size();
	float xigma_mc = 0;
	for (vector<int>::size_type i=0; i<histogram.size(); i++)
	{
		xigma_mc += (avg_mc - histogram[i])*(avg_mc - histogram[i]);
	}
	xigma_mc /= histogram.size() - 1;
	current_standard_deviation = sqrt(xigma_mc);

	// computing entropy 
	float N=0;
	for (vector<int>::size_type i=0; i<histogram.size(); i++)
	{
		N += histogram[i];
	}
	E=0;
	float pi;
	for (vector<int>::size_type i=0; i<histogram.size(); i++)
	{
		if (histogram[i] != 0)
		{
			pi = histogram[i]/N;
			E += -pi*log(pi)/log(2.0f);
		}
	}

	// LOG("Shannon Entropy II: "+QString::number(E));
	// LOG("Standard Deviation II: "+QString::number(current_standard_deviation));
	return E;
}
//
float MxStdModel::visible_tri_revised_entropy_II(int &visibleN, bool *pVis, int type, bool bVisibleFace, float &E, float &current_standard_deviation)
{
	return compute_shannon_entropy_II(visibleN, pVis, type, bVisibleFace,E, current_standard_deviation);
}
//{
//	//time cost
//	clock_t start_t = clock();
//
//	for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
//	{
//		//set auto_sampling_veiwpoint and compute the entropy
//		eye_pos_to_bsphere_center = viewpoint_candidates[i]->pos;
//		paintGL();
//		shannon_entropy[i] = compute_shannon_entropy_II();
//		standard_deviations[i] = current_standard_deviation; //computed in compute_shannon_entropy_II()
//	}
//
//	// get mean std_deviation
//	current_avg_stand_deviation=0;
//
//	for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
//	{
//		current_avg_stand_deviation += standard_deviations[i];
//	}
//	current_avg_stand_deviation/=VIEWPOINT_SAMPLES_NUM;
//
//	// get mean shannon entropy
//	current_avg_shannon_entropy=0;
//	for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
//	{
//		current_avg_shannon_entropy += shannon_entropy[i];
//	}
//	current_avg_shannon_entropy /= VIEWPOINT_SAMPLES_NUM;
//
//	// get rebalanced E
//	for (int i=0; i<VIEWPOINT_SAMPLES_NUM; i++)
//	{
//		standard_deviations[i] -= current_avg_stand_deviation;
//		standard_deviations[i] = abs(standard_deviations[i]);
//		revised_entropy[i] -= 3*standard_deviations[i]/(current_avg_stand_deviation/current_avg_shannon_entropy);
//	}
//
//	// find max revised entropy
//	float max_E = revised_entropy[0];
//	int max_E_pos = 0;
//	for (int i=1; i<VIEWPOINT_SAMPLES_NUM; i++)
//	{
//		if (max_E < revised_entropy[i])
//		{
//			max_E = revised_entropy[i];
//			max_E_pos = i;
//		}                                    
//	}
//
//	clock_t finish_t = clock();
//	secs = double(finish_t - start_t) / CLOCKS_PER_SEC;
//	int test = 1;// LOG("Total time elapsed: "+QString::number(secs)+" seconds\n");
//	eye_pos_to_bsphere_center = viewpoint_candidates[max_E_pos]->pos;
//}


inline void multMatirxVec4(float dst[4], const float srcm[16], const float srcv[4])
{
	float t[4];
	t[0] = srcv[0];
	t[1] = srcv[1];
	t[2] = srcv[2];
	t[3] = srcv[3];

	dst[0] = srcm[0]*t[0] + srcm[1]*t[1] + srcm[2]*t[2] + srcm[3]*t[3];
	dst[1] = srcm[4]*t[0] + srcm[5]*t[1] + srcm[6]*t[2] + srcm[7]*t[3];
	dst[2] = srcm[8]*t[0] + srcm[9]*t[1] + srcm[10]*t[2] + srcm[11]*t[3];
	dst[3] = srcm[12]*t[0] + srcm[13]*t[1] + srcm[14]*t[2] + srcm[15]*t[3];
}

inline void multMatirxVec3(float dst[3], const float srcm[16], const float srcv[3])
{
	float t[3];
	t[0] = srcv[0];
	t[1] = srcv[1];
	t[2] = srcv[2];

	dst[0] = srcm[0] *t[0] + srcm[1] *t[1] + srcm[2] *t[2] + srcm[3];
	dst[1] = srcm[4] *t[0] + srcm[5] *t[1] + srcm[6] *t[2] + srcm[7];
	dst[2] = srcm[8] *t[0] + srcm[9] *t[1] + srcm[10]*t[2] + srcm[11];
}

inline void transposeM(float m[16])
{
	float tmp[16];
	int i,j;
	for(i=0; i<16; i++)
		tmp[i] = m[i];

	for(i=0; i<4; i++)
		for(j=0; j<4; j++)
			m[i*4+j] = tmp[i+j*4];
}

typedef struct  
{
	MxStdModel* thisMesh;
	GLfloat *zBuffer;
	GLfloat *srcm;
	int idx;
	bool *visibleF;
	int visibleFN;
	int num;			// calc num of vertex [6/5/2012 Han]

} OccludDate;
// 利用多线程来加速计算,直接计算num个顶点的saliency
/*DWORD WINAPI*/unsigned __stdcall ThreadGetVisibleFaces(LPVOID param)
{
	OccludDate *sdP = (OccludDate*) param;

	float center_T[3];
	for (unsigned int i=sdP->idx; i<(sdP->idx+sdP->num); i++)
	{
		sdP->visibleF[i] = false;
		if (!sdP->thisMesh->face_is_valid(i))
			continue;

		// 处理背向面 [8/17/2012 Han]
		float fn[3], fnt[3];
		if(sdP->thisMesh->bFnormal)
			memcpy(fn, sdP->thisMesh->face(i).normal, 3*sizeof(float));
		else
		{
			sdP->thisMesh->compute_face_normal(i, fn);
			memcpy(sdP->thisMesh->face(i).normal, fn, 3*sizeof(float));
		}
		if (sdP->thisMesh->bRotated)
		{
			multMatirxVec3(fnt, sdP->srcm, fn);
			vec3 v = fnt;
			if(v.dot(vec3(0, 0, -1)) >0)
				continue;
		}
		else 
		{
			vec3 v = fn;
			if(v.dot(vec3(-sdP->thisMesh->eye_pos[0], -sdP->thisMesh->eye_pos[1], -sdP->thisMesh->eye_pos[2])) >0)
				continue;
		}

		// 三角形中点
		for(int n = 0 ; n < 3; n++)
		{
			center_T[n] = 0;
			for (int m = 0; m < 3; m++)
			{
				center_T[n] += 	sdP->thisMesh->vertex(sdP->thisMesh->face(i).v[m])[n];
			}
			center_T[n] /= 3.0f;
		}
		if (!sdP->thisMesh->IsOccluded(center_T, sdP->zBuffer)
			||!sdP->thisMesh->IsOccluded(sdP->thisMesh->vertex(sdP->thisMesh->face(i).v[0]), sdP->zBuffer)
			||!sdP->thisMesh->IsOccluded(sdP->thisMesh->vertex(sdP->thisMesh->face(i).v[1]), sdP->zBuffer)
			||!sdP->thisMesh->IsOccluded(sdP->thisMesh->vertex(sdP->thisMesh->face(i).v[2]), sdP->zBuffer))
			/*!IsOccluded(center_T, zBuffer)||(!IsOccluded(vertex(face(i).v[0]), zBuffer)&&!IsOccluded(vertex(face(i).v[1]), zBuffer)
			&&!IsOccluded(vertex(face(i).v[2]), zBuffer))*/
		{
			sdP->visibleFN++;
			sdP->visibleF[i] = true;
		}
	}

	return 0;
}
bool* MxStdModel::GetVisibleFaces(int &visibleFN, bool bMultThread)
{
	visibleFN = 0;
	//////////////////////////////////////////////////////////////////////////
	// 利用三角形三个顶点以及中点的深度和当前像素深度比较来决定三角形是否可见
	glGetIntegerv (GL_VIEWPORT, viewport);						//get actual viewport
	glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);			//get actual model view matrix
	glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);	//get actual projiection matrix
	int w = viewport[2], h = viewport[3];

	GLfloat *zBuffer = new GLfloat [w*h]; 

	// we read back one pixel from th depth buffer (exactly where our flare should be drawn)
	glReadPixels(0, 0,w,h,GL_DEPTH_COMPONENT, GL_FLOAT, zBuffer);

	// 不进行背向面剔除，因为现在使用的方法会使得计算出现失误，影响最终的视点选择结果 [8/14/2012 Han]
	GLfloat srcm[16];
	glGetMatrix(srcm, GL_MODELVIEW_MATRIX);
	transposeM(srcm);

	bool *visibleF = new bool[face_count()];

	// 进行逐面片的可见性比较 [1/27/2014 Han]
	if (bMultThread)
	{
		int nCPU = GetCpuNum();
		// 修改为使用多线程方式 [3/13/2012 Han]
		int nMaxThread = nCPU * 2;
		unsigned *dwThreadId = new unsigned [nMaxThread];
		HANDLE *hThread = new HANDLE [nMaxThread]; 
		OccludDate **od = new OccludDate* [nMaxThread];
		int numF = face_count() / nMaxThread;
		// for each vertices, compute 
		for (unsigned int i=0; i<nMaxThread; i++)
		{
			od[i] = new OccludDate;
			od[i]->srcm = srcm;
			od[i]->zBuffer = zBuffer;
			od[i]->visibleFN = 0;
			od[i]->visibleF = visibleF;
			od[i]->thisMesh = this;
			od[i]->idx = i*numF;
			od[i]->num = (i==nMaxThread-1)?(face_count() - numF*(nMaxThread-1)):numF;
			hThread[i] = (HANDLE)_beginthreadex( NULL, 0, ThreadGetVisibleFaces, od[i], 0, &dwThreadId[i] );
		}
		// Close all thread handles and free memory allocation.
		for (int i = 0; i < nMaxThread; i++)
		{
			if (WAIT_FAILED == WaitForSingleObject(hThread[i],INFINITE))
				continue;
			if (od[i] != NULL)
			{
				visibleFN += od[i]->visibleFN;
				delete od[i];
				od[i] = NULL;
			}
			CloseHandle(hThread[i]);
		}
		delete []od;
		delete []dwThreadId;
		delete []hThread;
	}
	else
	{
		OccludDate *ods = new OccludDate;
		ods->srcm = srcm;
		ods->zBuffer = zBuffer;
		ods->visibleFN = 0;
		ods->visibleF = visibleF;
		ods->thisMesh = this;
		ods->idx = 0;
		ods->num = face_count();
		ThreadGetVisibleFaces(ods);
		visibleFN = ods->visibleFN;
		delete ods;
	}
	bFnormal = true;
	delete []zBuffer;
	return visibleF;
}

bool* MxStdModel::GetVisibleVerts(int &visibleVN)
{
	visibleVN = 0;

	glGetIntegerv (GL_VIEWPORT, viewport);						//get actual viewport
	glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);			//get actual model view matrix
	glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);	//get actual projiection matrix
	int w = viewport[2], h = viewport[3];

	GLfloat *zBuffer = new GLfloat [w*h]; 

	// we read back one pixel from th depth buffer (exactly where our flare should be drawn)
	glReadPixels(0, 0,w,h,GL_DEPTH_COMPONENT, GL_FLOAT, zBuffer);

	bool *visibleV = new bool[vert_count()];
	
	float center_T[3];
	for (unsigned int i=0; i<vert_count(); i++)
	{
		visibleV[i] = false;
		if (!vertex_is_valid(i))
			continue;

		if (!IsOccluded(vertex_position(i), zBuffer))
		{
			visibleVN++;
			visibleV[i] = true;
		}
	}

	delete []zBuffer;
	return visibleV;
}

bool MxStdModel::InitSegData(bool bFace, bool bAverage, bool bLowest)
{
	if (nNumSeg == 0)
	{
		return false;
	}
	if (nNumSeg > 0 && saliencySegData.size() > 0)
		return true;

	int num = bFace?face_count():vert_count();
	float area_total = 0.0f;
	float area = 0.0f;
	for (int i = 0; i < num; i++)
	{
		bool bValid = bFace?face_is_valid(i):vertex_is_valid(i);
		if (!bValid)
			continue;
		if (bFace)
		{
			area = compute_face_area(i);
			area_total += area; 
		}

		int id = bFace?face(i).from:vertex(i).from;

		// 不应该出现，舍去 [7/11/2012 Han]
		if(id < 0)
			continue;
		float saliency = bFace?face(i).view_importance:vertex(i).view_importance;
		if (saliencySegData.count(id) > 0)
		{
			saliencySegData[id].num++;
			if (bAverage)
				saliencySegData[id].Saliency += saliency;
			else
			{
				if ((bLowest && saliency < saliencySegData[id].Saliency)
					||(!bLowest && saliency > saliencySegData[id].Saliency))
				{
					saliencySegData[id].Saliency =  saliency;
					saliencySegData[id].ID = i;
					for (int ttt = 0; ttt < 3; ttt++)
						saliencySegData[id].pos[ttt] = vertex(i)[ttt];
				}				
			}			
			if (bFace)
				saliencySegData[id].areaPer += area;
			else if (bAverage)
			{
				for (int ttt = 0; ttt < 3; ttt++)
					saliencySegData[id].pos[ttt] += vertex(i)[ttt];
			}

		}
		else
		{
			struct Seg seg;
			seg.num = 1;
			seg.ID = i;
			seg.Saliency = saliency;
			saliencySegData[id] = seg;
			if (bFace)
				saliencySegData[id].areaPer = area;
			else
			{
				for (int ttt = 0; ttt < 3; ttt++)
					saliencySegData[id].pos[ttt] = vertex(i)[ttt];
			}
			
		}
	}

	// test 修改面片的最大最小重要度值
	MinImportanceF = FLT_MAX, MaxImportanceF = 0;

	map<int,struct Seg>::iterator it;
	for ( it=saliencySegData.begin() ; it != saliencySegData.end(); it++ )
	{
		if (bAverage)
			it->second.Saliency /= it->second.num;
		it->second.area = it->second.areaPer;
		it->second.areaPer /= area_total;
		// test 修改面片的最大最小重要度值 [8/14/2013 Han]
		if(it->second.Saliency > MaxImportanceF)
			MaxImportanceF = it->second.Saliency;
		if (it->second.Saliency < MinImportanceF)
			MinImportanceF = it->second.Saliency;
		if (!bFace && bAverage)
		{
			for (int ttt = 0; ttt < 3; ttt++)
				it->second.pos[ttt] /= it->second.num;
		}		
	}
	return true;
}

float MxStdModel::compute_maxvisible_shannon(int &visibleN,bool *pVis,  bool bVisibleFace, float &E, float &current_standard_deviation)
{
	// 使用对象自定义的histo级别 [12/18/2012 Han]
	//const unsigned int HISTO_WIDTH = 32;			// 32级别结果较好，play with this value
	float hs = bVisibleFace?(MaxImportanceF - MinImportanceF)/HISTO_WIDTH:(MaxImportanceV - MinImportanceV)/HISTO_WIDTH; //histogram step width
	std::vector<float> histogram_projarea;
	histogram_projarea.resize(HISTO_WIDTH,0.f);
	float maxc=0,minc=FLT_MAX;

	int visV = 0;
	vec3 eye_normal = eye_pos;
	eye_normal.normalize();				// 归一化视点位置
	float total_area = 0.f;
	float proj_importance = 0.0f;
	float total_importance = 0.f;
	if (bVisibleFace)
	{
		for (unsigned int i=0; i<face_count(); i++)
		{
			if (face_is_valid(i))
				total_importance += face(i).view_importance;
			if (pVis[i])
			{
				float area = compute_face_area(i);
				// 计算投影到视平面的面积二维三角形的实际面积 [8/18/2012 Han]
				float fn[3];
				memcpy(fn, face(i).normal, 3*sizeof(float));
				vec3 v = fn;
				float fw = abs(v.dot(eye_normal));
				if (bNBestSelect && pFaceInfo != NULL)
					fw = fw > pFaceInfo[i] ? (fw-pFaceInfo[i]):0;

				proj_importance += face(i).view_importance * fw;	// 用于调整当前面片的重要度，越朝向摄像机，重要度加权越高

				area =abs(area * v.dot(eye_normal));

				//proj_importance += face(i).view_importance * area;	// 投影面积和重要度乘积作为重要度，这样考虑了三角形的大小区别

				total_area += area;

				visV += 1;

				unsigned int idx = int((face(i).view_importance - MinImportanceF)/hs);
				idx = idx >= HISTO_WIDTH? (HISTO_WIDTH-1):idx;
				histogram_projarea[idx] += area;										// 当前平均曲率范围内的投影面积增加

				if (face(i).view_importance > maxc)
				{
					maxc = face(i).view_importance;
				}
				if (face(i).view_importance < minc)
				{
					minc = face(i).view_importance;
				}
			}
		}
	}
	else
	{
		for (unsigned int i=0; i<vert_count(); i++)
		{
			if (pVis[i])
			{
				//curvatures[visV] = (type == 3)?vertex(i).KmaxCurv[0]:vertex(i).view_importance;
				proj_importance += vertex(i).view_importance;
				visV += 1;
			}
		}
	}

	// Adjust importance value [9/5/2012 Han]
	proj_importance = proj_importance / total_importance;
	//proj_importance = proj_importance / total_area;		// 用总体投影面积作为调整

	// calc entropy
	int non_empty = 0;
	float pi;
	E = 0.f;
	for (vector<float>::size_type i=0; i<histogram_projarea.size(); i++)
	{
		if (histogram_projarea[i] != 0)
		{
			pi = histogram_projarea[i]/total_area;
			E += -pi*log(pi)/*/log(2.0f)*/;
			non_empty++;
		}
	}

	// computing standard deviation
	float avg_mc=total_area / non_empty;
	float xigma_mc = 0;
	for (vector<float>::size_type i=0; i<histogram_projarea.size(); i++)
	{
		if (histogram_projarea[i] != 0)
			xigma_mc += (avg_mc - histogram_projarea[i])*(avg_mc - histogram_projarea[i]);
	}
	xigma_mc /= non_empty;
	current_standard_deviation = sqrt(xigma_mc);

	float importance = proj_importance * E/* *//* * visSeg.size()*/;
	return importance;
}

float MxStdModel::Segment_shannon_entropy_II(int &visibleN,bool *pVis,  bool bVisibleFace, float &E, float &current_standard_deviation)
{
	int visibleVN = 0;

	std::map<int, int> visSeg;
	std::map<int, float> visSegArea;	// 可见分块及其投影面积
	float total_area = 0.0f;			// 投影面积之和

	int num = bVisibleFace?face_count():vert_count();

	vec3 eye_normal = eye_pos;
	eye_normal.normalize();				// 归一化视点位置

	for (int i = 0; i < num; i++)
	{
		bool bValid = bVisibleFace?face_is_valid(i):vertex_is_valid(i);
		if (!pVis[i] || !bValid)
			continue;
		
		int segID = bVisibleFace?face(i).from:vertex(i).from;
		float area = bVisibleFace? compute_face_area(i):0.0f;
		// 计算投影到视平面的面积二维三角形的实际面积 [8/18/2012 Han]
		float fn[3];
		memcpy(fn, face(i).normal, 3*sizeof(float));
		vec3 v = fn;
		area =abs(area * v.dot(eye_normal));

		total_area += area;
		if (visSeg.count(segID) > 0)
		{
			visSeg[segID]++;
			if (bVisibleFace)
				visSegArea[segID] += area;			
		}
		else
		{
			visSeg[segID] = 1;
			if (bVisibleFace)
				visSegArea[segID] = area;			
		}
	}
	total_area = total_area ==0 ? 1:total_area;
	float importance = 0.0f;
	map<int,int>::iterator it;
	float avg_mc=0; // average mean curvature,

	float fl = bVisibleFace? MaxImportanceF-MinImportanceF:MaxImportanceV-MinImportanceV;
	fl = (fl == 0)?1:(1.0f/fl);

	for ( it=visSeg.begin() ; it != visSeg.end(); it++ )
	{
		avg_mc += float(it->second);
		// 将重要度扩所为0.1~1.1 [7/11/2012 Han]
		//importance += (0.01f+(saliencySegData[it->first].Saliency-MinImportance)*fl) * float(it->second)/ saliencySegData[it->first].num;
		importance += (0.01f+saliencySegData[it->first].Saliency) * float(it->second)/ saliencySegData[it->first].num;
	}
	//////////////////////////////////////////////////////////////////////////
	// !!! 应该将importance计算引入香农值当中，调整那些能看到更多重要信息的视点
	// 暂时没有实现 [7/12/2012 Han]
	importance *= visSeg.size();

	avg_mc/=visSeg.size();
	
	float xigma_mc = 0;
	for ( it=visSeg.begin() ; it != visSeg.end(); it++ )
	{
		xigma_mc += (avg_mc - it->second)*(avg_mc - it->second);
	}

	xigma_mc /= visSeg.size() - 1;
	current_standard_deviation = sqrt(xigma_mc);

	E=0;
	float pi;
	// computing entropy 
	//float N=0;
	//for ( it=visSeg.begin() ; it != visSeg.end(); it++ )
	//{
	//	N += it->second;
	//}
	//for ( it=visSeg.begin() ; it != visSeg.end(); it++ )
	//{
	//	if (it->second != 0)
	//	{
	//		//pi = it->second/N;
	//		// 修改为面积所占可见总面积的比例 [7/21/2012 Han]
	//		pi = visSegArea[it->first] / total_area;

	//		float eT = -pi*log(pi)/log(2.0f);
	//		// 加入saliency影响因子进去 [7/13/2012 Han]
	//		// 暂时不在此处进行entropy + semantic调整，放在外部整体candidate views中调整 [7/17/2012 Han]
	//		//eT *= (0.1f+saliencySegData[it->first].Saliency/255.0f) * float(it->second)/ saliencySegData[it->first].num;
	//		eT *= (0.01f+saliencySegData[it->first].Saliency) * float(it->second)/ saliencySegData[it->first].num;
	//		E += eT;
	//	}
	//}
	// //让返回值也等于E [7/14/2012 Han]
	//importance = E;
	//// 让可见面片等于可见的分块数 [7/26/2012 Han]
	//visibleN = visSeg.size();

	//////////////////////////////////////////////////////////////////////////
	// test delete :分别求shannon值和saliency 可见面积加权的总和，视点评分是这两者之间的成绩 [8/17/2012 Han]
	float fEntropy = 0.f, maxE = 0.f, minE = FLT_MAX;
	float fSaliency = 0.f, maxS = 0.f, minS = FLT_MAX;
	float fSemantic = 0.f;

	for ( it=visSeg.begin() ; it != visSeg.end(); it++ )
	{
		pi = visSegArea[it->first] / total_area;
		float eT = -pi*log(pi)/*/log(2.0f)*/;
		// 0 log 0 = 0
		if (pi == 0)
			eT = 0;
		fEntropy += eT;
		E += eT;

		//float sT = (/*0.01f+*/saliencySegData[it->first].Saliency) * float(it->second)/ saliencySegData[it->first].num;
		///* 按照2009 semantic中的计算公式计算得分 */
		fSemantic += saliencySegData[it->first].areaPer*(0.01f+(saliencySegData[it->first].Saliency-MinImportanceF)*fl)* visSegArea[it->first] / saliencySegData[it->first].area;
		float sT = (0.01f+(saliencySegData[it->first].Saliency-MinImportanceF)*fl)* visSegArea[it->first] / saliencySegData[it->first].area;		// 投影面积占原始面积的比例，越正向的分块，这个数值越大//float(it->second)/ saliencySegData[it->first].num ;
		fSaliency += sT;
	}
	/*添加了shannon的计算公式*/
	//importance = fEntropy * fSaliency/* * visSeg.size()*/;
	/*2009 semantic中的计算公式*/
	//importance = fSemantic * visSeg.size();
	/*我们的计算公式*/
	importance = fSaliency;
	//////////////////////////////////////////////////////////////////////////

	// LOG("Shannon Entropy II: "+QString::number(E));
	// LOG("Standard Deviation II: "+QString::number(current_standard_deviation));
	return importance;			// ！！！应该返回添加了传统importance调整的值，暂时未实现 [7/13/2012 Han]
}

//  [8/15/2012 Han]
float MxStdModel::SegmentImportance(int &visibleN,bool *pVis,  bool bVisibleFace)
{
	int visibleVN = 0;

	std::map<int, int> visSeg;
	std::map<int, float> visSegArea;
	float total_area = 0.0f;

	int num = bVisibleFace?face_count():vert_count();

	vec3 eye_normal = eye_pos;
	eye_normal.normalize();				// 归一化视点位置

	for (int i = 0; i < num; i++)
	{
		bool bValid = bVisibleFace?face_is_valid(i):vertex_is_valid(i);
		if (!pVis[i] || !bValid)
			continue;
		
		int segID = bVisibleFace?face(i).from:vertex(i).from;

		float area = bVisibleFace? compute_face_area(i):0.0f;
		// 计算投影到视平面的面积二维三角形的实际面积 [8/18/2012 Han]
		float fn[3];
		memcpy(fn, face(i).normal, 3*sizeof(float));
		vec3 v = fn;
		area =abs(area * v.dot(eye_normal));

		total_area += area;
		if (visSeg.count(segID) > 0)
		{
			visSeg[segID]++;
			if (bVisibleFace)
				visSegArea[segID] += area;			
		}
		else
		{
			visSeg[segID] = 1;
			if (bVisibleFace)
				visSegArea[segID] = area;			
		}
	}

// 严格按照文章09 semantic-driven的方法进行视点重要度计算 ：
	// 利用公式，可见分块数*每个分块可见性*每个分块面积占模型总面积的比例*每个分块重要度 [7/10/2012 Han]
	float importance = 0.0f;
	map<int,int>::iterator it;
	// test 将重要度进行归一化处理 [7/26/2012 Han]
	float fl = MaxImportanceF-MinImportanceF;
	fl = (fl == 0)?1:(1.0f/fl);
	for ( it=visSeg.begin() ; it != visSeg.end(); it++ )
	{
		importance += 
			//*面积比例*//*修改为严格按照面积的形式*/saliencySegData[it->first].areaPer * /*saliencySegData[it->first].num/float(num)	*/		
			//*/*重要度*/(0.1f+saliencySegData[it->first].Saliency/255.0f) * 	// 将重要度扩所为0.1~1.1 [7/11/2012 Han]
			///*重要度*/(0.01f+(saliencySegData[it->first].Saliency-MinImportance)*fl) *	// 重要度就按照原始重要度计算，考虑重要度为0的也要有一定贡献，所以加上0.01
			/*重要度*/(0.01f+saliencySegData[it->first].Saliency) *	// 重要度就按照原始重要度计算，考虑重要度为0的也要有一定贡献，所以加上0.01
			///*可见性*/float(it->second)/ saliencySegData[it->first].num;	// 这个分块可见的数目除以总的数目
			/*可见性*/visSegArea[it->first] / saliencySegData[it->first].area;		// 投影面积占原始面积的比例，越正向的分块，这个数值越大	
	}

	importance *= visSeg.size();											// 最后和总的可见分块相乘
	// 让可见面片等于可见的分块数 [7/26/2012 Han]
	visibleN = visSeg.size();
	return importance;
}

typedef struct  
{
	MxStdModel* thisMesh;
	int idx;
	int num;			// calc num of vertex [6/5/2012 Han]
	bool isGeod;
	double radius;
	float *originalImportance;
} ImportanceData;

// 利用多线程来加速计算
/*DWORD WINAPI*/unsigned __stdcall ThreadFilterImportance(LPVOID param)
{
	ImportanceData *sdP = (ImportanceData*) param;
	for (unsigned int i=sdP->idx; i<(sdP->idx+sdP->num); i++)
	{
		sdP->thisMesh->LaplacianFilterVertexImportance(sdP->isGeod, sdP->radius, i, sdP->originalImportance);
	}

	return 0;
}

void MxStdModel::LaplacianFilterVertexImportance(bool IsGeod,double radius, MxVertexID verID, float *origianlImportance)
{
	if(!vertex_is_valid(verID) )
		return;

	if (!IsGeod)
		// not implement [6/5/2012 Han]
		return;
	float fImportanceSum = 0.0f;

	MxVertex* pVertex = &vertex(verID);

	std::set<MxVertexID> verts ;
	vec3 O = vec3(vertex_position(verID));

	std::stack<MxVertexID> S ;

	S.push(verID) ;
	verts.insert(verID) ;
	fImportanceSum += origianlImportance[verID];
	int iter=1;
	while(!S.empty())
	{
		//Vertex* v = S.top() ;
		MxVertexID v = S.top();
		// 当前顶点 [4/18/2012 Han]
		MxVertex &vCurr = vertex(v);

		S.pop() ;

		vec3 P = vec3(vCurr.as.pos);

		MxVertexList star;
		star.reset();
		collect_vertex_star(v, star);

		for(int j = 0; j < star.length(); j++)
		{
			MxVertexID w = star(j);
			vec3 p1 = vec3(vCurr.as.pos);
			vec3 p2 = vec3(vertex_position(star(j)));

			vec3 V = (p2-p1);
			//if (v==verID || V.dot(P - O) > 0.0)
			//if (star(j) != v)
			// 没有处理过的邻接点用于计算滤波 [6/5/2012 Han]
			if(verts.find(w) == verts.end())
			{
				bool isect = sphere_clip_vector(O, radius, P, V) ;

				fImportanceSum += origianlImportance[w];

				//// delete [6/5/2012 Han]
				//v_infor_colors[w][0] = 0.0;
				//v_infor_colors[w][1] = 1.0;
				//v_infor_colors[w][2] = 0.0;

				iter++;
	
				if (!isect) {
					/*if (verts.find(w) == verts.end())*/ {
						verts.insert(w) ;
						S.push(w) ;
					}
				}
			}
		}
	}
	if (iter == 1)
		return;

	// 不计算当前顶点
	fImportanceSum -= origianlImportance[verID];
	iter -= 1;

	pVertex->view_importance = fImportanceSum / iter;
}

void MxStdModel::LaplacianFilterVertexImportance(bool IsGeod,double radius, FILE *out)
{
	float *originalImportances = new float[vert_count()];
	for (int i = 0; i < vert_count(); i++)
		originalImportances[i] = vertex(i).view_importance;

	int nCPU = GetCpuNum();
	// 修改为使用多线程方式 [3/13/2012 Han]
	int nMaxThread = nCPU * 2;
	unsigned *dwThreadId = new unsigned [nMaxThread];
	HANDLE *hThread = new HANDLE [nMaxThread]; 
	ImportanceData **sd = new ImportanceData* [nMaxThread];
	int numV = vert_count() / nMaxThread;

	double tm=get_cpu_time(); 

	// for each vertices, compute 
	for (unsigned int i=0; i<nMaxThread; i++)
	{
		sd[i] = new ImportanceData;
		sd[i]->thisMesh = this;
		sd[i]->idx = i*numV;
		sd[i]->num = (i==nMaxThread-1)?(vert_count() - numV*(nMaxThread-1)):numV;
		sd[i]->isGeod = IsGeod;
		sd[i]->radius = radius;
		sd[i]->originalImportance = originalImportances;
		hThread[i] = (HANDLE)_beginthreadex( NULL, 0, ThreadFilterImportance, sd[i], 0, &dwThreadId[i] );
	}
	// Close all thread handles and free memory allocation.
	for (int i = 0; i < nMaxThread; i++)
	{
		if (WAIT_FAILED == WaitForSingleObject(hThread[i],INFINITE))
			continue;
		if (sd[i] != NULL)
			delete sd[i];
		CloseHandle(hThread[i]);
	}

	tm=get_cpu_time() - tm;


	MinImportanceV = FLT_MAX, MaxImportanceV = 0;

	int vv = 0;
	for (int i = 0; i < vert_count(); i++)
	{
		if (!vertex_is_valid(i))
			continue;

		float s = vertex(i).view_importance;
		if (s > MaxImportanceV)
			MaxImportanceV = s;
		if (s < MinImportanceV)
			MinImportanceV = s;
		vv++;
	}

	if (out != NULL)
	{
		fprintf(out	, "Valid Vert Num:\t%d\t;Smooth Radius:\t%lf\t;Smooth Saliency Time:\t%lf\t;Max Sliency:\t%f\t, Min:\t%f\t\n"
			,vv, radius, tm, MaxImportanceV, MinImportanceV);
		// put buffer data to disk [6/8/2012 Han]
		fflush(out);
	}

	delete []sd;
	delete []dwThreadId;
	delete []hThread;

	delete []originalImportances;

	// 需要重新计算分块 [7/10/2012 Han]
	saliencySegData.clear();
	nNumSeg = 0;

	// 将得到的数值做扩缩 [7/10/2012 Han]
	IdentityVertexImportance(MinImportanceV, MaxImportanceV, 0.0f, 255.f);

	// 设置每个面片的重要度，重要度依据三个顶点重要度求平均得到 [7/10/2012 Han]
	UpdateFaceImportanceByVert();
}


void MxStdModel::IdentityVertexImportance()
{
	// 不对saliency值进行任何更改，所有更改都是临时的 [7/25/2012 Han]
	// test only [7/25/2012 Han]
	return;

	float fminGard = 0, fmaxGard = 255;
	float fmin = FLT_MAX, fmax = 0;

	int vv = 0;
	for (int i = 0; i < vert_count(); i++)
	{
		if (!vertex_is_valid(i))
			continue;

		float s = vertex(i).view_importance;
		if (s > fmax)
			fmax = s;
		if (s < fmin)
			fmin = s;
		vv++;
	}
	IdentityVertexImportance(fmin, fmax, fminGard, fmaxGard);
}

// 将每个顶点的重要度值扩缩为0~255 [7/10/2012 Han]
void MxStdModel::IdentityVertexImportance(float fmin, float fmax, float fminGard, float fmaxGard)
{
	// test only [7/25/2012 Han]
	return;

	MaxImportanceV = 0;
	MinImportanceV = FLT_MAX;

	for (int i = 0; i < vert_count(); i++)
	{
		if (!vertex_is_valid(i))
			continue;

		float vi = vertex(i).view_importance;
		vertex(i).view_importance =  fminGard + (fmaxGard - fminGard)*(vi - fmin)/(fmax - fmin);

		if (vertex(i).view_importance > MaxImportanceV)
			MaxImportanceV = vertex(i).view_importance;
		if (vertex(i).view_importance < MinImportanceV)
			MinImportanceV = vertex(i).view_importance;

	}
}

// 返回分块的数目 [7/11/2012 Han]
int MxStdModel::UpdateSaliencyBySegment(bool bFace, bool bAverage)
{
	InitSegData(bFace, bAverage);

	if (!bAverage)	// 当分块完成后，不改变观察球上视点的质量 [2/17/2014 Han]
		return saliencySegData.size();

	int num = bFace?face_count():vert_count();

	for (int i = 0; i < num; i++)
	{
		bool bValid = bFace?face_is_valid(i):vertex_is_valid(i);
		if (!bValid)
			continue;

		int id = bFace?face(i).from:vertex(i).from;

		float saliency = saliencySegData[id].Saliency;

		if (id < 0)
			saliency = 0.0f;

		if (bFace)
		{
			face(i).view_importance = saliency;
			for(MxVertexID v = 0; v < 3; v++)
				vertex(face(i)[v]).view_importance = saliency;
		}
		else
			vertex(i).view_importance = saliency;
	}
	return saliencySegData.size();
}

// 返回分块的数目 [7/11/2012 Han]
void MxStdModel::UpdateFaceImportanceByVert()
{
	MinImportanceF = FLT_MAX, MaxImportanceF = 0;

	for (MxFaceID f = 0; f < face_count(); f++)
	{
		face(f).view_importance = 0.0f;
		if (!face_is_valid(f))
			continue;

		float cfi = 0.0f;
		for(MxVertexID v = 0; v < 3; v++)
			cfi += vertex(face(f)[v]).view_importance;
		face(f).view_importance = cfi/3.0f;		
		// 同时更新最大最小importance [3/7/2013 Han]
		if(face(f).view_importance > MaxImportanceF)
			MaxImportanceF = face(f).view_importance;
		if (face(f).view_importance < MinImportanceF)
			MinImportanceF = face(f).view_importance;
	}
}

inline float dist2(const MxVertex &v1, const MxVertex &v2)
{
	float d2 = sqr(v2[0]-v1[0]);
	for (int i = 1; i < 3; i++)
		d2 += sqr(v2[i]-v1[i]);
	return d2;
}
inline float dist(const MxVertex &v1, const MxVertex &v2)
{
	return sqrt(dist2(v1,v2));
}

inline float MxStdModel::GaussianWeigtedAverage(unsigned int i, float delta)
{
	//  [7/14/2012 Han]
	vertex(i).view_importance = 0.0f;
	if (!vertex_is_valid(i))
		return 0.0f;

	// all variables are the same with the equation in the paper.
	std::set<std::vector<int>::size_type> N;
	std::set<std::vector<int>::size_type> pN; // potential N;
	std::set<std::vector<int>::size_type> dN; // discarded N;

	pN.insert(i);

	// compute N
	while (!pN.empty())
	{
		std::vector<int>::size_type v = *(pN.begin());
		pN.erase(v);
		if (dist(vertex(i), vertex(v)) <= 2*delta)
		{
			N.insert(v);
		}
		else
		{
			dN.insert(v);
			continue; // if the vertex is out of delta range, then his neighbor is also discarded.
		}

		MxVertexList star;
		star.reset();
		collect_vertex_star(v, star);
		for(uint j=0; j<star.length(); j++)
			//for (std::set<unsigned int>::iterator iter = adjacent[v].begin(); iter != adjacent[v].end(); ++iter)   
		{
			//if (N.count(*iter) == 0 && dN.count(*iter) == 0)
			//{
			//	pN.insert(*iter);
			//}
			if (N.count(star(j)) == 0 && dN.count(star(j)) == 0)
			{
				pN.insert(star(j));
			}

		}
	}

	// compute the result
	double numerator=0;
	double denominator=0;
	for (std::set<unsigned int>::iterator iter = N.begin(); iter != N.end(); ++iter)
	{
		std::vector<unsigned int>::size_type idx = *iter;
		//numerator += (curv1[idx] + curv2[idx])
		//	*exp(-dist2(vertices[i],vertices[idx])/(2*delta*delta));
		//denominator += exp(-dist2(vertices[i],vertices[idx])/(2*delta*delta));
		numerator += vertex(idx).KmaxCurv[0]*exp(-dist2(vertex(i),vertex(idx))/(2*delta*delta));
		denominator += exp(-dist2(vertex(i), vertex(idx))/(2*delta*delta));

	}
	return float(numerator/denominator);
}

typedef struct  
{
	MxStdModel* thisMesh;
	float *delta;
	int idx;
	std::vector<float> *saliencyScale;
	//std::set<std::vector<int>::size_type> *N;
	//std::set<std::vector<int>::size_type> *pN; // potential N;
	//std::set<std::vector<int>::size_type> *dN; // discarded N;

	int num;			// calc num of vertex [6/5/2012 Han]

} SaliencyData;


// 利用多线程来加速计算,直接计算num个顶点的saliency
/*DWORD WINAPI*/unsigned __stdcall ThreadNewComputeSaliencys(LPVOID param)
{
	SaliencyData *sdP = (SaliencyData*) param;

	for (unsigned int i=sdP->idx; i<(sdP->idx+sdP->num); i++)
	{
		// 如果此顶点已被删除，则直接下一个顶点 [7/15/2012 Han]
		sdP->thisMesh->vertex(i).view_importance = 0.0f;
		if (!sdP->thisMesh->vertex_is_valid(i))
			continue;

		float S[5];
		for (unsigned int j=0; j<5; j++)
		{
			S[j] = abs(sdP->thisMesh->GaussianWeigtedAverage(i,sdP->delta[j]) - sdP->thisMesh->GaussianWeigtedAverage(i,2*sdP->delta[j]));
			// 分别记录每个scale的saliency值 [3/14/2012 Han]
			sdP->saliencyScale[j][i] = S[j];
		}
		// then, combine using a linear interpolation
		//saliency[i] = (S[0]+S[1]+S[2]+S[3]+S[4])/5.0f;
	}

	return 0;
}

void MxStdModel::OutputImportances(FILE* outputSaliencyFile)
{
	if (outputSaliencyFile == NULL)
		return;
	float fmin = FLT_MAX, fmax = 0;

	fprintf(outputSaliencyFile, "%d\n",face_count());
	for (int i = 0; i < vert_count(); i++)
	{
		float ss = vertex(i).view_importance;
		if (ss > fmax)
			fmax = ss;
		if (ss < fmin)
			fmin = ss;

		fprintf(outputSaliencyFile, "%d\t%f\n", vertex_is_valid(i), ss);
	}

	fprintf(outputSaliencyFile, "max saliency:\t%f\t;min saliency:\t%f\t\n", fmax, fmin);

	// test 放在此处不合适，但暂时这样做 [7/15/2012 Han]
	IdentityVertexImportance(fmin, fmax);

}

void MxStdModel::compute_mesh_saliency(float model_radius, FILE *output, bool bMultiThread)
{
    if (vert_count()<=0)
    {
        return;
    }
    float epsilon = model_radius*2;
    // according to the paper , 0.3% of the diagonal of the BBox, here use diameter to approximate.
    epsilon *= 0.003f; 

	double t1, t2, t3;
    // compute curvature;
    //TIMING(t1, compute_curvatures());
    //TIMING(t2, compute_bsphere());// bounding box
    //TIMING(t3, compute_adjacent_vertices());
	// curvature [7/14/2012 Han]
	float size[3];
	GetBoundSize(size);
	float max = 0;
	for (int i = 0; i < 3; i ++)
		max = max < size[i]? size[i] : max;

	double RadiusCurvature=0.001/*MSDM2::mini_radius*/;
	// 采用的geo半径比较关键，应该采取比较小的半径 [7/15/2012 Han]
	TIMING(t1, principal_curvature(true,RadiusCurvature*max/*epsilon*//*使用meshsaliency中最小的长度计算*/, 0));
	MSDM2::KmaxKmean(this,max, 0);

	std::vector<float> saliencyScale[5];
	for (int i = 0; i < 5; i++)
		saliencyScale[i].resize(vert_count(),0);

    // 5 is the number of scales
    // float delta[5] = {1*epsilon, 2*epsilon, 3*epsilon, 4*epsilon, 5*epsilon};
	float delta[5] = {2*epsilon, 3*epsilon, 4*epsilon, 5*epsilon, 6*epsilon};
	if (bMultiThread)
//////////////////////////////////////////////////////////////////////////
// new multi thread method [6/5/2012 Han]
	{
		int nCPU = GetCpuNum();
		// 修改为使用多线程方式 [3/13/2012 Han]
		int nMaxThread = nCPU * 2;
		unsigned *dwThreadId = new unsigned [nMaxThread];
		HANDLE *hThread = new HANDLE [nMaxThread]; 
		SaliencyData **sd = new SaliencyData* [nMaxThread];
		int numV = vert_count() / nMaxThread;
		// for each vertices, compute 
		for (unsigned int i=0; i<nMaxThread; i++)
		{
			sd[i] = new SaliencyData;
			sd[i]->thisMesh = this;
			sd[i]->delta = delta;
			sd[i]->idx = i*numV;
			sd[i]->saliencyScale = saliencyScale;
			sd[i]->num = (i==nMaxThread-1)?(vert_count() - numV*(nMaxThread-1)):numV;
			hThread[i] = (HANDLE)_beginthreadex( NULL, 0, ThreadNewComputeSaliencys, sd[i], 0, &dwThreadId[i] );
		}
		// Close all thread handles and free memory allocation.
		for (int i = 0; i < nMaxThread; i++)
		{
			if (WAIT_FAILED == WaitForSingleObject(hThread[i],INFINITE))
				continue;
			if (sd[i] != NULL)
				delete sd[i];
			CloseHandle(hThread[i]);
		}
		delete []sd;
		delete []dwThreadId;
		delete []hThread;
	}
	else
	{
		for (unsigned int i=0; i<vert_count(); i++)
		{
			float S[5];
			for (unsigned int j=0; j<5; j++)
			{
				S[j] = abs(GaussianWeigtedAverage(i,delta[j]) - GaussianWeigtedAverage(i,2*delta[j]));
				// 分别记录每个scale的saliency值 [3/14/2012 Han]
				saliencyScale[j][i] = S[j];
			}
			// then, combine using a linear interpolation
			//saliency[i] = (S[0]+S[1]+S[2]+S[3]+S[4])/5.0f;
		}
	}

// test：计算saliency好像有误，原文中是利用高斯加权计算每个半径级别的saliency Li，然后每个级别的Li乘上最大值和平均值之差的平方 [3/14/2012 Han]
	float maxsS[5], minsS[5], meanS[5];
	memset(maxsS,0,5*sizeof(float));
	memset(minsS,1024, 5*sizeof(float));
	memset(meanS,0,5*sizeof(float));
	for ( int s = 0; s < 5; s++)
	{
		for (unsigned int i = 0; i < vert_count(); i++)
		{
			if (saliencyScale[s][i] > maxsS[s])
				maxsS[s] = saliencyScale[s][i];
			if (saliencyScale[s][i] < minsS[s])
				minsS[s] = saliencyScale[s][i];
			meanS[s] += saliencyScale[s][i];
		}
		meanS[s] /= vert_count();
	}
	
	for ( int s = 0; s < 5; s++)
	{
		meanS[s] = (meanS[s] - minsS[s]) / (maxsS[s] - minsS[s]);
		for (unsigned int i = 0; i < vert_count(); i++)
		{
			saliencyScale[s][i] = (saliencyScale[s][i] - minsS[s]) / (maxsS[s] - minsS[s]);
			// non-linear normalization [3/16/2012 Han]
			saliencyScale[s][i] *= (1-meanS[s])*(1-meanS[s]);
			// result saliency = adding all scales after non-linear normalization of suppression  [3/16/2012 Han]
			//saliency[i] += saliencyScale[s][i];
			vertex(i).view_importance += saliencyScale[s][i];
		}
	}

	if (output != NULL)
		OutputImportances(output);

	UpdateFaceImportanceByVert();
}

// 按照已经选定的所有最优视点的可见面片计算其占据模型面积的比例计算香农熵 [11/1/2012 Han]
// The relative entropy or Kullback-Leibler distance between two probability distributions p = {p(x)}
// and q = {q(x)} defined over X is given by
//	KL(p | q) = Σp(x) log(p(x)/q(x))
void MxStdModel::UpdateNBestInfo(float &E, float &klD, float &visArea, float meshArea, bool bProj)
{
	float hs = (MaxImportanceF - MinImportanceF)/HISTO_WIDTH; //histogram step width
	std::vector<float> histogram_projarea;
	histogram_projarea.resize(HISTO_WIDTH,0.f);
	std::vector<float> histogram_OrigMeshArea;
	histogram_OrigMeshArea.resize(HISTO_WIDTH,0.f);

	double importance = 0.0;
	int visibleN = 0;
	bool *pVis = GetVisibleFaces(visibleN);
	vec3 eye_normal = eye_pos;
	eye_normal.normalize();				// 归一化视点位置

	visArea = 0.f;

	if (visibleN > 0)
	{
		for (unsigned int i=0; i<face_count(); i++)
		{
			if (pVis[i])
			{
				float fn[3];
				memcpy(fn, face(i).normal, 3*sizeof(float));
				vec3 v = fn;
				float fw = abs(v.dot(eye_normal));	// 用于调整当前面片的重要度，越朝向摄像机，重要度加权越高
				if (fw > pFaceInfo[i])
					pFaceInfo[i] = fw;				// 更新可见面片的权重，取最大值
			}
			float area = compute_face_area(i);
			unsigned int idx = int((face(i).view_importance - MinImportanceF)/hs);
			idx = idx >= HISTO_WIDTH? (HISTO_WIDTH-1):idx;
			// 将面积计入原始模型面积统计当中，即q [12/29/2012 Han]
			histogram_OrigMeshArea[idx] += area;

			// 如果当前面片在N best view中可见，则计入面积统计 [12/29/2012 Han]
			if (pFaceInfo[i] > 0.f)
			{
				area = bProj ? area * pFaceInfo[i] : area;
				visArea += area;
				// 将面积计入投影面积统计当中，即p [12/29/2012 Han]
				histogram_projarea[idx] += area;										// 当前平均曲率范围内的投影面积增加
			}
		}
	}
	klD = 0.f;
	E = 0.f;
	for (vector<float>::size_type i=0; i<histogram_projarea.size(); i++)
	{
		if (histogram_projarea[i] != 0)
		{
			// KL(p | q) = Σp(x) log(p(x)/q(x))
			float pi = histogram_projarea[i]/ /*meshArea*/ visArea;
			float qi = histogram_OrigMeshArea[i]/meshArea;
			klD += pi*log(pi/qi);

			pi = histogram_projarea[i]/meshArea /*visArea*/;
			E += -pi*log(pi)/*/log(2.0f)*/;
		}
	}

	if (visibleN > 0)
		delete []pVis;
}

float MxStdModel::CalcMeshShannon(float &total_area)
{
	// Use class's  [12/18/2012 Han]
	//const unsigned int HISTO_WIDTH = 32;			// 32级别结果较好，play with this value
	float hs = (MaxImportanceF - MinImportanceF)/HISTO_WIDTH; //histogram step width
	std::vector<float> histogram_projarea;
	histogram_projarea.resize(HISTO_WIDTH,0.f);
	float maxc=0,minc=FLT_MAX;

	total_area = 0.f;
	for (unsigned int i=0; i<face_count(); i++)
	{
			if (face_is_valid(i))
			{
				float area = compute_face_area(i);
				total_area += area;

				unsigned int idx = int((face(i).view_importance - MinImportanceF)/hs);
				idx = idx >= HISTO_WIDTH? (HISTO_WIDTH-1):idx;
				histogram_projarea[idx] += area;										// 当前平均曲率范围内的投影面积增加

			}
	}
	float E = 0.f;
	for (vector<float>::size_type i=0; i<histogram_projarea.size(); i++)
	{
		if (histogram_projarea[i] != 0)
		{
			float pi = histogram_projarea[i]/total_area;
			E += -pi*log(pi)/*/log(2.0f)*/;
		}
	}
	return E;
}

float round(float a) {
	a += 0.5f;
	return floor(a);
}

// 修改为统计face的直方图 [3/8/2013 Han]
void MxStdModel::Histoeq(int band, float alpha) 
{
	bool bUseArea = false/*true*/; // 2013年7月15日，虽然使用面积更加符合实际情况，但使用面片得到的结果更好
	float area, totalA = 0.f;
	float imMin = /*image_in.min()*/ MinImportanceF;
	float imMax = /*image_in.max()*/ MaxImportanceF;

	//MinImportanceF = FLT_MAX;
	//MaxImportanceF = 0;

	if(band == -1)
		band = HISTO_WIDTH;

	//float range = (imMax - imMin)/* + 1.f*/; 
	//float binSize = /*ceil*/(range / band);

	//cout << "Bin size: " << binSize << endl;
	float hs = (imMax - imMin)/band; //histogram step width
	std::vector<float> binFreq;
	binFreq.resize(band,0.f);
	for (int i = 0; i < face_count(); i++)
	{
		unsigned int idx = int((face(i).view_importance - imMin)/hs);
		idx = idx >= band? (band-1):idx;
		//face(i).view_importance = idx;										// 当前平均曲率范围内的投影面积增加
		if (bUseArea)
		{
			area = compute_face_area(i);
			binFreq[idx] += area;
			totalA += area;
		}
		else
			binFreq[idx] += 1;

	}

	float mnCount = 0.f;
	for (int i = 0; i < band; i++) {
		mnCount += binFreq[i];
		binFreq[i] = binFreq[i] / (bUseArea?totalA:face_count());
	}

	//cout << "mnCount: " << mnCount << endl;

	for (int i = 0; i < face_count(); i++){
		float pixVal = face(i).view_importance;

		// alternate logic, takes less iterations 
		int binLoc =(int)((pixVal - imMin)/ hs/*binSize*/);

		// to accomodate the max value in the range
		if (binLoc == band)
			binLoc--;

		if (binLoc > (band - 1)) {
			//cout << "Error in calculating binLoc\n";
			exit(-1);
		}	

		float cdf = 0.f;
		for (int k = 0; k <= binLoc; k++) {
			cdf += binFreq[k];
		}

		float newPixVal = /*round*/(cdf * (imMax - imMin) + imMin); 
		// float newPixVal = round(((cdf - binFreq[0]) / (mnCount - binFreq[0]))* (max - min) + min); 

		// image_out->poke(i, j) = newPixVal; 

		// blending
		face(i).view_importance = alpha * newPixVal + (1.f - alpha) * pixVal;

		//  [3/7/2013 Han]
		//if(face(i).view_importance > MaxImportanceF)
		//	MaxImportanceF = face(i).view_importance;
		//if(face(i).view_importance < MinImportanceF)
		//	MinImportanceF = face(i).view_importance;
	}  
}		     
