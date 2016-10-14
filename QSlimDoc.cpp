#include "stdafx.h"
#include "MeshSimpDoc.h"
/************************************************************************

  This file implements the command line parsing interface for QSlim.

  Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.
  
  $Id: cmdline.cxx,v 1.17 1999/03/17 17:00:14 garland Exp $

 ************************************************************************/
#include <stdio.h>  
#include <stdmix.h>
#include <mixio.h>
#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#else
#  include <getopt.h>
#endif
#define CLEANUP(x)  if(x) { delete x; x=NULL; }

static char **global_argv;
static char *options = "O:B:W:t:Fo:m_pQslimMesh:c:rjI:M:qh";

static char *usage_string =
"-O <n>         Optimal placement policy:\n"
"                       0=endpoints, 1=endormid, 2=line, 3=optimal [default]\n"
"-B <weight>    Use boundary preservation planes with given weight.\n"
"-W <n>         Quadric weighting policy:\n"
"                       0=uniform, 1=area [default], 2=angle\n"
"-t <n>         Set the target number of faces.\n"
"-F             Use face contraction instead of edge contraction.\n"
"-o <file>      Output final model to the given file.\n"
"-I <file>      Deferred file inclusion.\n"
"-m_pQslimMesh <penalty>   Set the penalty for bad meshes.\n"
"-c <ratio>     Set the desired compactness ratio.\n"
"-r             Enable history recording.\n"
"-M <format>    Select output format:\n"
"                       {smf, iv, vrml, pm, mmf, log}\n"
"-q		Be quiet.\n"
"-j             Join only; do not remove any faces.\n"
"-h             Print help.\n"
"\n";




static
bool qslim_smf_hook(char *op, int, char *argv[], MxStdModel& pQslimMesh)
{
    if( streq(op, "e") )
    {
	if( !CMeshSimpDoc::target_edges )
	    CMeshSimpDoc::target_edges = new MxDynBlock<MxEdge>(pQslimMesh.vert_count() * 3);

	MxEdge& e = CMeshSimpDoc::target_edges->add();

	e.v1 = atoi(argv[0]) - 1;
	e.v2 = atoi(argv[1]) - 1;

	return true;
    }

    return false;
}

//bool (*unparsed_hook)(char *, int, char*[], MxStdModel&) = qslim_smf_hook;

//void slim_print_banner(ostream& out)
//{
//    out << "QSlim surface simplification software." << endl
//	<< "Version " << slim_version_string << " "
//	<< "[Built " << __DATE__ << "]." << endl
//	<< slim_copyright_notice << endl;
//}

void CMeshSimpDoc::slim_init()
{
	bool (*unparsed_hook)(char *, int, char*[], MxStdModel&) = qslim_smf_hook;
    if( !slim )
    {
	if( will_use_fslim )
	    slim = fslim = new MxFaceQSlim(*m_pQslimMesh);
	else
	    slim = eslim = new MxEdgeQSlim(*m_pQslimMesh);
    }
    else
    {
	if( will_use_fslim )
	    fslim = (MxFaceQSlim *)slim;
	else
	    eslim = (MxEdgeQSlim *)slim;
    }

    slim->placement_policy = placement_policy;
    slim->boundary_weight = boundary_weight;
    slim->weighting_policy = weighting_policy;
    slim->compactness_ratio = compactness_ratio;
    slim->meshing_penalty = meshing_penalty;
    slim->will_join_only = will_join_only;

    if( eslim && target_edges )
    {
	eslim->initialize(*target_edges, target_edges->length());
    }
    else
	slim->initialize();

    if( will_record_history )
    {
	if( !eslim )
	    mxmsg_signal(MXMSG_WARN,
			 "History only available for edge contractions.");
	else
	{
	    history = new QSlimLog(100);
	    //eslim->contraction_callback = slim_history_callback;
	}
    }
}


void CMeshSimpDoc::slim_cleanup()
{
    CLEANUP(smf);
    CLEANUP(m_pQslimMesh);
	CLEANUP(m_pOriginalQMesh);
    CLEANUP(slim);
    eslim = NULL;
    fslim = NULL;
    CLEANUP(history);
    CLEANUP(target_edges);
    if( output_stream != &cout )
    	CLEANUP(output_stream);
}

void CMeshSimpDoc::setup_output()
{
    if( !output_stream )
    {
	if( output_filename )
	    output_stream = new ofstream(output_filename);
	else
	    output_stream = &cout;
    }
}

bool CMeshSimpDoc::select_output_format(const char *fmt)
{
    bool h = false;

    if     ( streq(fmt, "mmf") ) { output_format = MMF; h = true; }
    else if( streq(fmt, "pm")  ) { output_format = PM;  h = true; }
    else if( streq(fmt, "log") ) { output_format = LOG; h = true; }
    else if( streq(fmt, "smf") ) output_format = SMF;
    else if( streq(fmt, "iv")  ) output_format = IV;
    else if( streq(fmt, "vrml")) output_format = VRML;
    else return false;

    if( h ) will_record_history = true;

    return true;
}

void CMeshSimpDoc::output_preamble()
{
    if( output_format==MMF || output_format==LOG )
	output_current_model();
}

void CMeshSimpDoc::output_current_model()
{
    setup_output();

    MxSMFWriter writer;
    writer.write(*output_stream, *m_pQslimMesh);

	CLEANUP(output_stream);
}

void CMeshSimpDoc::cleanup_for_output()
{
    // First, mark stray vertices for removal
    //
    for(uint i=0; i<m_pQslimMesh->vert_count(); i++)
	{
		if( m_pQslimMesh->vertex_is_valid(i) && m_pQslimMesh->neighbors(i).length() == 0 )
			m_pQslimMesh->vertex_mark_invalid(i);
	}

	// Compact vertex array so only valid vertices remain
    m_pQslimMesh->compact_vertices();
}

void CMeshSimpDoc::output_final_model()
{
    setup_output();

    switch( output_format )
    {
    case MMF:
	output_regressive_mmf(*output_stream);
	break;

    case LOG:
	output_regressive_log(*output_stream);
	break;

    case PM:
	output_progressive_pm(*output_stream);
	break;

    case IV:
	cleanup_for_output();
	output_iv(*output_stream);
	break;

    case VRML:
	cleanup_for_output();
	output_vrml(*output_stream);
	break;

    case SMF:
		// 暂时删除 [3/31/2012 Han]
	//cleanup_for_output();
	output_current_model();
	break;
    }


}

void CMeshSimpDoc::input_file(const char *filename)
{
    if( streq(filename, "-") )
	smf->read(cin, m_pQslimMesh);
    else
    {
	ifstream in(filename);
	if( !in.good() )
	    mxmsg_signal(MXMSG_FATAL, "Failed to open input file", filename);
	smf->read(in, m_pQslimMesh);
	in.close();
    }
}

static
MxDynBlock<char*> files_to_include(2);

void CMeshSimpDoc::defer_file_inclusion(char *filename)
{
    files_to_include.add(filename);
}

void CMeshSimpDoc::include_deferred_files()
{
    for(uint i=0; i<files_to_include.length(); i++)
	input_file(files_to_include[i]);
}

void CMeshSimpDoc::slim_history_callback(const MxPairContraction& conx, float cost)
{
    history->add(conx);
}

/************************************************************************

  Specialized QSlim output routines.  The bulk of the code here is to
  support multiresolution output formats.

  Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.
  
  $Id: output.cxx,v 1.4 1999/10/18 15:55:10 garland Exp $

 ************************************************************************/

void CMeshSimpDoc::output_ivrml(ostream& out, bool vrml)
{
    uint i;

    if( vrml )
	out << "#VRML V1.0 ascii" << endl;
    else
	out << "#Inventor V2.0 ascii" << endl;

    out << "Separator {" << endl
         << "Coordinate3 {" << endl
         << "point [" << endl;

    for(i=0; i<m_pQslimMesh->vert_count(); i++)
	if( m_pQslimMesh->vertex_is_valid(i) )
	    out << "   " << m_pQslimMesh->vertex(i)[0] << " "
		<< m_pQslimMesh->vertex(i)[1] << " "
		<< m_pQslimMesh->vertex(i)[2] << "," << endl;

    out << "]"<< endl << "}" << endl;
    out << "IndexedFaceSet {" << endl
         << "coordIndex [" << endl;

    for(i=0; i<m_pQslimMesh->face_count(); i++)
	if( m_pQslimMesh->face_is_valid(i) )
	    out << "   "
		<< m_pQslimMesh->face(i)[0] << ", "
		<< m_pQslimMesh->face(i)[1] << ", "
		<< m_pQslimMesh->face(i)[2] << ", "
		<< "-1," << endl;

    out << "]}}" << endl;
}

void CMeshSimpDoc::output_iv(ostream& out) { output_ivrml(out, false); }
void CMeshSimpDoc::output_vrml(ostream& out) { output_ivrml(out, true); }

void CMeshSimpDoc::output_regressive_mmf(ostream& out)
{
    if( !history ) return;

    out << "set delta_encoding 1" << endl;

    for(uint i=0; i<history->length(); i++)
    {
        MxPairContraction& conx = (*history)[i];

	// Output the basic contraction record
        out << "v% " << conx.v1+1 << " " << conx.v2+1 << " "
            << conx.dv1[X] << " " << conx.dv1[Y] << " " << conx.dv1[Z]
            << endl;

        // Output the faces that are being removed
        for(uint j=0; j<conx.dead_faces.length(); j++)
            out << "f- " << conx.dead_faces(j)+1 << endl;
    }
}

void CMeshSimpDoc::output_regressive_log(ostream& out)
{
    if( !history ) return;

    for(uint i=0; i<history->length(); i++)
    {
        MxPairContraction& conx = (*history)[i];

	// Output the basic contraction record
        out << "v% " << conx.v1+1 << " " << conx.v2+1 << " "
            << conx.dv1[X] << " " << conx.dv1[Y] << " " << conx.dv1[Z];

        // Output the faces that are being removed
        for(uint j=0; j<conx.dead_faces.length(); j++)
            out << " " << conx.dead_faces(j)+1;

        // Output the faces that are being reshaped
        out << " &";
        for(uint k=0; k<conx.delta_faces.length(); k++)
            out << " " << conx.delta_faces(k)+1;

        out << endl;
    }
}

void CMeshSimpDoc::output_progressive_pm(ostream& out)
{
    if( !history ) return;

    MxBlock<MxVertexID> vmap(m_pQslimMesh->vert_count());  // Maps old VIDs to new VIDs
    MxBlock<MxFaceID> fmap(m_pQslimMesh->face_count());    // Maps old FIDs to new FIDs
    uint i,k;

    MxVertexID next_vert = 0;
    MxFaceID   next_face = 0;

    ////////////////////////////////////////////////////////////////////////
    //
    // Output base mesh
    //
    for(i=0; i<m_pQslimMesh->vert_count(); i++)
	if( m_pQslimMesh->vertex_is_valid(i) )
	{
	    vmap(i) = next_vert++;
	    out << m_pQslimMesh->vertex(i) << endl;
	}
    
    for(i=0; i<m_pQslimMesh->face_count(); i++)
	if( m_pQslimMesh->face_is_valid(i) )
	{
	    fmap(i) = next_face++;
	    VID v1 = vmap(m_pQslimMesh->face(i)(0));
	    VID v2 = vmap(m_pQslimMesh->face(i)(1));
	    VID v3 = vmap(m_pQslimMesh->face(i)(2));

	    out << "f " << v1+1 << " " << v2+1 << " " << v3+1 << endl;
	}

    ////////////////////////////////////////////////////////////////////////
    //
    // Output mesh expansion
    //
    for(i=history->length()-1; i<=history->length(); i--)
    {	
	const MxPairContraction& conx = (*history)[i];
	SanityCheck( m_pQslimMesh->vertex_is_valid(conx.v1) );
	SanityCheck( !m_pQslimMesh->vertex_is_valid(conx.v2) );

	out << "v* " << vmap(conx.v1) + 1;
	out << "  "<<conx.dv1[X]<<" "<<conx.dv1[Y]<<" "<<conx.dv1[Z];
	out << "  "<<conx.dv2[X]<<" "<<conx.dv2[Y]<<" "<<conx.dv2[Z];
	out << " ";

	// Output new faces
	for(k=0; k<conx.dead_faces.length(); k++)
	{
	    FID fk = conx.dead_faces(k);
	    VID vk = m_pQslimMesh->face(fk).opposite_vertex(conx.v1, conx.v2);
 	    SanityCheck( m_pQslimMesh->vertex_is_valid(vk) );

 	    out << " ";
	    if( !m_pQslimMesh->face(fk).is_inorder(vk, conx.v1) ) out << "-";
 	    out << vmap(vk)+1;

	    fmap(conx.dead_faces(k)) = next_face++;
	    m_pQslimMesh->face_mark_valid(conx.dead_faces(k));
	}

	// Output delta faces
	out << " &";
	for(k=0; k<conx.delta_faces.length(); k++)
	{
	    out << " ";
	    FID fk = conx.delta_faces(k);
	    assert(m_pQslimMesh->face_is_valid(fk));
	    out << " " << fmap(fk)+1;
	}

	vmap(conx.v2) = next_vert++;
	m_pQslimMesh->vertex_mark_valid(conx.v2);
	out << endl;
    }
}


/************************************************************************

  QSlim command line program.  This provides a very simple interface to
  the underlying functionality.  Basically, it just reads in the input,
  simplifies it, and writes out the results.  It couldn't be simpler.

  Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.
  
  $Id: qslim.cxx,v 1.10 2000/11/20 20:52:41 garland Exp $

 ************************************************************************/

static ostream& vfcount(ostream& out, uint v, uint f)
{
    return out << "(" << v << "v/" << f << "f)";
}

//void CMeshSimpDoc::startup_and_input(int argc, char **argv)
//{
//    smf = new MxSMFReader;
//
//    process_cmdline(argc, argv);
//    if( m_pQslimMesh->face_count() == 0 )
//    {
//	smf->read(cin, m_pQslimMesh);
//    }
//
//    output_preamble();
//}


void CMeshSimpDoc::startup_and_input(char* filename, int otherType)
{
	smf = new MxSMFReader;

	int opt, ival;

	//global_argv = argv;

	smf->unparsed_hook = unparsed_hook;
	m_pQslimMesh = new MxStdModel(100, 100);

	if (otherType == 1)
	{
		ifstream in(filename);
		if( !in.good() )
			mxmsg_signal(MXMSG_FATAL, "Failed to open input file", filename);
		smf->read_off(in, m_pQslimMesh);
		in.close();
	}
	else// smf type
		input_file(filename);

	
	if( m_pQslimMesh->face_count() == 0 )
	{
		smf->read(cin, m_pQslimMesh);
	}
	if( m_pQslimMesh->normal_binding() == MX_UNBOUND )
	{
		m_pQslimMesh->normal_binding(MX_PERVERTEX/*MX_PERFACE*/);		// 逐顶点的normal计算
		m_pQslimMesh->synthesize_normals();
	}


	//output_preamble();
	//  [5/22/2012 Han]
	GetMiniBall(m_pQslimMesh);
	align_center(m_pQslimMesh);
}

void CMeshSimpDoc::GetMiniBall(MxStdModel *m)
{
	CLEANUP(m_miniBall);
	m_miniBall = new Miniball_C(3);


	// generate random points and check them in
	// ----------------------------------------
	Point p(3);
	for (int i=0; i<m->vert_count(); ++i) {
		for (int j=0; j<3; ++j)
			p[j] = m->vertex(i)[j];
		m_miniBall->check_in(p);
	}

	m_miniBall->build();
}

// 将所有顶点进行偏移，以便使模型中心和原点对齐
void CMeshSimpDoc::align_center(MxStdModel *m)
{
	for (unsigned int i = 0; i < m->vert_count(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			m->vertex(i)[j] -= m_miniBall->center()[j];
		}
	}
}
bool CMeshSimpDoc::UpdateCurModel()
{
	switch (m_mType)
	{
	case NORMAL:
		return false;	// not convert test
		break;
	case QSLIM:
		char fn[128];
		sprintf(fn, "QS%d.ply", rand());
		output_filename = fn;
		//output_filename = "QS.ply";
		output_final_model();
		//if (m_viewMyRender->load_model(std::wstring(output_filename, output_filename+6)))
		//{
		//	m_viewMyRender->initializeGL();
		//	m_mType = VIEW_SELECT;
		//	return true;
		//}
		
		if (LoadViewSelectModel(output_filename))
		{
			m_mType = VIEW_SELECT_OUR;
			return true;
		}

		break;
	default:
		return true;
		break;
	}
	
	return false;
}


void CMeshSimpDoc::SaveSMF2PLY(CArchive& ar, MxStdModel* m)
{
	unsigned int i;
	int * oldIndex2NewIndex = new int[m->vert_count()];
	memset(oldIndex2NewIndex, -1, m->vert_count()*sizeof(int));
	MxVertex *newIndexedVertex = new MxVertex[m->vert_count()];
	int nSimpedVN = 0;
	int nSimpedFN = 0;
	for (i=0;i <  m->face_count(); i++)
	{
		if (!m->face_is_valid(i))
			continue;
		nSimpedFN++;
		const MxFace& f = m->face(i);
		for(int j = 0; j < 3; j++)
		{
			if (oldIndex2NewIndex[f[j]] == -1)
			{
				oldIndex2NewIndex[f[j]] = nSimpedVN;
				newIndexedVertex[nSimpedVN] = m->vertex(f[j]);
				nSimpedVN++;
			}
		}
	}

	// header
	CString str;
	str.Format(
		"ply\r\nformat ascii 1.0\r\nelement vertex %d\r\nproperty float x\r\nproperty float y\r\nproperty float z\r\nelement face %d\r\nproperty list uchar int vertex_indices\r\nend_header\r\n"
		,nSimpedVN,nSimpedFN);
	ar.WriteString(str);

	// vetex
	for(unsigned int v=0; v<nSimpedVN; v++)
	{
		str.Format("%f %f %f\r\n", newIndexedVertex[v][0],newIndexedVertex[v][1],newIndexedVertex[v][2]);
		ar.WriteString(str);
	}

	// face
	for (i=0;i <  m->face_count(); i++)
	{
		if (!m->face_is_valid(i))
			continue;
		const MxFace& f = m->face(i);

		str.Format("3 %d %d %d\r\n", oldIndex2NewIndex[f[0]] ,oldIndex2NewIndex[f[1]] ,oldIndex2NewIndex[f[2]]);
		ar.WriteString(str);
	}

	delete[] oldIndex2NewIndex;
	delete[] newIndexedVertex;
		// 不保存其他信息 [5/3/2012 Han]
		//if( m.normal_binding() != MX_UNBOUND )
		//{
		//	ar << "bind n " << m.binding_name(m.normal_binding()) << "\r\n";
		//	for(i=0; i<m.normal_count(); i++)
		//		ar << "n " << m.normal(i)[0] << " "<<m.normal(i)[1] << " " << m.normal(i)[2] << "\r\n";
		//}

		//if( m.color_binding() != MX_UNBOUND )
		//{
		//	ar << "bind c " << m.binding_name(m.color_binding()) << "\r\n";
		//	for(i=0; i<m.color_count(); i++)
		//		ar<< "c " << m.color(i).R() << " " << m.color(i).G() << " " << m.color(i).B() << "\r\n" ;
		//}

		//if( m.texcoord_binding() != MX_UNBOUND )
		//{
		//	ar << "tex " << m.texmap_name() << "\r\n";
		//	ar << "bind r " << m.binding_name(m.texcoord_binding()) << "\r\n";
		//	for(i=0; i<m.texcoord_count(); i++)
		//		ar << "r " << m.texcoord(i)[0] << " " << m.texcoord(i)[1] << "\r\n";
		//}
}
//
//int main(int argc, char **argv)
//{
//    double input_time, init_time, slim_time, output_time;
//
//    // Process command line and read input model(s)
//    //
//    TIMING(input_time, startup_and_input(argc, argv));
//
//    if(!be_quiet) cerr << "+ Initial model    ";
//    if(!be_quiet) vfcount(cerr, m_pQslimMesh->vert_count(), m_pQslimMesh->face_count()) << endl;
//
//    // Initial simplification process.  Collect contractions and build heap.
//    //
//    TIMING(init_time, slim_init());
//
//    // Decimate model until target is reached
//    //
//    TIMING(slim_time, slim->decimate(face_target));
//
//    if(!be_quiet) cerr << "+ Simplified model ";
//    if(!be_quiet) vfcount(cerr, slim->valid_verts, slim->valid_faces) << endl;
//
//    // Output the result
//    //
//    TIMING(output_time, output_final_model());
//
//    if( !be_quiet )
//    {
//	cerr << endl << endl;
//	cerr << "+ Running time" << endl;
//	cerr << "    Setup      : " << input_time << " sec" << endl;
//	cerr << "    QSlim init : " << init_time << " sec" << endl;
//	cerr << "    QSlim run  : " << slim_time << " sec" << endl;
//	cerr << "    Output     : " << output_time << " sec" << endl;
//	cerr << endl;
//	cerr << "    Total      : "
//	     << input_time+init_time+slim_time+output_time <<endl;
//    }
//    else
//    {
//	cerr << slim->valid_faces << " " << init_time+slim_time << endl;
//    }
//
//    slim_cleanup();
//
//    return 0;
//}
